#!/usr/bin/perl


#eval 'exec perl -S $0 ${1+"$@"}'
#if 0;

use lib "$ENV{HOME}/lib/perl";

#require "hilbert_su.pl"; # Hilbert Transform
#require "my_math.pl";
#require "hilbert.pl"; # Hilbert Envelope
use Math::FFT;
use POSIX;
#
use strict 'vars';
use warnings;
use vars qw/%opt/;

my $progname="CSpolar.pl";

my $debug=0;
my $rpd= 0.01745329251994;
my $maxbaz=360; # search 0..maxbaz

# Parse commandline and read parmfile
&init();
if ($opt{v} ) {
	$debug=1;
	print STDOUT "Program set to grossly spew debugs\n";
}
if ($opt{d} ) {
	$debug=2;
	print STDOUT "Program set to grossly spew debugs and extra files\n";
}
my ($infile);
if (!$opt{f} ) {
	print STDOUT "ERROR: Need to specify an input file \n\n" ;
	&usage();
} elsif ($opt{f}) {
	$infile=$opt{f};
	die "$!: $infile\n" if (! -e $infile);
}
my ($outfile);
if (!$opt{o} ) {
	print STDERR "ERROR: Need to specify an outfile\n\n" ;
	&usage();
} elsif ($opt{o}) {
	$outfile=$opt{o};
}
my $writeout=0;
if ($opt{w} ) {
	$writeout=1;
}

#----------------------------------------------
my (@T,@Z,@N,@E);
my (@oT);
my ($nsamp);
my $Czr_max=-1e99;
my $Czr2_max=-1e99;
my ($Czr_az,$Czr2_az);
#----------------------------------------------
# Open data file and parse
open(IF,"<$infile") or die "$!: $infile\n";
my $n=0;
while (<IF>) {
	chomp;
	s/^\s+//;
	my $n=$. - 1;
	($T[$n],$Z[$n],$N[$n],$E[$n])=split(/\s+/);
}
close(IF);
$nsamp=scalar(@T);
my $samprate=1/($T[1]-$T[0]);

if ($debug>0) {
	print STDERR "Analysis window length: $nsamp samples \n";
}

my (@dz,@dn,@de);
@dz=@Z;
@dn=@N;
@de=@E;
my $ndz=scalar(@dz);
my $ndn=scalar(@dn);
my $nde=scalar(@de);
my @tt=@T;
my $myt=0;  # in here for legacy calling scripts TODO

#----------------------------------------------
my @hz=&hilbert_su(@dz); # do the hilbert dance, Hilbert Transform from SU
my $pow2=&nextpow2($ndz); # for zero padding
my $nsampfft=2*(2**$pow2);
for (my $n=$ndz; $n < $nsampfft; $n++) {
	push(@hz,0.0); # pad data vector
}
open(O,">$outfile") or die "$!: $outfile\n";
if ($writeout>0) {
	open(HZ,">hilbz.txt") or die "$!: hilb\n";
	for (my $n=0; $n < $nsampfft; $n++) {
		printf HZ ("%f\n", $hz[$n]);
	}
	close(HZ);
}
my $fftz = new Math::FFT(\@hz);
my $zzcorr = $fftz->correl($fftz); # Z-Z xcorr
my @zzxcor=@$zzcorr;
my $nzz=scalar(@zzxcor);
my ($zzmaxamp,$zzminamp,$zzmean,$zzabsmax) = &find_mmm(@zzxcor);
# Below is the zero-lag amplitude of the cross correlation
# correl spits out the correlation 
# +lag [0 .. nsampfft/2] then -lag [nsampfft/2 .. nsampfft]
# i think
my $zz0lag_amp=$zzxcor[0];
my @zzxcor_out=();
print "ZZXCOR $nzz $nsampfft -- $zzmaxamp,$zzminamp,$zzmean,$zzabsmax -\n" if $debug>0;
# HACK This is not a good way to reorder zzxcor
# Checked in octave, zzxcor_out is the same as xcorr(hz')
for (my $n=$nsampfft/2; $n < $nsampfft; $n++) {
	push(@zzxcor_out, $zzxcor[$n]);
}
for (my $n=0; $n < $nsampfft/2; $n++) {
	push(@zzxcor_out, $zzxcor[$n]);
}
if ($debug>1) {
	open(ZZ,">zzxcor.txt") or die "$!: zzcor\n";
	for (my $n=$nsampfft/2; $n < $nsampfft; $n++) {
		printf ZZ ("%f\n", $zzxcor[$n]);
	}
	for (my $n=0; $n < $nsampfft/2; $n++) {
		printf ZZ ("%f\n", $zzxcor[$n]);
	}
	close(ZZ);
}

my @baz=(0..$maxbaz);
my $ndeg=scalar(@baz);
if ($writeout>0) {
	open(R2,">radial_2.txt") or die "$!: radial.txt\n";
}
if ($debug>1) {
	open(C,">corr.txt") or die "$!: corr.txt\n";
	open(R,">radial.txt") or die "$!: radial.txt\n";
}
print STDOUT "Looping over back-azimuth ...\n" if $debug>0;
for (my $j=0; $j<$ndeg; $j++) { # Loop over backaz to compute radial comp.
	my $phi=$baz[$j];
	my ($drref,$dtref)=&rotate_horiz($phi,\@de,\@dn);
	my @data_radial=@$drref;
	my @data_tang=@$dtref;

	for (my $n=$ndz; $n < $nsampfft; $n++) {
		push(@data_radial,0.0); # pad data vector
	}
	if ($debug>1) {
		for (my $n=0; $n < $nsampfft; $n++) {
			printf R ("%f ", $data_radial[$n]);
		}
		printf R "\n";
	}
	if ($writeout>0) {
		my @rad_env=&hilbert(@data_radial); # Hilbert envelope
		for (my $n=0; $n < $ndz; $n++) {
			#printf R2 ("%f %f %f \n", $tt[$n],$phi,$data_radial[$n]);
			printf R2 ("%f %f %f \n", $tt[$n],$phi,$rad_env[$n]);
		}
	}
	my $fftr = new Math::FFT(\@data_radial);
	my $zrcorr = $fftz->correl($fftr); # Z-R xcorr
	my @zrxcor=@$zrcorr;
	my ($zrmaxamp,$zrminamp,$zrmean,$zrabsmax) = &find_mmm(@zrxcor);
	my $zr0lag_amp=$zrxcor[0];
	my $ncor=scalar(@zrxcor);
	if ($debug>1) {
		#print "$ncor Correlation pts\n";
		#for (my $k=0; $k<$ncor; $k++) {
		#	printf C ("%f ", $zrxcor[$k]);
		#}
		#printf C "\n";
		printf C ">\n";
		for (my $k=0; $k<$ncor; $k++) {
			printf C ("%f %f %f\n", $k,$phi,$zrxcor[$k]);
		}
		printf C ">\n";
	}
	my $rrcorr = $fftr->correl($fftr); # R-R xcorr
	my @rrxcor=@$rrcorr;
	my ($rrmaxamp,$rrminamp,$rrmean,$rrabsmax) = &find_mmm(@rrxcor);
	my $rr0lag_amp=$rrxcor[0];
	my $Czr=$zr0lag_amp/sqrt($zz0lag_amp*$rr0lag_amp); #Szr/sqrt(Szz*Srr)
	my $Czr2=$zr0lag_amp/$zz0lag_amp; #Szr/Szz
	if ($Czr>$Czr_max) { # Positive max of Czr is computed azimuth
		$Czr_max=sprintf("%.5f",$Czr);
		$Czr_az=$phi;
	}
	if ($Czr2>$Czr2_max) {
		$Czr2_max=sprintf("%.5f",$Czr2);
		$Czr2_az=$phi;
	}
	# time backaz Szz Szr Srr 
	#printf O ("%f %f %f %f %f\n", $myt, $phi, $zzabsmax, $zrabsmax, $rrabsmax);
	printf O ("%f %f %f %f %f %f %f\n", $myt, $phi, $zz0lag_amp, $zr0lag_amp, $rr0lag_amp,$Czr,$Czr2);
} # baz loop
if ($debug>1) {
	close(C);
	close(R);
}
if ($writeout>0) {
	close(R2);
}
close(O);
print STDOUT "Czr_max: $Czr_max Czr_az: $Czr_az Czr2_max: $Czr2_max Czr2_az: $Czr2_az \n";

# SUBROUTINES -------------------------------------------------

sub init()
{
	use Getopt::Std;
	my $opt_string = 'hvdwf:o:';
	getopts( "$opt_string", \%opt ) or usage();
	usage() if $opt{h};
}

sub usage()
{
print STDOUT << "EOF";
NAME
\t$progname - Perform 3-comp polarization analysis on surface wave
\t            via the Chael/Selby method (see refs below).
\t

SYNOPSIS
\t$progname [-h|-v|-d|-w] [-p out.ps] -f datafile -o junk.out 

DESCRIPTION
\t 
\t
\t

OPTIONS
\t-f data file (REQUIRED). Expect 4-colum ascii (T,Z,N,E)
\t-o output file (REQUIRED). Output fields:
\t   time(not used) azimuth Szz Szr Srr Czr Czr2
\t-w write out files for plotting (OPTIONAL).
\t-h print this usage (OPTIONAL).
\t-v spew some messages (OPTIONAL).
\t-d debug, write lots of files (OPTIONAL).

REFERENCES
\tChael, E. (1997), An automated Rayleigh-Wave Detection Algorithm,
\t\tBSSA, 87.
\tSelby, N.D. (2001), Association of Rayleigh Waves using backazimuth
\t\tmeasurements: Application to test ban verification, BSSA, 91.
\tBaker and Stevens (2004), Backazimuth estimation reliability using
\t\tsurface wave polarization, GRL, 31.
\tStachnik et al. (2012), Determination of New Zealand Ocean Bottom 
\t\tSeismometer Orientation via Rayleigh-wave Polarization, SRL, 83.

EOF
die "\n";
}

sub isint {
    my $val = shift;
    # "not an integer"       unless /^-?\d+$/;
    # "not an integer"       unless /^[+-]?\d+$/;
    return ($val =~ m/^\d+$/);
}

sub rotate_horiz {
	my $phi=shift;
	my @datae=@{ $_[0] };
	my @datan=@{ $_[1] };
	my $npts=scalar(@datae);
	my $mpi=3.14159265358979;
	my $rpd=$mpi/180;
	my $sinphi=sin($phi*$rpd);
	my $cosphi=cos($phi*$rpd);
	my @data_east=();
	my @data_north=();
	for (my $i=0; $i<$npts; $i++) {
		my $temp = $datan[$i]*$cosphi + $datae[$i]*$sinphi;
		$data_north[$i] = $temp;           #/* Radial */
		$data_east[$i] = $datae[$i]*$cosphi - $datan[$i]*$sinphi; #/* Changed 2nd term sign from RM's version */
		$data_east[$i] = $data_east[$i]*(-1.0);    # /* Tangential */
	}
	my @data_radial=@data_north;
	my @data_tang=@data_east;
	my @out=();
	$out[0]=[ @data_radial ];
	$out[1]=[ @data_tang ];

	return(\@data_radial,\@data_tang);
}
sub nextpow2 {
    my $n=shift(@_);
    my $i;
    for ($i=1; 2**$i < $n; $i++){};
    return ($i+1);
}

sub find_mmm {
    # my ($maxamp,$minamp,$mean,$absmax) = find_mmm(@dd);
    my @d=@_;
    my $i;
    my $sum=0;
    my $maxamp=-1e99;
    my $minamp=1e99;
    my $mean=0;
    my $absmax;
    for ($i=0; $i< $#d; $i++) {
        $maxamp=$d[$i] if ($d[$i] > $maxamp);
        $minamp=$d[$i] if ($d[$i] < $minamp);
        $sum+=$d[$i];
    }
    $mean=$sum/$#d;

    $mean=sprintf("%0.8f",$mean);
    $maxamp=sprintf("%0.8f",$maxamp);
    $minamp=sprintf("%0.8f",$minamp);
    if (abs($minamp) > $maxamp) {
        $absmax = abs($minamp);
    } else {
        $absmax = $maxamp;
    }
    return ($maxamp,$minamp,$mean,$absmax);
}

sub find_sum {
	my @d=@_;
	my $sum=0;
	my $i;
    for ($i=0; $i< $#d; $i++) {
        $sum+=$d[$i];
    }
	return $sum ;
}
sub hilbert_su {
	# these algorithms taken from CWP (seismic unix)
#/*****************************************************************************
#HILBERT - Compute Hilbert transform y of x
#
#hilbert     compute the Hilbert transform
#
#******************************************************************************
#Function Prototype:
#void hilbert (int n, float x[], float y[]);
#
#******************************************************************************
#Input:
#n       length of x and y
#x       array[n] to be Hilbert transformed
#
#Output:
#y       array[n] containing Hilbert transform of x
#
#******************************************************************************
#Notes:
#The Hilbert transform is computed by convolving x with a
#windowed (approximate) version of the ideal Hilbert transformer.
#    
#******************************************************************************
#Author:  Dave Hale, Colorado School of Mines, 06/02/89
#*****************************************************************************/

	use strict;
	use warnings;
	my $debug=0;
	my @x=@_;
	my $n=scalar(@x);
	print "input signal has $n samples\n" if ($debug > 0);
	my $LHHALF=50;   # half-length of Hilbert transform filter
					 # make smaller to run faster. CWP had set to 30
					 # I set it to 50 to more closely match matlab output
	my $LH=2*$LHHALF+1;   # filter length must be odd 
	my $PI=3.141592653589793;
	my @h; # size LH
	my $i;
	my $taper;
	# Hilbert transform filter; use Hamming window 
	$h[$LHHALF]=0.0;
	for ($i=1; $i<=$LHHALF; $i++) {
		$taper = 0.54+0.46*cos($PI*$i/$LHHALF);
		$h[$LHHALF+$i] = $taper*((-1*($i%2)*2.0)/($PI*$i));
		$h[$LHHALF-$i] = -1*$h[$LHHALF+$i];
	}
	# convolve Hilbert transform with input array 
	my $nh=scalar(@h);
	print "LH=$LH LHHALF=$LHHALF nh=$nh n=$n\n" if ($debug > 0);
	if ($debug > 1) {
		open(OF,">hamm.x");
		foreach my $i (@h) {
			print OF "$i\n";
		}
		close(OF);
	}
    my @y=&conv($LH,-1*$LHHALF,@h,$n,0,@x,$n,0);
	return @y;
}
sub conv {
#/*****************************************************************************
#CONVOLUTION - Compute z = x convolved with y
#
#conv    compute the convolution of two input vector arrays
#
#******************************************************************************
#Input:
#lx      length of x array
#ifx     sample index of first x
#x       array[lx] to be convolved with y
#ly      length of y array
#ify     sample index of first y
#y       array[ly] with which x is to be convolved
#lz      length of z array
#ifz     sample index of first z
#
#Output:
#z       array[lz] containing x convolved with y
#
#******************************************************************************
#Function Prototype:
#void conv (int lx, int ifx, float *x, int ly, int ify, float *y,
#    int lz, int ifz, float *z);
#
#******************************************************************************
#Notes:
#The operation z = x convolved with y is defined to be
#           ifx+lx-1
#    z[i] =   sum    x[j]*y[i-j]  ;  i = ifz,...,ifz+lz-1
#            j=ifx
#The x samples are contained in x[0], x[1], ..., x[lx-1]; likewise for
#the y and z samples.  The sample indices of the first x, y, and z values
#determine the location of the origin for each array.  For example, if
#z is to be a weighted average of the nearest 5 samples of y, one might
#use 
#    ...
#    x[0] = x[1] = x[2] = x[3] = x[4] = 1.0/5.0;
#    conv(5,-2,x,lx,0,y,ly,0,z);
#    ...
#In this example, the filter x is symmetric, with index of first sample = -2.
#
#This function is optimized for architectures that can simultaneously perform
#a multiply, add, and one load from memory; e.g., the IBM RISC System/6000.
#Because, for each value of i, it accumulates the convolution sum z[i] in a
#scalar, this function is not likely to be optimal for vector architectures.

#******************************************************************************
#Author:  Dave Hale, Colorado School of Mines, 11/23/91
#*****************************************************************************/

	my $lx=shift(@_); # LH
	my $ifx=shift(@_); # -LHHALF
	my @x=splice(@_,0,$lx); # transform filter
	my $ly=shift(@_); # size of signal to be filtered
	my $ify=shift(@_); # this should be zero
	my @y=splice(@_,0,$ly);
	my $lz=shift(@_); # size of signal to be filtered
	my $ifz=shift(@_); # should be zero
	my @z; # transformed output
	my $debug=0;
	my $nx=scalar(@x);
	my $ny=scalar(@y); 
	print "lx=$lx ifx=$ifx nx=$nx ly=$ly ify=$ify ny=$ny lz=$lz ifz=$ifz\n" if ($debug > 0);
	my $ilx=$ifx+$lx-1;
	my $ily=$ify+$ly-1;
	my $ilz=$ifz+$lz-1;
	my ($i,$j,$jlow,$jhigh);
	my $sum;

# this line in c resets pointers (array index)
##	x -= ifx;  y -= ify;  z -= ifz;
	for ($i=$ifz; $i<=$ilz; ++$i) {
		$jlow = $i-$ily;  
		$jlow = $ifx if ($jlow<$ifx);
		$jhigh = $i-$ify;
 		$jhigh = $ilx if ($jhigh>$ilx);
		$sum=0.0;
		for ($j=$jlow; $j<=$jhigh; ++$j) {
			#$sum += $x[$j]*$y[$i-$j];
			$sum += $x[$ifx - $j - 1 ]*$y[$i-$j]; # the $ifx-1 adjusts index
		}
		$z[$i] = $sum;
	}

	return(@z);

}

sub hilbert {
	# usage: @z=%hilbert(@dd);
	use strict;
	use warnings;
	use Math::FFT;
	my @dd=@_;
	my $nsamp=scalar(@dd);
	my $length2; # Next pow2
	my $nfft; # 2*$length2
	my @dd1; # zero-padded, imag. data array
	my ($m); # loop-increment 

	# Find nextpow2 
    for($m=1;2**$m <= $nsamp;$m++) {
        $length2=2**($m+1);
    }
	# Double the length, zero-pad, make array Imag.
	$nfft=2*$length2;
	my $nyq=$nfft/2;
	
	for ($m=0; $m < $nfft; $m++) {$dd1[$m]=0.0;}
    for ($m=0; $m < $nsamp; $m++) {$dd1[$m*2]=$dd[$m];}

    my $ddi=new Math::FFT(\@dd1); #initialize FFT 
    my $hb=$ddi->cdft($ddi);
	# Do hilbert
	for($m= 0; $m < $nfft;  $m++){
        if($m == 0) {
            $hb->[$m] *= 1.0;
        }  elsif($m == 1) {
            $hb->[$m] *= 0.0;
        } elsif($m == $nyq) {
            $hb->[$m] *= 1.0;
        } elsif($m < $nyq && $m > 0) {
            $hb->[$m] *= 2.0;
        } else {
            $hb->[$m] = 0.0;
        }
    }
    $ddi=$ddi->invcdft($hb);
	my (@tmp);
	for($m=0;$m<$nfft;$m+=2) {
		$tmp[$m]=sqrt(($ddi->[$m])**2+($ddi->[$m+1])**2);
	}
	my (@z);
    for($m=0;$m<2*$nsamp;$m+=2) {
            $z[$m/2]= $tmp[$m];
    }

	return @z;
}

