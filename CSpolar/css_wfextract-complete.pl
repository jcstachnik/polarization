#!/opt/antelope/5.0-64/bin/perl

eval 'exec perl -S $0 ${1+"$@"}'
if 0;

use lib "$ENV{ANTELOPE}/data/perl";
#use lib "$ENV{HOME}/lib/perl";
use Datascope;
use strict 'vars';
use warnings;
use vars qw/%opt/;
#require "costaper.pl";
#require "find_max.pl";
#require "my_math.pl";

my $progname="css_wfextract.pl";

# fields that are required in the parmfile
my @pfields=qw(dbname orid staname phase filter
				pretime posttime chans pertaper wffile normal);
my ($mypfile); 
my %pf; # hash to store pf key/values
my $debug=0;
my $correct_orient=0;
my $do_arrt=1;
my $opttime=0;
my $guesschan=0;
my $rotate=0;
my $rotate_angle=0;
my $baz=0;
my $inc=0;

# Parse commandline and read parmfile
&init();
if ($opt{v} ) {
	$debug=1;
	print STDOUT "Program set to grossly spew debugs\n";
}
if ($opt{o} ) {
	$correct_orient=1;
	print STDOUT "Program set to correct for hang/vang\n" if $debug==1;
}
if ($opt{g} ) {
	$guesschan=1;
	print STDOUT "Program set to guess Z,N,E based on input channel codes.\n" if $debug==1;
}
if ($opt{t} ) {
	$do_arrt=0;
	$opttime=$opt{t};
}
if ($opt{r} ) {
	$rotate=1;
	print STDOUT "Program set to rotate Z,N,E into Z,R,T based on SEAZ.\n" if $debug==1;
}
if ($opt{R} ) {
	$rotate=2;
	$rotate_angle=$opt{R};
	print STDOUT "Program set to rotate Z,N,E into Z,R,T based on $rotate_angle.\n" if $debug==1;
}
if ($opt{L} ) {
	$rotate=3;
	($baz,$inc) = split(/\,/,$opt{L});
	print STDOUT "Program set to rotate Z,N,E into L,Q,T based on baz=$baz inc=$inc.\n" if $debug==1;
}

if (!$opt{p} ) {
	print STDOUT "ERROR: Need to specify a parameter file \n\n" ;
	&usage();
} elsif ($opt{p}) {
	$mypfile=$opt{p};
	die "$!: $mypfile\n" if (! -e $mypfile);
	%pf=&read_parmfile($mypfile,@pfields);
}

if ($debug == 1) { # Check values read from parmfile
	my @keys=keys %pf;
	foreach my $k (@keys) {
		my $v=$pf{$k};
		print "PF: $k => $v\n";
	}
}
#----------------------------------------------
# DATABASE and TRACETABLE stuff
my @db=dbopen ($pf{dbname}, "r+");

# GET ARR TIME from arrival assoc origin join, or use time specified on cmd line
my ($arrt,$seaz);
if ($do_arrt==1) {
	$arrt=&get_arrt($pf{orid},$pf{staname},$pf{phase},@db);
	$seaz=&get_seaz($pf{orid},$pf{staname},$pf{phase},@db);
} else {
	$arrt=is_epoch_string($opttime);
	die "Input time $opttime is unrecognizable\n" if !defined $arrt ;
}

print "Using arrival time $arrt\n" if ($debug ==1);
my $t0=$arrt -  $pf{pretime}; 
my $t1=$arrt + $pf{posttime};
my $tmp=$t1-$t0;
print "Using time window $t0 to $t1, $tmp seconds\n";

# Load data into trace table for specified time window
# based on this subset
my $sta=$pf{staname};
my @chans=split(/\s+/,$pf{chans}); 

# set up filter to apply to trace data with trfilter
my $filt; # TODO set filter in pf i.e. DEMEAN; BW .3 4 8 4
$filt=$pf{filter};

my @data;
my ($t0got,$t1got,$nsamp,$samprate,$hang,$vang);
my %comps;
my ($chanN, $chanE, $chanZ);
my (@dNorth,@dEast,@dVert);
my ($maxN,$maxE,$maxV);
my $nc=1;
foreach my $ch (@chans) {
	my @tr=&get_tr($t0,$t1,$sta,$ch,@db);
	# Get actual time of data window -- assume this will be the
	# same for all channels 
	$tr[3]=0; # There should be 1 record in trace table--but force anyway
	($t0got,$t1got,$nsamp,$samprate)=
			dbgetv(@tr,qw(time endtime nsamp samprate ));
	trfilter(@tr,$filt); # $ret=trfilter(); if $ret==-2 then parsing error
	# Get hang/vang from sitechan does not exist in trace table
	my $chsub="sta=~/$sta/ && chan=~/$ch/ && ondate<=$t0got && (offdate>=$t1got || offdate==NULL)";
	my @dbs=dbprocess(@db,"dbopen sitechan",
						"dbsubset $chsub");
	my $ns=dbquery(@dbs,"dbRECORD_COUNT");
	if ($ns!=1) {
		print "Warning: $ns records after sitechan subset...\n";
		if ($correct_orient==1 || $rotate>0) {
			die("Quitting...because command line calls for horizontal component correction.\n")
		}
	}
	$dbs[3]=0;
	($hang,$vang)=dbgetv(@dbs,qw(hang vang));
	# - Process tr
	@data=trdata(@tr);
	@data=&costaper(@data,$pf{"pertaper"});
	trfree(@tr);
	$comps{$ch}{"samprate"}=$samprate;
	$comps{$ch}{"nsamp"}=$nsamp;
	$comps{$ch}{"time"}=$t0got;
	$comps{$ch}{"data"}=[ @data ];
	$comps{$ch}{"phslag"}=$t0got-$arrt;
	$comps{$ch}{"hang"}=$hang;
	$comps{$ch}{"vang"}=$vang;
	if ($hang == -999.9 || $vang == -999.9) {
		printf "Warning: hang ($hang) or vang ($vang) not defined in sitechan table for $ch !\n";
	}
	my @tt;
	for ( my $i=0; $i<$nsamp; $i++)  {
		$tt[$i]=$i/$samprate + $comps{$ch}{"phslag"};
	}
	$comps{$ch}{"tt"}=[ @tt ];
	if ($guesschan==1) {# Guess components
		my $ch3 = substr($ch,2,1); # orientation code, otherwise ch can match HNZ, etc.
		if ($ch3 =~/N/ || $ch3 =~ /R/ || $ch3=~/Y/ || $ch3=~/1/) {
			$chanN=$ch;
			@dNorth=@{$comps{$chanN}{data}};
			$comps{"datan"}=[ @dNorth ];
			$comps{"vangn"}=$vang;
			$comps{"hangn"}=$hang;
			print "North channel: $chanN nsamp $nsamp \n" if ($debug > 0);
		}
		if ($ch3 =~/E/ || $ch3 =~ /T/ || $ch3=~/X/ || $ch3=~/2/) {
			$chanE=$ch ;
			@dEast=@{$comps{$chanE}{data}};
			$nsamp=$comps{$chanE}{nsamp};
			$samprate=$comps{$chanE}{samprate};
			$comps{"datae"}=[ @dEast ];
			$comps{"vange"}=$vang;
			$comps{"hange"}=$hang;
			print "East channel: $chanE nsamp $nsamp \n" if ($debug > 0);
		}
		if ($ch3 =~/Z/ || $ch3=~/vert/) {
			$chanZ=$ch;
			@dVert=@{$comps{$chanZ}{data}};
			$comps{"dataz"}=[ @dVert ];
			$comps{"vangz"}=$vang;
			$comps{"hangz"}=$hang;
			$comps{"chanz"}=$chanZ;
			print "Vert channel: $chanZ nsamp $nsamp \n" if ($debug > 0);
		}
	} else { # Channels are ZNE in order specified in pf
		if ($nc==1) {
			$chanZ=$ch;
			@dVert=@{$comps{$chanZ}{data}};
			$comps{"dataz"}=[ @dVert ];
			$comps{"vangz"}=$vang;
			$comps{"hangz"}=$hang;
			$comps{"chanz"}=$chanZ;
			print "Vert channel: $chanZ nsamp $nsamp \n" if ($debug > 0);
		}
		if ($nc==2) {
			$chanN=$ch;
			@dNorth=@{$comps{$chanN}{data}};
			$comps{"datan"}=[ @dNorth ];
			$comps{"vangn"}=$vang;
			$comps{"hangn"}=$hang;
			print "North channel: $chanN nsamp $nsamp \n" if ($debug > 0);
		}
		if ($nc==3) {
			$chanE=$ch ;
			@dEast=@{$comps{$chanE}{data}};
			$nsamp=$comps{$chanE}{nsamp};
			$samprate=$comps{$chanE}{samprate};
			$comps{"datae"}=[ @dEast ];
			$comps{"vange"}=$vang;
			$comps{"hange"}=$hang;
			print "East channel: $chanE nsamp $nsamp \n" if ($debug > 0);
		}
	}
	$nc++;
} # End loop over channels
# Do some checks
if ( ! defined $chanN || ! defined $chanE || ! defined $chanZ ) {
	die "Big Problem: Missing one of the components\nchanN $chanN chanE $chanE chanZ $chanZ\n";
}
if ( $comps{$chanN}{nsamp} != $comps{$chanE}{nsamp} ) {
	my $tmp0=$comps{$chanN}{nsamp};
	my $tmp1=$comps{$chanE}{nsamp};
	die "Fatal Error: E($tmp1) and N($tmp0) nsamps differ!\n";
}
if ( $comps{$chanN}{nsamp} != $comps{$chanZ}{nsamp} ) {
	my $tmp0=$comps{$chanN}{nsamp};
	my $tmp1=$comps{$chanZ}{nsamp};
	die "Fatal Error: N($tmp0) and Z($tmp1) nsamps differ!\n";
}
if ( $comps{$chanE}{nsamp} != $comps{$chanZ}{nsamp} ) {
	my $tmp0=$comps{$chanE}{nsamp};
	my $tmp1=$comps{$chanZ}{nsamp};
	die "Fatal Error: E($tmp0) and Z($tmp1) nsamps differ!\n";
}

# Correct for hang/vang
if ($correct_orient == 1) {
	my @dummy=();
	$comps{"dz_true"}=[ @dummy ];
	$comps{"dn_true"}=[ @dummy ];
	$comps{"de_true"}=[ @dummy ];
	&correct_orientation(\%comps);
	$comps{$chanZ}{"data"}=[ @{$comps{"dz_true"}} ];
	$comps{$chanN}{"data"}=[ @{$comps{"dn_true"}} ];
	$comps{$chanE}{"data"}=[ @{$comps{"de_true"}} ];
	@dVert=@{$comps{$chanZ}{data}};
	@dNorth=@{$comps{$chanN}{data}};
	@dEast=@{$comps{$chanE}{data}};
}
if ($rotate>0) {
	my @rotdata=();
	my $rotang;
	if ($rotate==3) {
	# 3 comp rotation
		@rotdata=&rotate_ZNE_LQT($nsamp,$baz,$inc,@dVert,@dNorth,@dEast);
		$comps{$chanZ}{"data"}=[ @rotdata[0..$nsamp-1] ]; # L, longitudinal
		$comps{$chanN}{"data"}=[ @rotdata[0..$nsamp-1] ]; # Q, 
		$comps{$chanE}{"data"}=[ @rotdata[$nsamp .. $#rotdata] ]; # T, Tangential
		@dVert=@{$comps{$chanZ}{data}}; # L
		@dNorth=@{$comps{$chanN}{data}}; # Q
		@dEast=@{$comps{$chanE}{data}}; # T
	} elsif ($rotate==1) {
		$rotang=$seaz;
		@rotdata=&rot2d($nsamp,$rotang,@dNorth,@dEast);
		$comps{$chanN}{"data"}=[ @rotdata[0..$nsamp-1] ]; # Radial
		$comps{$chanE}{"data"}=[ @rotdata[$nsamp .. $#rotdata] ]; # Tangential
		@dNorth=@{$comps{$chanN}{data}};
		@dEast=@{$comps{$chanE}{data}};
	} elsif ($rotate==2) {
		$rotang=$rotate_angle;
		@rotdata=&rot2d($nsamp,$rotang,@dNorth,@dEast);
		$comps{$chanN}{"data"}=[ @rotdata[0..$nsamp-1] ]; # Radial
		$comps{$chanE}{"data"}=[ @rotdata[$nsamp .. $#rotdata] ]; # Tangential
		@dNorth=@{$comps{$chanN}{data}};
		@dEast=@{$comps{$chanE}{data}};
	}
}
# --- Data extracted now ampltude scale and print
my ($ampN,$ampE,$ampV);
$nsamp=$comps{$chanZ}{nsamp};
$maxV=&find_max(@dVert);
$maxN=&find_max(@dNorth);
$maxE=&find_max(@dEast);
print "Vert channel: $chanZ nsamp $nsamp maxamp $maxV\n" if ($debug > 0);
print "North channel: $chanN nsamp $nsamp maxamp $maxN\n" if ($debug > 0);
print "East channel: $chanE nsamp $nsamp maxamp $maxE\n" if ($debug > 0);
if ($pf{normal} == 0) { # max Z
	$ampV=$maxV;
	$ampN=$maxV;
	$ampE=$maxV;
} elsif ($pf{normal} == 1) { # Each to 1
	$ampV=$maxV;
	$ampN=$maxN;
	$ampE=$maxE;
} elsif ($pf{normal} == 3) { # max of 3 comps
	my @tmp=($maxV, $maxN, $maxE);
	my $amax=&find_amax(@tmp);
	$ampV=$amax;
	$ampN=$amax;
	$ampE=$amax;
} else { # do nothing
	$ampV=1;
	$ampN=1;
	$ampE=1;
}

open(OF,">$pf{wffile}") or die "$!: $pf{wffile}\n";

for (my $i=0; $i<$nsamp; $i++ ) {
	my $z=$ampV ? $dVert[$i]/$ampV: $dVert[$i];
	# Assume time vector same for all comps
	my $zt=$comps{$chanZ}{"tt"}[$i];
	my $n=$ampN ? $dNorth[$i]/$ampN: $dNorth[$i];
	my $nt=$comps{$chanN}{"tt"}[$i];
	my $e=$ampE ? $dEast[$i]/$ampE: $dEast[$i];
	my $et=$comps{$chanE}{"tt"}[$i];
	printf OF ("%0.4f %0.4f %0.4f %0.4f\n",$zt,$z,$n,$e) 
	;
}
close(OF);

#
# SUBROUTINES -------------------------------------------------

sub init()
{
	use Getopt::Std;
	my $opt_string = 'hvogrt:p:R:L:';
	getopts( "$opt_string", \%opt ) or usage();
	usage() if $opt{h};
}

sub usage()
{
print STDOUT << "EOF";
NAME
\t$progname - Extract 3-component data from a database and write to 
\ta file.

SYNOPSIS
\t$progname [-h|-v|-o|-g|-r] [-R angle] [-t time] -p parameter_file.pf

DESCRIPTION

\tExtract 3-comp data with time window around a pick 
\tin the arrival time. Need arrival assoc origin and wfdisc
\ttable. Initial join is origin->assoc->arrival. Returns t,Z,R,T?

\tChannels should be specified in ZNE order in the parameter file
\tunless -g option is used. Then the N and E components will be 
\tguessed according to the third letter in the channel code.
\tN,1,R,Y= North; E,2,T,X=East.  Z is only option for vertical.

\tOnly one of -o, -r, -R may be specified. -r and -R also take into
\taccount the hang value of the horizontal components. 

\tDefault filter is DEMEAN.

OPTIONS
\t-p parameter_file (REQUIRED; see example below)
\t-t time (OPTIONAL)
\t\tsimulated arrival time where pre/post time in pf relative to this
\t-o correct ZNE for hang/vang in sitechan table (OPTIONAL)
\t-g Guess which channels are Z,N,E based on chan code  (OPTIONAL)
\t-r Rotate into ZRT based on event SEAZ (OPTIONAL)
\t-R angle(deg) Rotate horizontals by angle (OPTIONAL)
\t-L baz,inc Rotate into LQT based on baz,inc
\t-h print this usage (OPTIONAL)
\t-v spew some debug messages (OPTIONAL)

EXAMPLE

\t$progname -v -o -p kkar.pf

BUGS
\tProbably creepy crawly ones.

EXAMPLE PARAMETER FILE
CUT HERE ---------------
dbname testdb   	# input database name
orid    1       	# orid number
staname KKAR    	# station name
phase   P       	# phase name (arrival/assoc/origin join)
chans   SHZ SHN SHE # chan to do
filter DEMEAN; BWZ 0.01 4 1.0 4 # Use trfilter
pretime  5      	# pre-singal time length (sec)
posttime 25     	# post-signal time length (sec)
pertaper 0.2  		#  (fractional) percent taper
normal   0          # normalize amps (0=max(Z),1=max(ZNE),2=no
wffile out.file	    # output waveform data
CUT HERE ---------------
EOF
die "\n";
}

sub read_parmfile {
	my $pfname=shift;
	my @keys=@_;
	my %pf;
	foreach my $k (@keys) {
		my $val=pfget($pfname,$k);
		if ( ! defined $val) {
			die "ERROR: need to set $k in $pfname\n";
		} else {
			$pf{$k}=$val;
		}
 	}
	return(%pf);
}

sub get_tr {
	my $t0=shift;
	my $t1=shift;
	my $sta=shift;
	my $chan=shift;
	my @db=@_;
	my $applycal=1; # applies calib if it's non-zero
	my @tr=trloadchan(@db,$t0,$t1,$sta,$chan,$applycal);
	my $ntr=dbquery(@tr,"dbRECORD_COUNT");
	if ($ntr == 0) {
		die "No records after trloadchan\n $sta $chan $t0 $t1\n";
	} 
	elsif  ($ntr > 1) {
		print STDOUT "1 records after trloadchan...trying trsplice\n";
		trsplice(@tr);
		$ntr=dbquery(@tr,"dbRECORD_COUNT");
	}
	die "ERROR: too many records($ntr) after trsplce\n" if ($ntr > 1);
	return(@tr);
}
sub get_arrt {
	my $orid=shift;
	my $sta=shift;
	my $phs=shift;
	my @db=@_;
	# Get start time for data subset
	my $sub="orid == $orid && sta=~/$sta/ && iphase=~/$phs/";
	print "SUB: $sub\n" if ($debug ==1);
	my @dbproc=dbprocess(@db, "dbopen arrival",
					"dbjoin assoc",
					"dbjoin origin",
					"dbsubset $sub");
	my $nrec=dbquery(@dbproc,"dbRECORD_COUNT");
	if ($nrec == 0) {
		die "No records after dbprocess\n" 
	} elsif ($nrec > 1) {
		print STDOUT "More then 1 record after dbprocess for start time--\n";
		print STDOUT "-- Not good, but proceeding anyway with 1st record.\n";
	}

	# Get arrival time for zero-offset and time sub window
	$dbproc[3]=0;
	my $arrt=dbgetv(@dbproc,qw(arrival.time));
	dbfree(@dbproc);
	return($arrt);
}
sub get_seaz {
	my $orid=shift;
	my $sta=shift;
	my $phs=shift;
	my @db=@_;
	# Get start time for data subset
	my $sub="orid == $orid && sta=~/$sta/ && iphase=~/$phs/";
	print "SUB: $sub\n" if ($debug ==1);
	my @dbproc=dbprocess(@db, "dbopen arrival",
					"dbjoin assoc",
					"dbjoin origin",
					"dbsubset $sub");
	my $nrec=dbquery(@dbproc,"dbRECORD_COUNT");
	if ($nrec == 0) {
		die "No records after dbprocess\n" 
	} elsif ($nrec > 1) {
		print STDOUT "More then 1 record after dbprocess for start time--\n";
		print STDOUT "-- Not good, but proceeding anyway with 1st record.\n";
	}

	# Get arrival time for zero-offset and time sub window
	$dbproc[3]=0;
	my $seaz=dbgetv(@dbproc,qw(assoc.esaz));
	dbfree(@dbproc);
	return($seaz);
}

sub normal1 {
	my @data=@_;
	my $n=scalar(@data);
	my $zmax=-10E10;
	for (my $i=0; $i<$n; $i++) {
		$zmax=$data[$i] if ($data[$i] > $zmax);
	}
	for (my $i=0; $i<$n; $i++) {
		$data[$i]=$data[$i]/$zmax;
	}
	return(@data);
}

sub correct_orientation {
	my (%comps) = %{$_[0]};
	my $mpi=3.1415926535;
	my $rpd=$mpi/180;
	my $ch=$comps{"chanz"};
    my $samprate=$comps{$ch}{"samprate"};
    my $nsamp=$comps{$ch}{"nsamp"};
    my @dataz=@{$comps{"dataz"}};
    my $hangz=$comps{"hangz"};
    my $vangz=$comps{"vangz"};
    my @datan=@{$comps{"datan"}};
    my $hangn=$comps{"hangn"};
    my $vangn=$comps{"vangn"};
    my @datae=@{$comps{"datae"}};
    my $hange=$comps{"hange"};
    my $vange=$comps{"vange"};
	my @dn_true = ();
	my @de_true = ();

	# From rot.c
	my ($cosz,$cosn,$cose,$sinn,$sine);
    if ($vangz < -360.) {
		$cosz = 1.;
	} else {
		$cosz = cos($vangz*$rpd);
	}
	if ($hangn < -360.) {
		$sinn = 0.;
		$cosn = 1.;
	} else {
		$sinn = sin($hangn*$rpd);
		$cosn = cos($hangn*$rpd);
	}
	if ($hange < -360.) {
		$sine = 0.;
		$cose = 1.;
	} else {
		$sine = sin($hange*$rpd);
		$cose = cos($hange*$rpd);
	}

	#print "hange $hange hangn $hangn \n";
	#print "sine $sine cose $cose sinn $sinn cosn $cosn \n";
	for ( my $i=0 ; $i<$nsamp ; $i++ ){
		$dataz[$i] *= $cosz;
		$dn_true[$i] = $cosn*$datan[$i] + $cose*$datae[$i];
		$de_true[$i] = $sinn*$datan[$i] + $sine*$datae[$i];
		#print "$dataz[$i] $dn_true[$i] $de_true[$i] \n";
	}

	#my $nz=scalar(@dataz);	
	#my $nn=scalar(@dn_true);	
	#my $ne=scalar(@de_true);	
	#print "ZZZ ----  $nz NNN $nn EEE $ne\n"; 
	# this syntax makes sure it's an output in the original hash
	@{$_[0]}{"dz_true"}=[ @dataz ];
	@{$_[0]}{"dn_true"}=[ @dn_true ];
	@{$_[0]}{"de_true"}=[ @de_true ];
	#$comps{"dz_true"}=[ @dataz ];
	#$comps{"dn_true"}=[ @dn_true ];
	#$comps{"de_true"}=[ @de_true ];

	printf "Done correcting for hang/vang...\n";
}

sub rot2d {
	use POSIX;
    # usages (@data)=rot2d(nsamp,angle,dNorth,dEast)
    # rotate 2 orthogonal components
    # Assumes data are North and East component data where
    # the vertical angle is +90 
    # and horizontal angle is +90 for East on 0 for North
    # returns the radial then tangential components
    my $nsamp=shift;
    my $ang=shift;
    $ang = $ang/57.296;
    my @data=@_;
    my (@N,@E);
    for (my $i=0; $i<$nsamp; $i++) {
        push(@N,$data[$i]);
        push(@E,$data[$i+$nsamp]);
    }
    my $rpd = 0.017453293;
    my $cosphi = cos($ang);
    my $sinphi = sin($ang);
    my $ehang=90;
    my $nhang=0;
    my $sine = sin($ehang*$rpd);
    my $cose = cos($ehang*$rpd);
    my $sinn = sin($nhang*$rpd);
    my $cosn = cos($nhang*$rpd);

    for (my $i=0 ; $i<$nsamp ; $i++ ) {
        my $dn_true = $cosn*$N[$i] + $cose*$E[$i];
        my $de_true = $sinn*$N[$i] + $sine*$E[$i];
        my $temp = $dn_true*$cosphi + $de_true*$sinphi;
        $E[$i] = $de_true*$cosphi - $dn_true*$sinphi;
        $N[$i] = $temp;          # Radial 
        $E[$i] = $E[$i]*(-1.0);   # Tangential
    }
    return(@N,@E);

}

sub costaper {
	use strict;
	use warnings;
	use POSIX;
	my $perc=pop;
	my @dd=@_;
	my $N=scalar(@dd);
	$perc=$perc/100.0;
	my $ntap=floor($perc*$N/2);
	my $pi=3.14159265358979;
	for (my $i=0; $i<$ntap; $i++) {
	    my $wt=0.5*(1-cos($i*$pi/$ntap));
		$dd[$i]*= $wt;
		$dd[$N-($i+1)]*= $wt;
	}
	return(@dd);
}

sub find_max {
    my $max=shift(@_);
    foreach (@_) {
        $max=$_ if ($max < $_);
    }
    return($max);
}

sub find_amax { # absolute max
    my $max=abs(shift(@_));
    foreach (@_) {
		my $a=abs($_);
        $max=$a if ($max < $a);
    }
    return($max);
}

sub rotate_ZNE_LQT {
	#@rotLQT=rotate_ZNE_LQT($nsamp,$baz,$inc,@Z,@N,@E);
    
    my $nsamp=shift;
    my $baz=shift;
    my $inc=shift;
    my @data=@_;
    my (@Z,@N,@E);
    my (@L,@Q,@T);
    for (my $i=0; $i<$nsamp; $i++) {
        push(@Z,$data[$i]);
        push(@N,$data[$i+$nsamp]);
        push(@E,$data[$i+2*$nsamp]);
    }
    my $rpd = 0.017453293;
    $baz = $baz*$rpd;
    $inc = $inc*$rpd;
    my $cosinc = cos($inc);
    my $sininc = sin($inc);
    my $cosbaz = cos($baz);
    my $sinbaz = sin($baz);

    for (my $i=0 ; $i<$nsamp ; $i++ ) {
        $L[$i] = $Z[$i] * $cosinc - $N[$i] * $sininc * $cosbaz - $E[$i] * $sininc * $sinbaz;
        $Q[$i] = $Z[$i] * $sininc + $N[$i] * $cosinc * $cosbaz + $E[$i] * $cosinc * $sinbaz;
        $T[$i] = $N[$i] * $sinbaz - $E[$i] * $cosbaz;
    }
    return(@L,@Q,@T);
}
