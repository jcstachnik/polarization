eval 'exec perl -S $0 ${1+"$@"}'
if 0;

use lib "$ENV{ANTELOPE}/data/perl";
use lib "$ENV{HOME}/lib/perl";
use Datascope;
use strict 'vars';
use warnings;
use vars qw/%opt/;
require "wfextract_module.pl";
#require "find_max.pl";
#require "my_math.pl";
use POSIX;

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
my %comps=();
my %info_hash;

%info_hash = (
'opt' => \%opt,
'pf' => \%pf, 
'debug' => $debug,
'correct_orient' => $correct_orient,
'do_arrt' => $do_arrt,
'opttime' => $opttime,
'guesschan' => $guesschan,
'rotate' => $rotate,
'rotate_angle' => $rotate_angle,
'baz' => $baz,
'inc' => $inc
			);

if ($debug>2) {
print " Premodule -------------------- \n";
while ( my ($key, $value) = each(%info_hash) ) {
	print "$key => $value \n";
}
}


&wfextract_module(\%info_hash);


if ($debug>2) {
print " Postmodule ------------------- \n";
while ( my ($key, $value) = each(%info_hash) ) {
	print "$key => $value \n";
}
print "Postmodule -------------------- \n";
}
# At this point, one could change the rotate parameters and 
# outfile name, call the module again, and go along your merry way.
# comps should be filled out with the data arrays as well.
#
# This is an example of how to get the data arrays out of 
# the wfextract_module
my $new_comp=$info_hash{comps}; # reference to hash comps
my $chanZ=$new_comp->{chanz};
my @T=@{$new_comp->{$chanZ}->{"tt"}};
my @Z=@{$new_comp->{$chanZ}->{"data"}};

if ($debug>2) {
while ( my ($key, $value) = each(%$new_comp) ) {
	print "new_comp  $key => $value \n";
}
}
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
filter DEMEAN; BWZ 1.0 4 0 0 
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

