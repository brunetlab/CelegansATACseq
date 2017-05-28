#! /usr/bin/perl

use warnings;
use strict;
use IO::CaptureOutput qw(capture_exec);

# a script to downsample to desired size
unless (scalar @ARGV == 3) {
	die "\nNot enough command line arguments.\n".
	"Usage : random_split_fastq.pl <original bed file> <output bed 1> <output bed 2>\n";
}

# create an array that contains the list of files to be treated
my $original = shift @ARGV ;
die "cannot open $original.\n" unless(open(ORIG, "<$original"));

my $out1 = shift @ARGV;
my $out2 = shift @ARGV;

die "cannot open $out1.\n" unless(open(OUT1, ">$out1"));
die "cannot open $out2.\n" unless(open(OUT2, ">$out2"));

srand(time|$$);

while(my $line = <ORIG>) {

	my $random = rand(10);
	
	if ($random < 5) {
		print OUT1 $line;		
	} else {
		print OUT2 $line;
	}
}


close ORIG;
close OUT1;
close OUT2;

exit;
