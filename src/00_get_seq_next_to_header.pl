#!/usr/bin/perl

use strict;
use warnings;

# Open the input file(s)
open (ALN, "<$ARGV[0]") or die "Cannot open $ARGV[0] file: $!";
# Get the name of the string from the input file
my $string = $ARGV[0];
$string =~ s/\..*/\_tab.aln/;
chomp($string);
# Open the output file to save the sequences aligned to their tags
open(OUT, ">$string") or die "Cannot create the output file: $!";

# Save input file(s) in array 
my @aln_inp = <ALN>;

# Loop over each line in file to 
# align each seq to each header
foreach my $line (@aln_inp){
	if ($line =~ /^>/){
		$line =~ s/\n/\t/;
		$line =~ s/>//;
		print OUT $line;
	}
	else{
		print OUT $line;
	}
}

# Close file(s) 
close(ALN);
close(OUT);