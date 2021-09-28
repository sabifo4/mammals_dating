#!/usr/bin/perl

use strict;
use warnings;
 
## Script that reads a FAST alignment and generates an output file 
## called <name_inp>_seqs.txt. This contains only the sequences 
## without the header (i.e., no ">species" tag before the sequence).
## This file is used to generate later the concatenated file.
##
## Usage:
## <path_to_script>/Get_partitions.pl <path_to_aln1> 
## Contact information: <sandra.ac93@gmail.com>

#IMPORTANT! CREATE ".TXT" FILES IN UNIX, NOT WINDOWS!!!!! IT THEN ADDS THE "\r" CARRIAGE AND IT MAKES THE SCRIPT CRASH...
		
## Open the input files
open (INFILE1, "<$ARGV[0]") or die "Cannot open $ARGV[0] file: $!";
## Get second argument to set output file name
my $outname = $ARGV[0];
chomp($outname);
$outname =~ s/\.fasta//;
$outname =~ s/..*\///;
print "Parsing alignment ".$outname." ... ... \n";
my $outname2 = $outname."_seqs.txt";
## Open the output file to save the alignment file 
open(OUT, ">$outname2") or die "Cannot create the output file: $!";

## Get lines of input file
my @aln1 = <INFILE1>;

## Create variables 
my $count = 0;

## Loop over all the lines in input file 
foreach my $line (@aln1){

	chomp($line);
	
	# Ommit header
	if( $line =~ /^>/ ){
		$count += 1;
		print OUT "\n";
	}
	else{
		# Get sequence 
		print OUT $line."\n";
	}
}
					
print "Total no. species parsed: ".$count."\n";

## Close files
close(OUT);
close(INFILE1);