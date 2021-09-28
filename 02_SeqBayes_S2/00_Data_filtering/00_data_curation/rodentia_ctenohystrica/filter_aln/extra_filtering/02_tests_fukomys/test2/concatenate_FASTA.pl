#!/usr/bin/perl

use strict;
use warnings;
 
## Script that reads a file with several sequences in FASTA format for the same 
## taxa and then concatenates them in a unique line.
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
my $outname2 = $outname."_out.fasta";
## Open the output file to save the alignment file 
open(OUT, ">$outname2") or die "Cannot create the output file: $!";

## Get lines of input file
my @aln1 = <INFILE1>;

## Create variables
my $species = "";
my $count = 0;
my $count_sp = 0;

## Loop over all the lines in input file 
foreach my $line (@aln1){

	chomp($line);
	$count += 1;
	
	# Print header 
	if( $line =~ /^>/ ){
		# Write header if line == 1 
		if ( $count == 1 ){
			$species = $line;
			$species =~ s/,.*//;
			$species =~ s/..*[0-9] //;
			$species =~ s/ /\_/g;
			print OUT ">".$species."\n";
		}
		$count_sp += 1;
	}
	else{
		# Get sequences and print OUT in a line
		print OUT $line;
	}

}

print "A total of ".$count_sp." fragments for this species have been parsed!\n";
# Add new line at the end of the sequence		
print OUT "\n\n";

## Close files
close(OUT);
close(INFILE1);
