#!/usr/bin/perl

use strict;
use warnings;
 
## Script that reads a file in FASTA format and extracts 12 and 3 CPs.
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
$outname =~ s/\_one\_line\.fasta//;
$outname =~ s/..*\///;
print "Parsing alignment ".$outname." ... ... \n";
my $outname2 = $outname."_out12CP.fasta";
my $outname3 = $outname."_out3CP.fasta";
## Open the output file to save the alignment file 
open(OUT, ">$outname2") or die "Cannot create the output file: $!";
open(OUT2, ">$outname3") or die "Cannot create the output file: $!";

## Get lines of input file
my @aln1 = <INFILE1>;

## Create variables
my $species = "";
my $seq_cp12 = "";
my $seq_cp3 = "";
my $count = 0;
my $count_sp = 0;
my @bases = "";

my $count12cp = 0;
my $length12cp = 0;
my $length3cp = 0;
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
			$species =~ s/\n//;
			$species =~ s/\r//;
			chomp( $species );
			# In the script inside `test3 I did not have to 
			# add a new line. NOW I have to...
			# Lines 63 and 64 are then different from the other
			# script version
			print OUT $species."\n";
			print OUT2 $species."\n";
		}
		$count_sp += 1;
	}
	else{
		# Get sequences and print OUT in a line
		@bases = split( //, $line );
		foreach my $base (@bases){
			chomp( $base );
			#print $base;
			$count12cp += 1;
			if ( $count12cp < 3 ){
				$seq_cp12 .= $base;
			}else{
				# Restart counter 
				$count12cp = 0;
				$seq_cp3 .= $base;
			}
		}
		$seq_cp12 =~ s/\n//;
		$seq_cp12 =~ s/\r//;
		$seq_cp3 =~ s/\n//;
		$seq_cp3 =~ s/\r//;
		$length12cp = length( $seq_cp12 );
		$length3cp = length( $seq_cp3 );
		$species =~ s/>//;
		print OUT $seq_cp12."\n";
		print OUT2 $seq_cp3."\n";
		print "\nSpecies parsed: ".$species."\n"."Length of sequence with 12CP: ".$length12cp."\n";
		print "Length of sequence with 3CP: ".$length3cp."\n\n";
	}
	
	

}

print "This is the 12cp partition: ".$seq_cp12."\n";
print "This is the 3cp partition: ".$seq_cp3."\n";
print "A total of ".$count_sp." species have been parsed!\n";
# Add new line at the end of the sequence		

## Close files
close(OUT);
close(INFILE1);
