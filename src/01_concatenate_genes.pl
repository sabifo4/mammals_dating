#!/usr/bin/perl

use strict;
use warnings;
 
## Script that concatenates the genes of different species 
## in a unique line
##
## Usage:
## <path_to_script>/concatenate_genes.pl <alignment_file>
## Contact information: <sandra.ac93@gmail.com>

#IMPORTANT! CREATE ".TXT" FILES IN UNIX, NOT WINDOWS!!!!! IT THEN ADDS THE "\r" CARRIAGE AND IT MAKES THE SCRIPT CRASH...
		
# Open the input file(s)
open (ALIGNMENTS, "<$ARGV[0]") or die "Cannot open $ARGV[0] file: $!";

# Get the name of the string from the input file
my $string = $ARGV[0];
$string =~ s/\_.*//;
chomp($string);

# Open the output file to save the sequences divided according to the partitions 
#open(OUT1, ">concatenated_tab_$string.aln") or die "Cannot create the output file: $!";
open(OUT2, ">$string.aln") or die "Cannot create the output file: $!";
open(OUT3, ">$string.fasta") or die "Cannot create the output file: $!";
open(OUT4, ">$string.log.txt") or die "Cannot create the output file: $!";

# Get lines of the file
my @alignments = <ALIGNMENTS>;

# Create variables 
my $species = "";
my $sequence = "";
my $count = 0;
my %species_hash = ();
my %species_hash2 = ();

foreach my $line (@alignments){

	chomp($line);
	# Get name of species and sequence 
	# The format of a line in this file is "<SPECIES>\t<SEQUENCE>"
	my @split_line = split(/\t/, $line);
	$species = $split_line[0];
	# Remove part of string that contains the specific chromosome
	# it was gotten from 
	$species =~ s/\-chr..*//g;
	chomp($species);
	$sequence = $split_line[1];
	chomp($sequence);
	
	# If species name not in hash, add it and save 
	# first gene sequence 
	if ( ! exists $species_hash{$species} ){
		$species_hash{$species} = $sequence."\n";
		$species_hash2{$species} = $sequence;
		$count += 1;
		#print $species_hash{$species} = $sequence."\n";
	}
	else{
		# Append the sequence to the corresponding 
		# species
		$species_hash{$species} .= $sequence."\n";
		$species_hash2{$species} .= $sequence;
	}
	
}												   

print OUT4 "Number of species: ".$count."\n";

# Get header for alignment in phylip format with concatenated genes 
# and print out
my $length_aln = 0;
my $firstkey = ( sort( keys %species_hash2) )[0];
$length_aln = length( $species_hash2{$firstkey} );

print OUT2 $count." ".$length_aln."\n";
print OUT2 " \n";

# Print out file(s) with sequences 
foreach $species (sort keys %species_hash){
	#print $species."\n";
	#print OUT1 "$species\t$species_hash{$species}\n";
	print OUT2 "$species      $species_hash2{$species}\n";
	print OUT3 ">$species\n$species_hash{$species}";
	print OUT4 "Length for species ".$species." is: ".length( $species_hash2{$species} )."\n";
	
} 

# Close file(s)
#close(OUT1);
close(OUT2);
close(OUT3);
close(OUT4);
close(ALIGNMENTS);
