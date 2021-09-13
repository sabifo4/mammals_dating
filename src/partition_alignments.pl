#!/usr/bin/perl

use strict;
use warnings;
 
## Script that outputs the alignment of X species partitioned 
## by codon positions. The first partition will contain 1st and 
## 2nd CPs while the second partition will have the 3rd CPs.
##
## Usage:
## <path_to_script>/partition_alignments.pl <alignment_file> <lines_to_skip> <separator>
## Contact information: <sandra.ac93@gmail.com>

#IMPORTANT! CREATE ".TXT" FILES IN UNIX, NOT WINDOWS!!!!! IT THEN ADDS THE "\r" CARRIAGE AND IT MAKES THE SCRIPT CRASH...
		
# Open the input files
open (ALIGNMENTS, "<$ARGV[0]") or die "Cannot open $ARGV[0] file: $!";
# Get the name of the string from the input file
my $string = $ARGV[0];
$string =~ s/\..*//;
chomp($string);
# Get the number of lines to skip that do not contain sequences
my $skip_lines = $ARGV[1];
chomp( $skip_lines );
# Get the kind of separator 
my $separator = $ARGV[2];
$separator =~ s/\"//;
$separator =~ s/\'//;
chomp( $separator );
# Print message with options to keep track
print "Your alignment file is ".$ARGV[0]."\n";
print "You want to skip ".$skip_lines." lines from your alignment file.\n";
print "You separate species from sequences using the separator ".$separator."\n";
# Open the output file to save the sequences divided according to the partitions 
open(OUT, ">partitions12.3_$string.aln") or die "Cannot create the output file: $!";
open(OUT2, ">partitions12_$string.aln") or die "Cannot create the output file: $!";

# Get lines of the file
my @alignments = <ALIGNMENTS>;

# Create variables 
my $species = "";
my $sequence = "";
my @nucleotides = ();
my $part1 = "";
my $part2 = "";
my $count = 0 ;
my $count_sp = 0;
my $length_seq = 0;

foreach my $line (@alignments){
	
	chomp($line);
	# Initialize counter
	$count += 1;
	# Skip $skip_lines:
	if ( $count > $skip_lines ){
	
		# Skip blank lines
		if( $line =~ /^\s*$/ ){
			next;
		}
		
		# Count species 
		$count_sp += 1;
		
		# Get name of species and sequence 
		# The format of a line in this file is "<SPECIES>$separator<SEQUENCE>"
		#my @split_line = split(/\t/, $line);
		my @split_line = split(/$separator/, $line);
		$species = $split_line[0];
		chomp($species);
		$sequence = $split_line[1];
		chomp($sequence);
		@nucleotides = split("", $sequence);
		#print $nucleotides[0];
		#print length($sequence);
		
		# Append the species name at the begining of the sequence
		# If it is the first sequence...
		if ( $count == 1 ){
			$part1 = $species."\t";
			$part2 = $species."\t";
			$length_seq = length($sequence);
		}
		# If not, just concatenate to the previous string
		else{
			$part1 .= $species."\t";
			$part2 .= $species."\t";
		}
		# Go through each nucleotide and find the 1st, 2nd, and 3rd CPs 
		# The 1st and 2nd CPs will be separated into one partition ($part1)
		# The 3rd CP will be separated in a second partition ($part2)
		my $count_nucs = 0;
		for ( my $i = 0; $i < (length($sequence)); $i += 1 ){
			# 1. Start counting nucleotides
			$count_nucs += 1;

			# 2.1. If $count_name == 1 | 2 ( 1st & 2nd CPs )
			if ( $count_nucs == 1 || $count_nucs == 2 ){
				$part1 .= $nucleotides[$i];
			} 
			# 2.2. If $count_name == 3 (3rd CP), then restart $count_nucs
			#      and put nucleotide in second partition
			elsif ( $count_nucs == 3 ){
				$part2 .= $nucleotides[$i];	
				$count_nucs = 0;
			}
			
			# 3. If it is the last position and this is a 1st or 2nd CP...
			#    NOTE: The nucleotide has already been added in previous conditions,
			#    therefore only a space is added as the sequence is all read
			if ( $i == (length($sequence)-1) ){
				$part1 .= "\n";
				$part2 .= "\n";
			}
		
		}
	
	}
	
}												   

# Get length of the sequences in each partition
my @part1_div = split('\t', $part1);
my @part1_div2 = split('\n', $part1_div[1]);
my $seq_length_part1 = length($part1_div2[0]);

my @part2_div = split('\t', $part2);
my @part2_div2 = split('\n', $part2_div[1]);
my $seq_length_part2 = length($part2_div2[0]);

# Replace tab that separates species name from 
# sequence with 6 spaces to meet PAML format 
$part1 =~ s/\t/      /g;
$part2 =~ s/\t/      /g;

# Print the content of the partitions in the output file
# in PAML format
print OUT $count_sp." ".$seq_length_part1."\n".$part1."\n\n".$count_sp." ".$seq_length_part2."\n".$part2."\n";
print OUT2 $count_sp." ".$seq_length_part1."\n".$part1;

# Close files
close(OUT);
close(OUT2);
close(ALIGNMENTS);
