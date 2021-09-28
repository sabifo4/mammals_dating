#!/usr/bin/perl

use strict;
use warnings;
 
## Script that reads one alignment and finds which sequencs
## are to be kept if compared to a list of given sequences
##
## Usage:
## <path_to_script>/Get_partitions.pl <text_file> <alignment_file> <num_taxa>
## Contact information: <sandra.ac93@gmail.com>

#IMPORTANT! CREATE ".TXT" FILES IN UNIX, NOT WINDOWS!!!!! IT THEN ADDS THE "\r" CARRIAGE AND IT MAKES THE SCRIPT CRASH...
		
## Open the input files
open (INFILE1, "<$ARGV[0]") or die "Cannot open $ARGV[0] file: $!";
open (INFILE2, "<$ARGV[1]") or die "Cannot open $ARGV[1] file: $!";
## Get second argument to set output file name
my $outname = $ARGV[1];
chomp($outname);
$outname =~ s/\.aln//;
$outname =~ s/..*\///;
print "Parsing alignment ".$outname." ... ... \n";
my $outname2 = $outname."_out.fasta";
$outname = $outname."_out.aln";
## Open the output file to save the alignment file 
open(OUT, ">$outname") or die "Cannot create the output file: $!";
open(OUT2, ">>sum_len_aln.txt") or die "Cannot create the output file: $!";
open(OUT3, ">$outname2") or die "Cannot create the output file: $!";
open(OUT4, ">>total_sp_per_aln.csv") or die "Cannot create the output file: $!";

## Get lines of input file
my @names = <INFILE1>;
my @alignment = <INFILE2>;

## Create variables 
my %species_hash = ();
my $count = -1;
my $species = "";
my $count_sp = 0;
my $count_sp2 = 0;
my $count_sp3 = 0;
my $length_seq = 0;
my $fasta_seq = "";

## Loop over text file and save name of species in the hash 
print "Parsing text file with names...\n";
foreach my $line (@names){
	
	chomp( $line );
	$line =~ s/\r//;
	$count_sp += 1;
	# Create entry in hash with $line as key and "1" as entry
	$species_hash{$line} = 1;
	#print "Species ".$line." parsed!\n";
	
}
print "A total of ".$count_sp." have been parsed in the text file!\n";



## Loop over all the lines in input file 
foreach my $line (@alignment){

	chomp($line);
	$count += 1;
	
	# Get info of length and save it 
	if( $count == 0 ){
		my @split_seq = split(/  /, $line);
		$length_seq = $split_seq[1];
		## Print header 
		print OUT $ARGV[2]."  ".$length_seq."\n";
	}
	
	# Skip 0 as this is the header
	else{
		# Get name species in alignment 
		my @split_line = split(/   /, $line);
		$species = $split_line[0];
		$fasta_seq = $split_line[1];
		# Remove gaps as we are going to realign everything again
		$fasta_seq =~ s/\-//g;
		chomp( $species );
		chomp( $fasta_seq );
		
		# If species name in hash, add it and write 
		# the line in output file 
		## NOTE> "exists" works with arrays and hashes only
		##       For strings use "defined"
		if ( exists $species_hash{$species} ){
			print OUT $line."\n";
			if( length( $fasta_seq ) > 1 ){
				print OUT3 ">".$species."\n";
				print OUT3 $fasta_seq."\n";
				$count_sp3 += 1;
			}
			#print "Species ".$species." added!\n";
			$count_sp2 += 1;
		}
	}
	
}
								   
## Print info about sequences parsed
print "\nTotal of sequences parsed: ".$count_sp2."\n";
print OUT2 "Output file:".$ARGV[1]."\tNum.taxa: ".$count_sp2."\n";

# Add the 4 new taxa added in the previous loop
$count_sp3 += 4;
print OUT4 $ARGV[1].",".$count_sp3."\n";

## Close files
close(OUT);
close(OUT2);
close(OUT3);
close(OUT4);
close(INFILE1);
close(INFILE2);
