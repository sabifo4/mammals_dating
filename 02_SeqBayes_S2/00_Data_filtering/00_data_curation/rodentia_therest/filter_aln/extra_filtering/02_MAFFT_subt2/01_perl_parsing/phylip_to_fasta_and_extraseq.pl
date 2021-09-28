#!/usr/bin/perl

use strict;
use warnings;
 
## Script that reads alignment 1 (main), alignment 2 (secondary), and a file with 
## the list of taxa to be extracte from alignmen 2. Then, it attaches 
## taxon/taxa from aln2 to aln1 and outputs a FASTA file keeping the gaps,
## another without gaps. Additional output files are:
##  - sequences.txt: file with only the sequences with gaps, no species 
##  - species.txt: files with the taxa names in the alignment 
##  - sum_len_aln.txt: file with a summary of the length of the alignment 
##  - taxa_without_sequences.txt: file with a list of taxa which sequence was only 
##                                gaps, i.e., no data
##
## Usage:
## <path_to_script>/Get_partitions.pl <path_to_aln1> <path_to_aln2> <list_taxa_to_add>
## Contact information: <sandra.ac93@gmail.com>

#IMPORTANT! CREATE ".TXT" FILES IN UNIX, NOT WINDOWS!!!!! IT THEN ADDS THE "\r" CARRIAGE AND IT MAKES THE SCRIPT CRASH...
		
## Open the input files
open (INFILE1, "<$ARGV[0]") or die "Cannot open $ARGV[0] file: $!";
open (INFILE2, "<$ARGV[1]") or die "Cannot open $ARGV[1] file: $!";
open (INFILE3, "<$ARGV[2]") or die "Cannot open $ARGV[2] file: $!";
## Get second argument to set output file name
my $outname = $ARGV[0];
chomp($outname);
$outname =~ s/\.aln//;
$outname =~ s/..*\///;
print "Parsing alignment ".$outname." ... ... \n";
my $outname2 = $outname."_out.fasta";
my $outname3 = $outname."_outnogaps.fasta";
my $outname4 = $outname."_outforMAFFT.fasta";
if (! -e 'forMAFFT') {
 mkdir('forMAFFT') or die "Can't create forMAFFT:$!\n";
}
## Open the output file to save the alignment file 
open(OUT, ">$outname2") or die "Cannot create the output file: $!";
open(OUT2, ">$outname3") or die "Cannot create the output file: $!";
open(OUT3, ">>sum_len_aln.txt") or die "Cannot create the output file: $!";
open(OUT4, ">>taxa_without_sequences.txt") or die "Cannot create the output file: $!";
open(OUT5, ">>species.txt") or die "Cannot create the output file: $!";
open(OUT6, ">>sequences.txt") or die "Cannot create the output file: $!";
open(OUT7, ">>forMAFFT/$outname4") or die "Cannot create the output file: $!";
open(OUT8, ">>forMAFFT/taxatoadd_outforMAFFT.fasta") or die "Cannot create the output file: $!";

## Get lines of input file
my @aln1 = <INFILE1>;
my @aln2 = <INFILE2>;
my @add_fasta = <INFILE3>;

## Create variables 
my $species = "";
my $count = 0;
my $count_sp = 0;
my $length_seq = 0;
my $fasta_seq = "";
my $fasta_seq_nogaps = "";

my $lenseq1 = 0;
my $lenseq2 = 0;
## Loop over all the lines in input file 
foreach my $line (@aln1){

	chomp($line);
	$count += 1;
	$line =~ s/\r//;
	
	# Print header 
	if( $line =~ /^[0-9]/ ){
		my @split_seq = split(/  /, $line);
		$length_seq = $split_seq[1];
		## Print header 
		print "Length file 1 = ".$length_seq."\n";
	}
	else{
		# Get name species in alignment 
		my @split_line = split(/   /, $line);
		$species = $split_line[0];
		$fasta_seq = $split_line[1];
		$fasta_seq_nogaps = $split_line[1];
		# Remove gaps as we are going to realign everything again
		$fasta_seq_nogaps =~ s/\-//g;
		chomp( $species );
		chomp( $fasta_seq );
		chomp( $fasta_seq_nogaps );
		#print length( $fasta_seq_nogaps )."\n";
		
		$lenseq1 = length( $fasta_seq );
		$lenseq2 = length( $fasta_seq_nogaps );
		print OUT3 "Species :".$species." Length: ".$lenseq1."\n";
		print OUT3 "Species :".$species." Length: ".$lenseq2."\n\n";
		
		# Print OUT in FASTA format. Add only taxa in "no gaps"
		# output file if there is at least one gene available
		print OUT ">".$species."\n";
		print OUT $fasta_seq."\n";
		print OUT7 ">".$species."\n"; # for MAFFT 
		print OUT7 $fasta_seq."\n";   # for MAFFT
		
		print OUT5 $species."\n";
		print OUT6 $fasta_seq."\n";
		# If length == 0 or == 1, that is <2
		if( length( $fasta_seq_nogaps ) < 2 ){
			print OUT4 $species."\n";
			# We will use OUT4 later to append gaps once the 
			# alignment has been performed
		}else{ 
			print OUT2 ">".$species."\n";
			print OUT2 $fasta_seq_nogaps."\n";
		}
		
		$count_sp += 1;
	}

}
					
print "Total no. species parsed: ".$count_sp."\n";
					
# Parse second alignment 
my $name_taxa = "";
my $match_taxa = "";
my $counter_tree2 = -scalar(@add_fasta); # Do not count headers for each taxon added!
foreach my $taxa (@add_fasta){
	chomp( $taxa );
	$name_taxa = $taxa;
	foreach my $seq (@aln2){
		chomp( $seq );
		$seq =~ s/\r//;
		$counter_tree2 += 1;
		my @split_aln2 = split( /   /, $seq );
		$match_taxa = $split_aln2[0];
		chomp( $match_taxa );
		if( $name_taxa =~ /\Q$match_taxa\E/ ){
			
			$species = $split_aln2[0];
			$fasta_seq = $split_aln2[1];
			$fasta_seq_nogaps = $split_aln2[1];
			# Remove gaps as we are going to realign everything again
			$fasta_seq_nogaps =~ s/\-//g;
			chomp( $species );
			chomp( $fasta_seq );
			chomp( $fasta_seq_nogaps );
			
			$lenseq1 = length( $fasta_seq );
			$lenseq2 = length( $fasta_seq_nogaps );
			print OUT3 "Species :".$species." Length: ".$lenseq1."\n";
			print OUT3 "Species :".$species." Length: ".$lenseq1."\n\n";
			print OUT5 $species."\n";
			print OUT6 $fasta_seq."\n";
			
			# Print OUT in FASTA format
			print OUT ">".$species."\n";
			print OUT $fasta_seq."\n";
			print OUT8 ">".$species."\n"; # for MAFFT
			print OUT8 $fasta_seq_nogaps."\n";   # for MAFFT
			
			print OUT2 ">".$species."\n";
			print OUT2 $fasta_seq_nogaps."\n";
			$count_sp += 1;
			
		}
	}
}

# Adjust counter so only the num of 
# true iterations per taxa to add remains
if ( scalar( @add_fasta ) > 1 ){
	$counter_tree2 = $counter_tree2/scalar( @add_fasta );
}
print "Total no. species in tree 2: ".$counter_tree2."\n";
print "Total no. species parsed after adding new taxa: ".$count_sp."\n";


## Close files
close(OUT);
close(OUT2);
close(OUT3);
close(OUT4);
close(OUT5);
close(OUT6);
close(INFILE1);
close(INFILE2);
close(INFILE3);
