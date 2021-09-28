#!/usr/bin/perl

use strict;
use warnings;
 
## Script that reads the concatenated sequences generated with 
## UNIX paste command in dir 1 and then generates a FASTA file 
## with the species names
## Usage:
## <path_to_script>/Get_partitions.pl <path_to_conc_file> <path_to_sp>
## Contact information: <sandra.ac93@gmail.com>

#IMPORTANT! CREATE ".TXT" FILES IN UNIX, NOT WINDOWS!!!!! IT THEN ADDS THE "\r" CARRIAGE AND IT MAKES THE SCRIPT CRASH...
		
## Open the input files
open (INFILE1, "<$ARGV[0]") or die "Cannot open $ARGV[0] file: $!";
open (INFILE2, "<$ARGV[1]") or die "Cannot open $ARGV[1] file: $!";
## Get second argument to set output file name
my $outname = $ARGV[0];
chomp($outname);
$outname =~ s/\.txt//;
$outname =~ s/..*\///;
print "Parsing alignment ".$outname." ... ... \n";
my $outname2 = $outname.".fasta";
my $outname3 = $outname."_outnogaps.fasta";
## Open the output file to save the alignment file 
open(OUT, ">$outname2") or die "Cannot create the output file: $!";
open(OUT2, ">$outname3") or die "Cannot create the output file: $!";

## Get lines of input file
my @seq = <INFILE1>;
my @sp = <INFILE2>;

## Create variables 
my $tmpseq = "";
my $nogapsseq = "";
my $count = -1;
my $count_sp = 0;
my $count_sp_nogaps = 0;
## Loop over all the lines in input file 
foreach my $line (@sp){

	chomp($line);
	$count += 1;
	$count_sp += 1;
	# Access array and extract the specific sequence
	$tmpseq = $seq[$count];
	chomp( $tmpseq );
	$nogapsseq = $tmpseq;
	$nogapsseq =~ s/\-//g;
	print OUT ">".$line."\n";
	print OUT $tmpseq."\n";
	#print length( $nogapsseq )."\n";
	if( length( $nogapsseq ) == 0 ){
		print "I should not be here, all taxa should have at least one gene!\n";
		# We will use OUT4 later to append gaps once the 
		# alignment has been performed
		$count_sp_nogaps += 1;
	}else{ 
		print OUT2 ">".$line."\n";
		print OUT2 $nogapsseq."\n";
	}
	

}
					
print "Total no. species parsed: ".$count_sp."\n";
print "Total no. species parsed without gaps: ".$count_sp_nogaps."\n";

## Close files
close(OUT);
close(OUT2);
close(INFILE1);
close(INFILE2);
