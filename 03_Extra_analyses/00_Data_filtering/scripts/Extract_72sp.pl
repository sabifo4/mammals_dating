#!/usr/bin/perl

use strict;
use warnings;
use Cwd qw(getcwd);
use List::MoreUtils qw(uniq);

## Script that reads a text file with the names of the taxa 
## that need to be extracted from an alignment and generates 
## an output alignment only with them. It will generate a tmp 
## file with those taxa that are left to find in other 
## directories. See the for loop that you will need to use to run 
## this script properly below.
##
## Usage:
## <path_to_script>/Extract_72sp.pl <path_to_aln1> <path_to_text> <out_dir>
##
## NOTE: If the user wants to run this in a loop for several
##       directories, an example is given below: 
##
##       ```  
##       # Run from <subtree> directory
##       cd <main_dir>
##       home_dir=$( pwd )
##       cp ../scripts/taxa_72sp.txt ../scripts/tmp_taxa_72sp.txt
##       for i in */<NUM>
##       do 
##       cd $i
##       printf "<< DIR "$i" >>\n"
##       perl ../../../scripts/Extract_72sp.pl *aln ../../../scripts/tmp_taxa_72sp.txt
##       mv tmp_taxa_72sp.txt ../../../scripts
##       cd $home_dir 
##       done
##       cd $home_dir
##       # Remove tmp file 
##       rm ../scripts/tmp_taxa_72sp.txt
##       ```
##
##       This code snippet assumes that you have a PHYLIP file as an input 
##       file.
##
## Contact information: <sandra.ac93@gmail.com>

#IMPORTANT! CREATE ".TXT" FILES IN UNIX, NOT WINDOWS!!!!! IT THEN ADDS THE "\r" CARRIAGE AND IT MAKES THE SCRIPT CRASH...
		
## Open the input files
open (ALN, "<$ARGV[0]") or die "Cannot open $ARGV[0] file: $!";
open (TXT, "<$ARGV[1]") or die "Cannot open $ARGV[1] file: $!";
## Get second argument to set output file name
my $outname = $ARGV[0];
chomp($outname);
$outname =~ s/..*\///;
$outname =~ s/\.aln/\_out_72sp\.aln/;
my $outname2 = $outname;
$outname2 =~ s/\.aln/\_nogaps\.fasta/;
print "Parsing alignment ".$outname." ... ... \n\n";
my $outdir = "out_aln";
## Open the output file to save the alignment file 
unless( mkdir $outdir ){
	die "Cannot create the output file: $outdir\n";
}
open(OUT, ">$outdir/$outname") or die "Cannot create the output file: $!";
open(OUT2, ">$outdir/log_taxa.txt") or die "Cannot create the output file: $!";
open(OUT3, ">tmp_taxa_72sp.txt") or die "Cannot create the output file: $!";
open(OUT4, ">$outdir/$outname2") or die "Cannot create the output file: $!";
print OUT2 "Species\tSeq_length\n";

## Get lines of alignment and txt file
my @aln = <ALN>;
my @taxa = <TXT>;

## Create global variables
my $tmp_taxon = "";
my $tmp_seq = "";
my $tmp_seq2 = "";
my $tmp_seq_fasta = "";
my $tmp_nucs = "";
my $tmp_seqlength = 0;
my $count_sp = 0;
my $count_it = 0;
my $to_check = 0;
my @taxa2 = @taxa;
## Loop over all the 72sp names and try to find 
## them in the taxa available for every subtree.
## This order is preferred as some trees are > 72 
## so, if not found, search will finish earlier
foreach my $line (@taxa){
	
	chomp($line);
	$line =~ s/\r//;
	$line =~ s/\n//;
	#print $line."\n";
	foreach my $aln_line (@aln){
		
		chomp( $aln_line );
		# Skip blank lines!
		if( $aln_line =~ /^$/ ){
			next;
		}
		# Now just keep the first bit 
		# with taxon name
		$tmp_taxon = $aln_line; 
		$tmp_taxon =~ s/ ..*//;
		#print $tmp_taxon."\n";
		$tmp_taxon =~ s/\r//;
		$tmp_taxon =~ s/\n//;
		#if ( $line =~ /\Q$tmp_taxon\E/ ){ #Issues with exact match, hence not used
		if ( $line eq $tmp_taxon ){
			print "Species found: ".$tmp_taxon."\n";
			$tmp_seq .= $aln_line."\n";
			$tmp_seq2 = $aln_line;
			$tmp_seq2 =~ s/-//g;
			$tmp_seq2 =~ s/..* //;
			$tmp_seq_fasta .= ">".$tmp_taxon."\n".$tmp_seq2."\n";
			$count_sp += 1;
			# Get the seq length only with first taxon
			if ( $count_sp == 1 ){
				$tmp_nucs = $aln_line;
				$tmp_nucs =~ s/..* //;
				$tmp_nucs =~ s/\n//;
				$tmp_nucs =~ s/\r//;
				$tmp_seqlength = length( $tmp_nucs );
				print OUT2 $tmp_taxon."\t".$tmp_seqlength."\n";
			}
			#print "Species before grep: ".scalar(@taxa2)."\n";
			@taxa2 = grep {$_ !~ $line} @taxa2;
			#print "Species after grep: ".scalar(@taxa2)."\n";
		}
		
	}
	
}
			
## Print alignments
print OUT $count_sp."  ".$tmp_seqlength."\n";
print OUT $tmp_seq;
print OUT4 $tmp_seq_fasta;

print "\nThere are ".$count_sp." taxa with length ".$tmp_seqlength."\n";
print "Output file generated and saved as ".$outdir."/".$outname."\n\n";

## Output the tmp list with the remaining ones
#my @taxa2 = uniq @taxa2;
#print "There are: ".scalar(@uniquealn2)." remaining taxa\n\n";
foreach my $line (@taxa2){
	chomp( $line );
	print OUT3 $line."\n";
}

## Close files
close( OUT );
close( OUT2 );
close( OUT3 );
close( OUT4 );
close( TXT );
close( ALN );