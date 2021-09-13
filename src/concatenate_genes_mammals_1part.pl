#!/usr/bin/perl

use strict;
use warnings;
 
## Script that concatenates the genes of different species 
## in a unique line
##
## Usage:
## <path_to_script>/concatenate_genes.pl <alignment_file> <species_list> <separator>
## Contact information: <sandra.ac93@gmail.com>

#IMPORTANT! CREATE ".TXT" FILES IN UNIX, NOT WINDOWS!!!!! IT THEN ADDS THE "\r" CARRIAGE AND IT MAKES THE SCRIPT CRASH...
		
# 1. Open the input file(s)
open ( ALIGNMENTS, "<$ARGV[0]" ) or die "Cannot open $ARGV[0] file: $!";
open ( SPECIES, "<$ARGV[1]" ) or die "Cannot open $ARGV[1] file: $!";

# 2. Get the name of the string from the input file
my $string = $ARGV[0];
$string =~ s/\_.*//;
chomp( $string );

my $separator = $ARGV[2];
$separator =~ s/\"//;
$separator =~ s/\'//;
chomp( $separator );

# 3. Open the output file to save the sequences divided according to the partitions 
#open(OUT1, ">concatenated_tab_$string.aln") or die "Cannot create the output file: $!";
open(OUT2, ">$string.aln") or die "Cannot create the output file: $!";
open(OUT3, ">$string.fasta") or die "Cannot create the output file: $!";
open(OUT4, ">$string.log.txt") or die "Cannot create the output file: $!";
open(OUT5, ">$string.log2.txt") or die "Cannot create the output file: $!";

# Get output for log2 
print OUT5 "Gene"."\t"."Number of species"."\t"."Number of missing species"."\t"."Total genes = 72"."\t"."Name of missing sp"."\n";

# 4. WORKING WITH FILE WITH SPECIES LIST

# 4.1. Get lines of file with list of species
#      and num of species  
my @list_species = <SPECIES>;
my $num_sp = scalar( @list_species );
#print "My number of species is: ".$num_sp."\n";

# 4.2. Create variables 
my %species_hash = ();
my %species_hash2 = ();
my %visited_species = ();

# 4.3. Declare keys in hash with undefined values
#      and get amount of species
for my $sp ( @list_species ){
	# Add the f* chomp that caused me so 
	# many problems...
	chomp( $sp );
	$species_hash{ $sp } = undef;
	$species_hash2{ $sp } = undef;	
	$visited_species{ $sp } = 0; 
}
#print "I have ".scalar(@list_species)." species in the list\n";

# 5. PROCESS ALIGNMENT 

# 5.1. Get lines of the alignment
my @alignments = <ALIGNMENTS>;

# 5.2. Create variables 
my $species = "";
my $sequence = "";
my $gaps = "";

my $count = 0;
my $count_sp = 0;
my $count_genes = 0;

my $length_aln = 0;
my $length_gaps = 0;
my @missing_sp = ();

# 5.3. Get concatenated data 
foreach my $line (@alignments){

	# Start counting lines
	$count += 1;
	
	# Check that it is not a blank line.
	# If it is a blank line, then this means we have finished to 
	# read one gene, thus increase by one gene counter
	# and add gaps to the species in which gene is missing
	if( $line =~ /^\s*$/ ){
		# We have one more gene!
		print "I am here at line ".$count."\n";
		$count_genes += 1;
		
		# foreach my $element ( keys %visited_species ) {
			# print $visited_species{$element}
		# }
		# Check which species did not have that gene in the alignment 
		# if there are not 72 sp and add gaps with same length
		if( $count_sp != $num_sp ){
			
			## NOTE: The next commented 8 lines werejust written when
			## I forgot a chomp above. Just keeping the code for fun!
			# Find which species have not been visited (they have value = 0 )
			#foreach my $val ( sort keys %visited_species ){
			#	print $val."  ".$visited_species{$val}."\n\n";
				# if( $visited_species{$val} == 0 ){
					# push @missing_sp, $val;
					# print $val."  ".$visited_species{$val}."\n";
				# }
			#}
			# Use the grep that was not working when I forgot the chomp...
			@missing_sp = grep { $visited_species{$_} == 0 } keys %visited_species;
			# foreach my $val ( @missing_sp ){
				# print $val."";
			# }
			print "For gene ".$count_genes." we have ".$count_sp." species\n";
			print "There are ".scalar(@missing_sp)." genes missing!\n\n";
			print OUT5 $count_genes."\t".$count_sp."\t".scalar(@missing_sp)."\t".($count_sp + scalar(@missing_sp))."\t".join(", ", @missing_sp)."\n";
			# Add gaps to this species entry
			for my $spnogene ( @missing_sp ){
				$species_hash{ $spnogene } .= $gaps."\n";
				$species_hash2{ $spnogene } .= $gaps;				
			}			
		}
		else{
			print OUT5 $count_genes."\t".$count_sp."\t".scalar(@missing_sp)."\n";
		}
		
		# Reset value of species & @missing_sp
		$count_sp = 0 ;
		@missing_sp = ();
		# Set values of keys in %visited_species to 0 (restart for new gene)
		for my $sp ( @list_species ){ 
			$visited_species{ $sp } = 0; 
		}
	}
	else{
	
		# Remove any weird character
		chomp( $line );
		
		# Count me this line as a species sequence!
		$count_sp += 1;
		
		# Get name of species and sequence 
		# The format of a line in this file is "<SPECIES>$separator<SEQUENCE>"
		my @split_line = split( /$separator/, $line );
		$species = $split_line[0];
		chomp( $species );
		$sequence = $split_line[1];
		chomp( $sequence );
		
		#print "I am fucking here!\n";
		# Get length of sequence only once and generate an 
		# object with "n" gaps 
		if ( $count_sp == 1 ){
			# Get length of sequence 
			$length_gaps = length( $sequence );
			#print $length_gaps."\n\n";
			$gaps = "-" x $length_gaps;
			#print $gaps."\n\n\n";
			#$length_aln = length( $species_hash2{ $species } );
		}
		# Append the sequence to the corresponding species
		# and count one visit for that species!
		$species_hash{ $species } .= $sequence."\n";
		$species_hash2{ $species } .= $sequence;
		## NOTE: I needed to do this because of the 
		## chomp I forgot above... Keeping it as it might 
		## be useful for future scripts
		#delete $visited_species{$species};
		$visited_species{ $species } = 1;
		
	}
}												   

# Get header for alignment in phylip format with concatenated genes 
# and print out
my $firstkey = ( sort( keys %species_hash2) )[0];
$length_aln = length( $species_hash2{$firstkey} );

print OUT2 $num_sp." ".$length_aln."\n";
print OUT2 " \n";

# Paste sequences next to species for this partition
foreach my $sp_in_hash (sort keys %species_hash){
	#print $sp_in_hash."\n";
	#print OUT1 "$sp_in_hash\t$species_hash{$sp_in_hash}\n";
	print OUT2 "$sp_in_hash      $species_hash2{$sp_in_hash}\n";
	print OUT3 ">$sp_in_hash\n$species_hash{$sp_in_hash}";
	print OUT4 "Length for species ".$sp_in_hash." is: ".length( $species_hash2{$sp_in_hash} )."\n";
	print "Length for species ".$sp_in_hash." is: ".length( $species_hash2{$sp_in_hash} )."\n";
			
}

# Close file(s)
#close(OUT1);
close(OUT2);
close(OUT3);
close(OUT4);
close(OUT5);
close(ALIGNMENTS);
