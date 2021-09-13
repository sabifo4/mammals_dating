#!/bin/bash

# =============================================================================== #
# This script goes through all the `partitions12_*aln` files saved in             #
# the `filtered_genes_step2` directory during the previous R filtering step.      #
# Then, it carries out the following tasks:                                       #
#                                                                                 #
#  * Skips first line (phylip header)                                             #
#  * Prints the alignment to `combined_all_mammals_for_perl_partX.aln`,           #
#    where X goes from 1 to 4 as there are 4 partitions.                          #
#  * Then adds an extra blank line which will be used by the perl script          #
#    that is subsequently run to indicate the end of one gene alignment.          #
#                                                                                 #
#                                                                                 #
# You are expected to run this script from main dir,                              #
# i.e., `00_Gene_filtering` as:                                                   #
#                                                                                 #
# ./scripts/04.1_Concatenate_genes_separated_for_partition_subsamples.sh <nlim> \ #
#  <namedir> <numparts>                                                           #
# =============================================================================== #
# Any questions/doubts/bugs, please send a message to:                            #
# Sandra Alvarez-Carretero, <s.alvarez-carretero@ucl.ac.uk>                       #
# =============================================================================== #


# 1. Set counters and vars
i=0
jinc=$1 # E.g., 3859
j=$1    # E.g., 3859
count=0
count_per_parts=0
namedir=$2
numparts=$3 

# 2. Move to main dir, where this script saved in `00_Gene_filtering/scripts`
#    is run, i.e., `00_Gene_filtering`
curr_dir=$( pwd )
cd $curr_dir

# 3. Create `00_alignments_subsamples` if not there 

if [ ! -d 00_alignments_subsamples ] 
then 
mkdir 00_alignments_subsamples 
fi 

if [ ! -d 00_alignments_subsamples/$namedir/parts_$numparts/ ]
then 
mkdir -p 00_alignments_subsamples/$namedir/parts_$numparts/
fi 


# 4. Combine all the different fasta alignment in a unique
#    fasta file for each partition
for partition in `seq 1 $numparts`
do 
	# Create folder for each partition 
	mkdir -p 00_alignments_subsamples/$namedir/parts_$numparts/part$partition
	
	# If it is not the last partition
	if [ $partition != $numparts ]
	then
																		
		while [ $i -lt $j ]
		do
			i=$(( i + 1 ))
			count=$(( count + 1 ))
			count_per_parts=$(( count_per_parts + 1 ))
			name_gene=$( echo `ls $namedir/${i}/*aln` | sed 's/..*\/partitions12\_//' | sed 's/\.aln//'  )
			cp $namedir/${i}/*aln 00_alignments_subsamples/$namedir/parts_$numparts/part$partition/gene$count"_"$name_gene.fasta
			sed -i -e '1,1d' 00_alignments_subsamples/$namedir/parts_$numparts/part$partition/gene$count"_"$name_gene.fasta          # Remove phylip header
			sed -i -e 's/^/\>/' 00_alignments_subsamples/$namedir/parts_$numparts/part$partition/gene$count"_"$name_gene.fasta       # Add ">" at the beginning of each sequence 
			sed -i -e 's/\s\{6\}/\n/' 00_alignments_subsamples/$namedir/parts_$numparts/part$partition/gene$count"_"$name_gene.fasta # Replace spaces between taxa names and seqs with new line
			echo Gene $count visited
			echo Gene $count visited >> 00_alignments_subsamples/$namedir/parts_$numparts/log_$namedir.txt
		done 
		
		# Update counters for genes
		printf "Total of genes for partition "$partition": "$count_per_parts"\n"
		printf "Total of genes for partition "$partition": "$count_per_parts"\n" >> 00_alignments_subsamples/$namedir/parts_$numparts/log_$namedir.txt
		i=$( echo $count )
		j=$(( i + jinc ))
		count_per_parts=$( echo 0 )
	
	# Otherwise, just reach until the end of genes available 
	# in the directory
	else
		tot=$( ls $namedir | wc -l )
		k=$(( tot - 2 )) # Now I need to subtract 2 because I have 2 txt files!
		# Set i as the last gene !
		i=$( echo $count )
		while [ $i -lt $k ]
		do
			i=$(( i + 1 ))
			count=$(( count + 1 ))
			name_gene=$( echo `ls $namedir/${i}/*aln` | sed 's/..*\/partitions12\_//' | sed 's/\.aln//'  )
			cp $namedir/${i}/*aln 00_alignments_subsamples/$namedir/parts_$numparts/part$partition/gene$count"_"$name_gene.fasta
			# Convert to fasta-phylip-partitions needed input format (fasta format)
			sed -i -e '1,1d' 00_alignments_subsamples/$namedir/parts_$numparts/part$partition/gene$count"_"$name_gene.fasta          # Remove phylip header
			sed -i -e 's/^/\>/' 00_alignments_subsamples/$namedir/parts_$numparts/part$partition/gene$count"_"$name_gene.fasta       # Add ">" at the beginning of each sequence 
			sed -i -e 's/\s\{6\}/\n/' 00_alignments_subsamples/$namedir/parts_$numparts/part$partition/gene$count"_"$name_gene.fasta # Replace spaces between taxa names and seqs with new line
			echo Gene $count visited
			echo Gene $count visited >> 00_alignments_subsamples/$namedir/parts_$numparts/log_$namedir.txt
			count_per_parts=$(( count_per_parts + 1 ))
		done 
		printf "Total of genes for partition "$partition": "$count_per_parts"\n"
		printf "Total of genes for partition "$partition": "$count_per_parts"\n" >> 00_alignments_subsamples/$namedir/parts_$numparts/log_$namedir.txt
	
	fi
	
	printf "Partition "${partition}" finishes !\n"

done 

# Now generate additional directory for concatenated 
if [ ! -d 00_alignments_subsamples/$namedir/conc ]
then 
mkdir -p 00_alignments_subsamples/$namedir/conc
for partition in `seq 1 $numparts`
do 
	cp -R 00_alignments_subsamples/$namedir/parts_$numparts/part$partition/*fasta 00_alignments_subsamples/$namedir/conc
done 
fi 

