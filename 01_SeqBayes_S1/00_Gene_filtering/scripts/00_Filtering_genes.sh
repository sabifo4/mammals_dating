#!/bin/bash

# ===================================================================== #
# This script goes through each fasta file within the `genes` directory #
# to generate a summary file with the following information:            #
#                                                                       #
#   * Gene name                                                         #
#   * Number of sequences included in the alignment                     #
#   * Presence of _Mus musculus_ in the alignment                       # 
#   * Presence of _Homo sapiens_ in the alignment                       #
#   * Number of codons the sequences in the alignment have              #
#   * Number of nucleotides the sequences in the alignment have         #
#                                                                       #
# You are expected to run this script from main dir,                    #
# i.e., `00_Gene_filtering` as:                                         #
#                                                                       #
# ./scripts/00_Gene_filtering.sh                                        #
#                                                                       #
# ===================================================================== #
# Any questions/doubts/bugs, please send a message to:                  #
# Sandra Alvarez-Carretero, <s.alvarez-carretero@ucl.ac.uk>             #
# ===================================================================== #


# Start counter
count=0

# Move to main dir, where this script saved in `00_Gene_filtering/scripts`
# is run, i.e., `00_Gene_filtering`
curr_dir=$( pwd )
cd $curr_dir
cd ../../src
src=$( pwd )
cd $curr_dir

# Create log file
printf "Gene,Num_seqs,Mouse_presence,Human_presence,Num_codons_hum,Length_seq\n" > out_data/genes.csv

# Loop over genes (15,904 genes) and filter
for i in genes/*
do

# Start general counter
count=$(( count + 1 ))

# If the cds.aln.fasta exists...
if [ -f $i/*cds.aln.fasta ]
then

# Get gene name
gene=$( echo $i | sed 's/..*\///')
# Check num sequences
num_seq=$( grep '>' $i/*cds.aln.fasta | wc -l )
# Check mus musculus is there 
mouse=$( grep '>mus_musculus' $i/*cds.aln.fasta | wc -l )
# Check homo sapiens is there 
homo=$( grep '>homo_sapiens' $i/*cds.aln.fasta | wc -l )
# Check length of codons (1 codon = 3 bps)
# The perl scripts gets the fasta sequences in one line 
# so it is easier to count the amount of nucleotides,
# then this file is removed.
$src/one_line_fasta.pl $i/*cds.aln.fasta
num_nucs=$( wc -L $i/*one_line.fa | awk '{print $1}' )
num_codons=$( grep -A1 'homo_sapiens' $i/*one_line.fa | sed 's/\-//g' | sed -n '2,2p' | wc -L | awk '{print $1/3}' )
rm $i/*_one_line.fa

else 
gene=$( echo NULL )
num_seq=$( echo NULL )
mouse=$( echo NULL )
homo=$( echo NULL )
num_nucs=$( echo NULL )

fi 

echo There are $count genes visited to generate summary statistic
printf $gene","$num_seq","$mouse","$homo","$num_codons","$num_nucs"\n" >> out_data/genes.csv
printf $gene","$num_seq","$mouse","$homo","$num_codons","$num_nucs"\n"

done











