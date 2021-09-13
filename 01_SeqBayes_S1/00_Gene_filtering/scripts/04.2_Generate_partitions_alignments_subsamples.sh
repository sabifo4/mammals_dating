#!/bin/bash

# ============================================================================= #
# This script partitions the alignments using an in-house pipeline.             #
#                                                                               #
# You are expected to run this script from main dir,                            #
# i.e., `00_Gene_filtering` as:                                                 #
#                                                                               #
# ./scripts/04.2_Generate_partitions_alignments_subsamples.sh <namedir>         #                                                                 
# ============================================================================= #
# Any questions/doubts/bugs, please send a message to:                          #
# Sandra Alvarez-Carretero, <s.alvarez-carretero@ucl.ac.uk>                     #
# ============================================================================= #


# 1. Set counters and vars
curr_dir=$( pwd )
namedir=$1
cd 00_alignments_subsamples/$namedir 
home_dir=$( pwd )

# 2. Move to 01_SeqBayes_S1/00_Gene_filtering/00_alignments_subsamples
#    and run the code

## 2.1. Concatenated alignment 
printf "We are working with the concatenated alignment ... ... \n"
cd $home_dir/conc
cp $curr_dir/00_alignments/part1/species_names.txt .
Run_tasks.sh . mammals_$namedir"_conc" partN

## 2.2. Two partitions 
printf "We are working with the 2 partitions alignment ... ... \n"
cd $home_dir/parts_2

for i in `seq 1 2` 
do 
cd part$i 
cp $curr_dir/00_alignments/part1/species_names.txt .
Run_tasks.sh . mammals_$namedir"_2parts" partN
cd ..
done

cd $home_dir/parts_2
for i in `seq 1 2`
do
cd part$i/phylip_format/02_concatenated_alignments 
cat *aln >> $home_dir/parts_2/$namedir"_2parts.aln"
cd $home_dir/parts_2
done 
 
## 2.3. Four partitions 
printf "We are working with the 4 partitions alignment ... ... \n"
cd $home_dir/parts_4

for i in `seq 1 4` 
do 
cd part$i 
cp $curr_dir/00_alignments/part1/species_names.txt .
Run_tasks.sh . mammals_$namedir"_4parts" partN
cd ..
done

cd $home_dir/parts_4
for i in `seq 1 4`
do
cd part$i/phylip_format/02_concatenated_alignments 
cat *aln >> $home_dir/parts_4/$namedir"_4parts.aln"
cd $home_dir/parts_4
done 

printf "FINISHED!\n\n"
