#!/bin/bash

# ===================================================================== #
# This script goes through each directory in the                        #
# `filtered_genes_step2_all72sp_all72sp` directory and finds out which  #
# of these genes have also been included in the data set that is to be  #
# analysed in the second step of this sequential Bayesian dating        #
# analysis.                                                             #
#                                                                       #
# You are expected to run this script from main dir,                    #
# i.e., `00_Gene_filtering` as:                                         #
#                                                                       #
# ./scripts/02.2_Remove_nucmit_genes_2ndstep_common72sp.sh              #
#                                                                       #
# ===================================================================== #
# Any questions/doubts/bugs, please send a message to:                  #
# Sandra Alvarez-Carretero, <s.alvarez-carretero@ucl.ac.uk>             #
# ===================================================================== #

# Start counter
count_out=0

# Find main dir from which this script is 
# run, i.e., `00_Gene_filtering`
# and move to `filtered_genes_step2_all72sp`
curr_dir=$( pwd )
cd $curr_dir/filtered_genes_step2_all72sp

# Check if gene included in second data set to analyse 
# is present and, if so, mark to be deleted
genes_check=$( cat $curr_dir/scripts/genes_to_check.txt )
i=$( echo $i | sed 's/\r//' | sed 's/\n//' )
for j in */*tree 
do
name=$( echo $j | sed 's/..*\///' | sed 's/\.tree//' )
dir=$( echo $j | sed 's/\/..*//' )
for i in $genes_check
do
if [[ $i =~ $name ]]
then 
printf "Gene "$name" found in dir "$dir"! This will be deleted in step 3\n" >> $curr_dir/out_logs/log_13_removed_genes_common72sp.txt
count_out=$(( count_out + 1 ))
printf $count_out": gene in dir "$dir", "$name", will be deleted!\n"
elif [[ ! $i =~ $name ]]
then
cp -R $curr_dir/filtered_genes_step2_all72sp/$dir $curr_dir/filtered_genes_step3_all72sp
fi
done
done

printf "Filtering step finished!\n\n"
printf "Removed genes: "$count_out"\n\n"