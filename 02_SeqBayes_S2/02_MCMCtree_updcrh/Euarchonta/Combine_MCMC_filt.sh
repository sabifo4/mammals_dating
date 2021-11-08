#!/bin/bash

curr_dir=$( pwd )


for dat in 7
do

	mkdir -p $dat/mcmc_files_filt 

	for i in 2 4 5 7 8 10 11 12 13 14 15 16 17 18 19 22 23 26 27 28 30 
	do 
		if [[ ! -f $dat/run$i/mcmctree_GBM/mcmc.txt ]]
		then 
			printf "Sorry, no samples for run"$i"/mcmc"$j"/mcmctree_GBM ...\n"
			printf "Data_"$dat"\trun"$i"\n" >> "Not_collected_samples.tsv"
		else 
			printf "Parsing dat "$dat" and run"$i"/mcmctree_GBM ... ... \n"
			end=$( wc -l $dat/run$i/mcmctree_GBM/mcmc.txt | sed 's/ ..*//' )
			if [[ $i -eq 1 ]]
			then
				begin=1
			else 
				begin=2 
			fi
			sed -n ''${begin}','${end}'p' $dat/run$i/mcmctree_GBM/mcmc.txt >> $dat/mcmc_files_filt/mcmc_tracer.txt
		fi
	done 

done 


# NOTE:
# After this script, I copied dummy_aln, ctl file, and tree file 
# to generate the FigTree file using the option -1 !