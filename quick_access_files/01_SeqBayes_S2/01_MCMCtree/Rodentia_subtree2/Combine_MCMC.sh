#!/bin/bash

curr_dir=$( pwd )


for dat in 7
do

	mkdir -p $dat/mcmc_files 

	for i in `seq 1 32` 
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
			sed -n ''${begin}','${end}'p' $dat/run$i/mcmctree_GBM/mcmc.txt >> $dat/mcmc_files/mcmc_tracer.txt
			sed -n '1,'${end}'p' $dat/run$i/mcmctree_GBM/mcmc.txt >>  $dat/run$i/mcmctree_GBM/mcmc_clean.txt 
		fi
	done 

done 


# NOTE:
# After this script, I copied dummy_aln, ctl file, and tree file 
# to generate the FigTree file using the option -1 !