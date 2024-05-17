#!/bin/bash

curr_dir=$( pwd )
mkdir -p mcmc_files

for i in `seq 1 5`
do

	cd mcmc$i
	printf "Parsing mcmc"$i"/mcmc.txt ... ... \n"
	end=$( wc -l mcmc.txt | sed 's/ ..*//' )
	if [[ $i -eq 1 ]]
		then
			begin=1
		else 
			begin=2 
		fi
	sed -n ''${begin}','${end}'p' mcmc.txt >> $curr_dir/mcmc_files/mcmc_tracer.txt
	sed -n '1,'${end}'p' mcmc.txt >>  mcmc_clean.txt
	cd $curr_dir
done


# NOTE:
# After this script, I copied dummy_aln, ctl file, and tree file 
# to generate the FigTree file using the option -1 !