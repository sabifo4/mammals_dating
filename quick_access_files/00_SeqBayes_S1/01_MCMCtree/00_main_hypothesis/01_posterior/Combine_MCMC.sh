#!/bin/bash

curr_dir=$( pwd )

mkdir mcmc_files
for i in 1 3 4 
do 
	printf "Parsing run"$i"/mcmctree_GBM ... ... \n"
	end=$( wc -l run$i/mcmctree_GBM/mcmc.txt | sed 's/ ..*//' )
	if [[ $i -eq 1 ]]
	then
		begin=1
	else 
		begin=2 
	fi
	sed -n ''${begin}','${end}'p' run$i/mcmctree_GBM/mcmc.txt >> mcmc_files/mcmc_tracer.txt
	sed -n '1,'${end}'p' run$i/mcmctree_GBM/mcmc.txt >>  run$i/mcmctree_GBM/mcmc_clean.txt
done 
 
 


# NOTE:
# After this script, I copied dummy_aln, ctl file, and tree file 
# to generate the FigTree file using the option -1 !