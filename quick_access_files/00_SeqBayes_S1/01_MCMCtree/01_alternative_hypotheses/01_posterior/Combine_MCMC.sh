#!/bin/bash

curr_dir=$( pwd )

for dat in `seq 1 10`
do

	mkdir -p mcmc_files 
	cd mcmc$dat/mcmctree_GBM/
	for tt in *
	do 
		if [[ ! -d $curr_dir/mcmc_files/$tt ]] 
		then 
			mkdir -p $curr_dir/mcmc_files/$tt
		fi
		printf "Parsing mcmc"$dat"/mcmctree_GBM/"$tt"/mcmc.txt ... ... \n"
		end=$( wc -l $tt/mcmc.txt | sed 's/ ..*//' )
		if [[ $dat -eq 1 ]]
			then
				begin=1
			else 
				begin=2 
			fi
		sed -n ''${begin}','${end}'p' $tt/mcmc.txt >> $curr_dir/mcmc_files/$tt/mcmc_tracer.txt
		# sed -n '1,'${end}'p' $dat/run$i/mcmc$j/mcmctree_GBM/mcmc.txt >>  $dat/run$i/mcmc$j/mcmctree_GBM/mcmc_clean.txt
	done
	
	cd $curr_dir

done 


# NOTE:
# After this script, I copied dummy_aln, ctl file, and tree file 
# to generate the FigTree file using the option -1 !