#!/bin/bash

# =========================================================================== #
# This script goes through every gene alignment inside the `baseml` directory #
# and extracts the sequences for mouse and human. Then, it saves them in      #
# a new file which will be later used by an R script to calculate the         #
# sequence distance                                                           #
#                                                                             #
# You are expected to run this script from main dir,                          #
# i.e., `00_Gene_filtering` as:                                               #
#                                                                             #
# ./scripts/01.2_Get_mouse_human_aln.sh                                       #
#                                                                             #
# =========================================================================== #
# Any questions/doubts/bugs, please send a message to:                        #
# Sandra Alvarez-Carretero, <s.alvarez-carretero@ucl.ac.uk>                   #
# =========================================================================== #

# Go through each alignment, extract mouse and human sequences 
# and create file to be read in R  
for i in baseml/*
do

# 1. Get gene name 
gene_name=$( ls $i/*tree | sed 's/..*\///' | sed 's/\..*//' )
echo Parsing $gene_name
echo Parsing $gene_name >> out_logs/log_03_get_mouse_human.txt

#> # 2. Create the mouse_human file to be later read in R 
grep 'mus_musculus' $i/partitions12_*aln > $i/$gene_name"_mouse_human.aln"
grep 'homo_sapiens' $i/partitions12_*.aln >> $i/$gene_name"_mouse_human.aln"
sed -i 's/      /\t/g' $i/$gene_name"_mouse_human.aln"

done 

## NOTE: This script might take a while if it is not run as a job array in a cluster.
##       You might want to copy lines XXX and paste them in a bash script to be 
##       submitted in a cluster as a job array. 
##
## E.g. of bash script to be submited to the Apocrita HPC (QMUL)
##
## >> SCRIPTS GOES FROM LINES 54-89. REMOVE "#> " AT THE BEGINNING OF EACH LINE, COPY
## >> THE CODE SNIPPET, PASTE IT IN A NEW FILE, SAVE IT AS A BASH SCRIPT, AND 
## >> READY TO BE SUBMITTED TO THE APOCRITA HPC.

#> #!/bin/bash
#> #$ -cwd                     # Run the code from the current directory
#> #$ -V                       # Export environment to job
#> #$ -j y                     # Merge the standard output and standard error
#> #$ -l h_rt=240:00:00        # Limit each task to 10 days
#> #$ -l h_vmem=1G             # Request 1GB RAM
#> #$ -t 1-15569
#> 
#> # NOTE: Remember to save this script in the folder `scripts` as explained 
#> #       in the tutorial as well as the scripts in the `src` directory.
#> #       This way it follows the file architecture described in the tutorial 
#> #       and the subsequent code will work without any issues. You might want 
#> #       to transfer the whole GitHub repository to the cluster to ensure 
#> #       this happens.
#> 
#> # 1. Get path to main dir,`00_Gene_filtering/scripts`,
#> #    and then get the path to `src`. Move to baseml
#> #    Note that $curr_dir would be `00_Gene_filtering/scripts`.
#> curr_dir=$( pwd )
#> cd $curr_dir/baseml
#>
#> # 2. Go through each alignment, extract mouse and human sequences 
#> #    and create file to be read in R. We access each gene by using the 
#> #    $SGE_TASK_ID 
#> i=$( echo $SGE_TASK_ID )
#> 
#> # 2.1. Get gene name 
#> gene_name=$( ls $i/*tree | sed 's/..*\///' | sed 's/\..*//' )
#> echo Parsing $gene_name >> out_logs/tmp3_$SGE_TASK_ID".txt"
#> 
#> # 2.2. Create the mouse_human file to be later read in R 
#> grep 'mus_musculus' $i/partitions12_*aln > $i/$gene_name"_mouse_human.aln"
#> grep 'homo_sapiens' $i/partitions12_*.aln >> $i/$gene_name"_mouse_human.aln"
#> sed -i 's/      /\t/g' $i/$gene_name"_mouse_human.aln"
#> 
#> # 3. Check if all genes have been visited and, if so, 
#> #    remove all the tmp files after concatenating them into a single file
#> num_genes=(`ls out_logs/tmp3*`)
#> visited_genes=$( echo ${#num_genes[@]} )
#> if [ $visited_genes -eq 15569 ]
#> then 
#>
#> for f in out_logs/tmp3*
#> do 
#> cat $f >> out_logs/log_03_get_mouse_human.txt 
#> printf "\n" >> out_logs/log_03_get_mouse_human.txt 
#> rm $f 
#> done 
#>
#> fi
