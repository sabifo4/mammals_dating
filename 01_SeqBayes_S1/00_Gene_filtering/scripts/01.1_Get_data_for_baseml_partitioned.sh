#!/bin/bash

# ========================================================================== #
# This script goes through every directory inside `00_Gene_filtering/baseml` #                             #
# and partitions the corresponding gene alignment according to gene codon    #
# partition scheme.                                                          #
#                                                                            #
# You are expected to run this script from main dir,                         #
# i.e., `00_Gene_filtering` as:                                              #
#                                                                            #
# ./scripts/01.1_Get_data_for_baseml_partitioned.sh                          #
#                                                                            #
# This scripts outputs two alignment files,`partitions12.3_*aln` and         #
# `partitions12_*aln`, in the corresponding directory for each gene.         #
#                                                                            #
# ========================================================================== #
# Any questions/doubts/bugs, please send a message to:                       #
# Sandra Alvarez-Carretero, <s.alvarez-carretero@ucl.ac.uk>                  #
# ========================================================================== #

# 1. Get path to main dir,`00_Gene_filtering/scripts`,
#    and then get the path to `src`. Move to baseml
curr_dir=$( pwd )
cd $curr_dir
cd ../../src
src=$( pwd )
cd $curr_dir/baseml

# 2. Run `src/partition_alignment` for each of the genes 
#    in the `baseml` directory
for i in *
do
cd $i
printf "Parsing "$i" ...\n\n"
# The usage of this perl script is the following:
#   <path_to_script>/partition_alignments.pl <alignment_file> <lines_to_skip> <separator>
$src/partition_alignments.pl ENS*[0-9].aln 2 "\s{6}"
cd ..
done 

## NOTE: This script might take a while to run if it is not run as a job
##       array in a cluster, so you might want to adapt this code accordingly.
##
## Example of bash script to be submited to the Apocrita HPC (QMUL) below.
##
## >> SCRIPT GOES FROM LINES 50-83. REMOVE "#> " AT THE BEGINNING OF EACH LINE, COPY
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
#> cd $curr_dir
#> cd ../../src
#> src=$( pwd )
#> cd $curr_dir/baseml
#> 
#> # 2. Run `src/partition_alignment` for each of the genes 
#> #    in the `baseml` directory. We can access each file 
#> #    by assigning each directory name to the $SGE_TASK_ID
#> i=$( echo $SGE_TASK_ID )
#> cd $i
#> printf "Parsing "$i" ...\n\n"
#> # The usage of this perl script is the following:
#> #   <path_to_script>/partition_alignments.pl <alignment_file> <lines_to_skip> <separator>
#> $src/partition_alignments.pl ENS*[0-9].aln 2 "\s{6}"
#> 