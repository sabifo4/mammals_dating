#!/bin/bash

# ======================================================================== #
# This script goes through every directory for all the genes filtered      #
# in the previous step and does the following:                             #
#                                                                          #
#   1. Within a folder called `baseml`, it creates directories labelled    #
#      as `1`, `2`, `n`; where `n` equals to the number of directories     #
#      with filtered genes. There, the corresponding gene and tree files   #
#      will be saved.                                                      #
#   2. It gets the gene name and number of sequences and saves them in a   #
#      variable.                                                           #
#   3. It gets the fasta alignments in one line format using               #
#      `../../src/one_line_fasta.pl` script.                               #
#   4. It uses the script `../../src/00_get_seq_next_to_header.pl` to      #
#      reformat the output from step 3 so the sequences are next to the    #
#      species name tab separated.                                         #
#   5. It uses the script `../../src/01_concatenate_genes.pl` to now       #
#      generate an alignment in PHYLIP format.                             #
#   6. It removes unnecessary files generated during steps 3-5.            #
#   7. It gets the RAxML best-scoring ML tree and saves it in the          #
#      corresponding `baseml/$num` directory in PAML format.               #
#                                                                          #
#In the end, each `baseml/$SGE_TASK_ID` directory has 4 files:             #
#   1. *aln: The alignment in PHYLIP format.                               #
#   2. *fasta: The alignment in FASTA format.                              # 
#   3. *log.txt: A log file with the length of each sequence in the        #
#      alignment.                                                          #
#   4. *tree: The tree file in PHYLIPormat.                                #
#                                                                          #
# You are expected to run this script from main dir,                       #
# i.e., `00_Gene_filtering` as:                                            #
#                                                                          #
# ./scripts/01.1_Get_data_for_baseml.sh                                    #
#                                                                          #
# ======================================================================== #
# Any questions/doubts/bugs, please send a message to:                     #
# Sandra Alvarez-Carretero, <s.alvarez-carretero@ucl.ac.uk>                #
# ======================================================================== #

# 1. Start counter
count=0

# 2. Move to main dir, where this script saved in `00_Gene_filtering/scripts`
#    is run, i.e., `00_Gene_filtering`
curr_dir=$( pwd )
cd $curr_dir
cd ../../src
src=$( pwd )
cd $curr_dir

# 3. Loop over filtered genes and prepare them for baseml 
#    analyses
for i in filtered_genes/*
do

	# 3.0. Start general counter & create dir
	count=$(( count + 1 ))
	mkdir -p baseml/$count

	# 3.1 Get gene name  
	gene=$( echo $i | sed 's/..*\///')
	# 3.2 Check num sequences
	num_seq=$( grep '>' $i/*cds.aln.fasta | wc -l )

	echo Parsing gene $count":" $gene 
	echo Parsing gene $count":" $gene >> out_logs/log_02_baseml_getdataformat.txt

	# 3.3 Get one line sequences
	$src/one_line_fasta.pl $i/*cds.aln.fasta
	mv $i/*one_line.fa baseml/$count
	# 3.4. Get sequence next to header 
	$src/00_get_seq_next_to_header.pl baseml/$count/*one_line*
	# 3.5. Get alignment in PHYLIP format 
	$src/01_concatenate_genes.pl baseml/$count/*_tab.aln
	# 3.6. Remove unnecessary files 
	rm baseml/$count/*one_line*

	# 3.7. Get tree in PHYLIP format
	printf $num_seq" 1\n" > baseml/$count/$gene".tree"
	cat $i/"RAxML_bestTree"* >> baseml/$count/$gene".tree"

done

echo There are $count genes visited to generate summary statistic
echo There are $count genes visited to generate summary statistic >> out_logs/log_02_baseml_getdataformat.txt

## NOTE: This script might take a while to run if it is not run as a job
##       array in a cluster, so you might want to adapt this code accordingly.
##
## Example of bash script to be submited to the Apocrita HPC (QMUL) below.
##
## >> SCRIPT GOES FROM LINES 99-178. REMOVE "#> " AT THE BEGINNING OF EACH LINE, COPY
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
#> # 1. Get counter equal to $SGE_TASK_ID 
#> count=$( echo $SGE_TASK_ID ) 
#> 
#> # 2. Move to main dir, where this script saved in `00_Gene_filtering/scripts`
#> #    is run, i.e., `00_Gene_filtering`
#> curr_dir=$( pwd )
#> cd $curr_dir
#> cd ../../src
#> src=$( pwd )
#> cd $curr_dir
#> 
#> # 3. Save all the dirnames inside an array and then use 
#> #    $SGE_TASK_ID to access each gene name. Then, go back to 
#> #    $curr_dir
#> cd filtered_genes
#> filt_genes=(`ls`)
#> # Note: In order to access each name in this bash array, we
#> #       use "${filt_genes[${SGE_TASK_ID}]}"
#> i="${filt_genes[${SGE_TASK_ID}]}"
#> cd $curr_dir
#> 
#> # 4. Go through each filtered gene and prepare it for baseml 
#> #    analyses
#> 
#> # 4.0. Create dir
#> mkdir -p baseml/$count
#> 
#> # 4.1 Get gene name  
#> gene=$( echo $i | sed 's/..*\///')
#> # 4.2 Check num sequences
#> num_seq=$( grep '>' $i/*cds.aln.fasta | wc -l )
#> 
#> echo Parsing gene $count":" $gene >> out_logs/tmp_$SGE_TASK_ID".txt"
#> 
#> # 4.3 Get one line sequences
#> $src/one_line_fasta.pl $i/*cds.aln.fasta
#> mv $i/*one_line.fa baseml/$count
#> # 4.4. Get sequence next to header 
#> $src/00_get_seq_next_to_header.pl baseml/$count/*one_line*
#> # 4.5. Get alignment in PHYLIP format 
#> $src/01_concatenate_genes.pl baseml/$count/*_tab.aln
#> # 4.6. Remove unnecessary files 
#> rm baseml/$count/*one_line*
#> 
#> # 4.7. Get tree in PHYLIP format
#> printf $num_seq" 1\n" > baseml/$count/$gene".tree"
#> cat $i/"RAxML_bestTree"* >> baseml/$count/$gene".tree"
#> 
#> echo We have visited gene $count to generate summary statistic >> out_logs/tmp_$SGE_TASK_ID".txt"
#> 
#> # 5. Check if all genes have been visited and, if so, 
#> #    remove all the tmp files after concatenating them into a single file
#> num_genes=(`ls baseml/`)
#> visited_genes=$( echo ${#num_genes[@]} )
#> if [ $visited_genes -eq 15569 ]
#> then 
#>
#> for f in out_logs/tmp*
#> do 
#> cat $f >> out_logs/log_02_baseml_getdataformat.txt 
#> printf "\n" >> out_logs/log_02_baseml_getdataformat.txt 
#> rm $f 
#> done 
#>
#> fi
#>