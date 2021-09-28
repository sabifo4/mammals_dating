#!/bin/bash
#$ -cwd                     # Run the code from the current directory
#$ -V                       # Export environment to job
#$ -j y                     # Merge the standard output and standard error
#$ -l h_rt=240:00:00        # Limit each task to 10 days
#$ -l h_vmem=10G            # Request 10GB RAM
#$ -l node_type=nxv         # Request fastest nodes
#$ -t 2-3

# Create log file
exec 3>&1> >(while read line; do echo "$line" >> log.mammals.mafft.ctenohystrica.aln$SGE_TASK_ID.txt; done;) 2>&1

#############################################
## START SETTING FILES/FOLDER ARCHITECTURE ##
#############################################

echo The analyses will be carried out in the directory $SGE_TASK_ID
printf "\n"

# Move to the home directory
cd $SGE_TASK_ID

# Get aln file name
aln_name=$( ls cteno*fasta )
out_name=$( echo $aln_name | sed 's/\.fasta/\_out\.fasta/' )
new_seq=$( ls taxa*fasta )

# Run PRANK
printf "\nRunning MAFFT ...\n"
mafft --add $new_seq $aln_name >$out_name

printf "\n"
echo MAFFT finished"!"
printf "\n"
