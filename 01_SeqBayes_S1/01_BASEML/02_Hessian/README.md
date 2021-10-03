# Running `BASEML`
Now that we have the partitioned alignments and the rooted trees for the 7 tree hypothesis with the 72
mammal taxa, we are ready to run `BASEML`.

The aim of this step is to calculate the Hessian and the gradient for each of the 4 partitioned 
alignments with genes ordered from slow- to fast-evolving under each of the 7 tree hypothesis. 
Both the Hessian and the gradient will be needed when we run `MCMCtree` to estimate 
the divergence times using the approximate likelihood. Note that both `MCMCtree` and `BASEML` are 
two software included in the `PAML` suite [Yang 2007](https://academic.oup.com/mbe/article/24/8/1586/1103731),
which you can download [here](http://abacus.gene.ucl.ac.uk/software/paml4.9j.tgz) (if you had issues with this link,
just click here and download the PAML. 
Note that the latest version is `v4.9j`, although we have inclued `PAML` v4.9h in the [src](https://github.com/sabifo4/mammals_dating/tree/main/src) 
directory because this is the one that was used in this study. If you are to run `BASEML` and `MCMCtree` 
to reproduce the results we have generated, make sure that you download and extract the [paml4.9h.zip](https://github.com/sabifo4/mammals_dating/tree/main/src/paml4.9h.zip) file,
install it on your PC, and export the path to the executable files.

All the tasks that are next described were carried out in the [Queen Mary's Apocrita HPC](https://zenodo.org/record/438045#.XhW1eUf7RPY).
The file architecture that was followed is the following:   

```
MAMMALS_Hessian/
             |- alignments/ 
             |         |- alignment_4parts.aln 
             |         |- mammals_concat_part1.aln 
             |         |- mammals_concat_part2.aln 
             |         |- mammals_concat_part3.aln 
             |         |- mammals_concat_part4.aln 
             |
             |- baseml_method1/
             |            |- 0X_..*/         # Note: There is one directory for each tree hypothesis 
             |                  |            # Each directory follows the same file architecture, only 
             |                  |            # the directory name changes 		 
             |                  |- p01/
             |                  |- p02/
             |                  |- p03/ 
             |                  |- p04/
             |                  |- in.BV*    # The in.BV file that will be later used by MCMCtree
             |                               # It contains the gradient and the Hessian of each partition 
             |
             |- pipelines_method1/ 
             |              |- 0X_..*tree/     # Note: There is one directory for each tree hypothesis
             |                  |              # Each directory follows the same file architecture, only 
             |                  |              # the directory name changes and the file name of the bash script
             |                  |              # saved in it			 
             |                  |            
             |                  |- pipeline_tree0[0-9].sh            
             |
             |- trees/ 
                  |- 0X_..*/      # Note: There is one directory for each tree hypothesis
                       |          # Each directory follows the same file architecture, only 
                       |          # the directory name changes and the file name corresponding
                       |          # tree hypothesis saved in it			 
                       |            
                       |- 72sp_MCMCtree..*.tree   
```

Keeping the large amount of data generated during the steps that will be next described 
is not possible in GitHub. 
Therefore, we zipped the `MAMMALS_Hessian` directory that we had in the Apocrita HPC once 
all the jobs were done and then uploaded it [here](https://www.dropbox.com/s/gxy13v0y7fjoefo/SeqBayesS1_BASEML_method1_6treehyp.zip?dl=0). You can just download the 
zipped file and go through each directory and pipelines used in this project. We recommend you do this once you finish to read the step-by-step tutorial provided below to check you have been able to reproduce our results.
>>**NOTE 1**: This zip file does not contain the `alignments` directory (see "IMPORTANT NOTE" below).

### IMPORTANT NOTE
The results for the main tree hypothesis (T2) were repeated as we found out that 11 nuclear genes out of the final 
15,268 genes included in data set 1 were also present in data set 2 (the data set that we were going to use in
the second step of the sequential Bayesian dating analysis).
While this does not affect the other 6 tree hypotheses, we had to make sure that there were no overlapping 
genes between the two data sets that are going to be used in the sequential Bayesian dating approach to
infer the divergence times, which affects the main tree hypothesis (T2) because this is the tree topology
used in this analysis. Consequently, we repeated the `BASEML` analysis with the alignment of data set 1 
without these 11 nuclear genes for the main tree because the posterior estimated times
obtained with the main tree were going to be then used to fit the skew-*t* distributions used as a prior
distributios in the second step of the sequential Bayesian approach. Note that:   

   * The alignments used to infer the divergence times for the other 6 tree hypotheses contain these 11 nuclear genes.
   These alignments can be downloaded from [here](https://www.dropbox.com/s/1z95mvjw0djgnqo/SeqBayesS1_6treehyp_alignments.zip?dl=0). These alignments were **not** used to infer the posterior divergence times with the 
   main tree hypothesis (T2) (we used the alignments in the [`000_alignments`](https://www.dropbox.com/s/mrvzzvd4o6qqyqk/000_alignments.zip?dl=0)
   that we describe in the previous gene filtering
   step), they were only used to infer the divergence times when using the other 6 tree hypotheses -- note that
   there is no conflict with these 6 tree hypotheses because they are not used in the subsequent steps and there is no need 
   to avoid overlapping genes.   
   * Last, please note that these alignments are not the ones provided in the `000_alignments` directory as describe above.
   The alignments in the `000_alignments` we describe in the gene filtering step and that can be downloaded [here](https://www.dropbox.com/s/mrvzzvd4o6qqyqk/000_alignments.zip?dl=0) are only
   used to infer the divergence times with the main tree hypothesis (T2).   
   
The zip file with the `BASEML` results with the main tree hytpohesis (T2) can be downloaded from
[here](https://www.dropbox.com/s/w6xnoleo4ssjalv/SeqBayesS1_BASEML_method1_T2hyp.zip?dl=0).
It follows the same file architecture as described above.

Below, we give you a detailed summary of how we carried out this analysis in the Apocrita HPC.

## 1. Loading data to the HPC and set working environment
First, we loaded both the concatenated and partitioned alignments to the HPC and saved them in the `alignments`
directory (or `00_alignments` when carrying out the analysis when using the T2 main tree). We have provided the 
links above to download these directories with the files you need to reproduce the analyses.

Then, we created one directory for each tree hypothesis where `BASEML` was going to be run 
for each of the partitions. We had 4 partitions, so we created four subdirectories within each of these 
directories as it follows:    

```
baseml_method1/
           |- 01_atlantogenata_scandentia_primates_tarver2016/   
           |     | 	        # Note: There is one directory for each tree hypothesis 
           |     |          # Each directory follows the same file architecture, only 
           |     |          # the directory name changes 		 
           |     |- p01/
           |     |- p02/
           |     |- p03/ 
           |     |- p04/
           |
           | .
           | .
           | .
           | 
           | 
           |- 07_afrotheria_tarver2016/   
                 | 	        # Note: There are 7 directories for each tree hypothesis 
                 |          # Each directory follows the same file architecture, only 
                 |          # the directory name changes 		 
                 |- p01/
                 |- p02/
                 |- p03/ 
                 |- p04/				
``` 

Then, we soft linked each partitioned alignment to its corresponding subdirectory 
(i.e., `p01` for first partition, `p02` for second partition, `p03` for third partition, 
and `p04` for fourth partition) in each of the seven directories for each tree hypothesis,
as well as the corresponding file with the tree hypothesis and the control file to run 
`MCMCtree`.

We used the next bash code for that purpose: 

```sh
# Bash code to be run in Apocrita HPC once you have logged in 
# an interactive session

# 1. Start interactive session 
qlogin 

# 2. Move to home directory and get pwd 
cd ~/MAMMALS_Hessian/
wd=$( pwd )

# 2. Run a `for` loop to go through so each of the directories 
#    and soft link each of the partitioned alignments, 
#    the corresponding tree hypothesis and the control file 
#    to run `MCMCtree`
cd baseml_method1 
for i in */
do
for j in `seq 1 4`
do

# Link tree and aln files. Copy ctl file and edit 
# accordingly
ln -s $wd/alignments/mammals_concat_part$j".aln" $i/p0$j
ln -s $wd/trees/$i/*tree $i/p0$j
cp $wd/control_file/mcmctree_MAMMALS.ctl $i/p0$j
cd $i/p0$j
tree_name=$( ls *tree )
aln_name=$( ls *aln )
sed -i -e 's/ALN/'${aln_name}'/' *ctl 
sed -i -e 's/TREE/'${tree_name}'/' *ctl

# Go back to baseml_dir for next iteration
cd $wd/baseml_method1

done

echo Visited directory $i

done
```

## 2. Generate `BASEML` files 
Now, within an interactive session in the HPC, we call `MCMCtree` within each subdirectory. In the control file used by
this dating software, we set the option
`usedata = 3`, but `MCMCtree` is not to be run until the task is completed. We will "kill" the process,
i.e., stop the `MCMCtree` run, once the files that are needed to run `BASEML` are generated in the correct format. 
This means that, once we see that `tmp000?..*` files generated (i.e., this is how the file names start for the 
alignment, tree, and control files we will need to later run `BASEML`)
and the message `XXXXX site patterns read, XXXXX sites` printed on the screen,
we can kill the job from the terminal (i.e., press `ctrl+C`). 

You just need to type the following code within each of the directories: 

```sh 
# 1. You should run this code within an interactive session 
qlogin 

# 2. Now move to each of the directories where the alignment, tree, and control 
#    files were previously saved and run the following command
#    NOTE: It assumes that you have exported the PATH to where you have installed 
#    the PAML suite (i.e., the `src` directory where the executable files of the 
#    software of this suite are).
mcmctree *ctl 

# 3. Once you see that the tmp000?.ctl, tmp000?.trees, tmp000?.txt files 
#    are generated and the message 
#    "XXXXX site patterns read, XXXXX sites" appears on the screen, 
#    kill the job by typing `ctrl+C`.
```

Then, you might want to remove output files we will not be using in the subsequent steps. In total, 
there are 8 output files that are not needed, which can be removed with the command below.
Note that the `SeedUsed` file should not be removed if you want to reproduce the exact same run again.
There are two ways of doing this, just pick the option you prefer from the ones given below: 

```sh 
# Option A: This command should be run inside each of the `p0[0-9]` 
# directories 
rm rst* out.* 2base.t rub tmp*out 

# Option B: If you want to run the command above in a `for` loop, 
# you can run the next code from the `baseml_method1` directory 
for i in */
do
for j in `seq 1 4`
do
cd $i/p0$j
rm rst* out.* 2base.t rub tmp*out
cd ../../
done
done
```

Last, you should update the control file. We need to change one of the options specified in 
the control file so we can speed up the computation of the Hessian and the 
gradient with a new approach implemented in `BASEML`. Pick one of the two options 
detailed below to carry out this task:

```sh
# Option A: This command should be run inside each of the `p0[0-9]` 
# directories 
sed -i 's/method\ \=\ 0/method\ \=\ 1/' tmp.*ctl

# Option B: If you want to run the command above in a `for` loop, 
# you can run the next code from the `baseml_method1` directory
for i in */
do
sed -i 's/method\ \=\ 0/method\ \=\ 1/' $i/p0*/tmp*ctl
done 
```

Once you have run the above code for each of the subdirectories inside the tree
directories, you can check that step 5 worked (i.e., the option `method` should 
be equal to 1, not to 0) by running the next command: 

```sh 
# Code to be run inside `baseml_method1` directory
for i in */
echo Visited tree directory $i
do grep 'method' $i/p0*/tmp*ctl
done
```

**NOTE**: In case this job could not be run in your HPC because of running out of memory, just increase 
the memory requirements of the interactive session, if allowed. Otherwise, you can always copy the 
`baseml_method1` directory to your desktop PC and follow the same instructions provided above.
Note that we assume that you have installed the PAML suite in your PC (as well as in your HPC!) as well as
exported the path to the executable files. Otherwise, change the commands to 
run `MCMCtree` in your PC accordingly. Note that you might need to wait ~2-5 minutes, depending on your RAM.
When these jobs were run on a desktop PC with i7-8750H, CPU 2.20GHz, 16GB RAM (OS: Windows 10, Linux subsystem used);
it took between 5-35 minutes, depending on the partition and alternative processes that were running on the PC 
at the same time. 

## 3. Preparing bash scripts to run `BASEML` 
When running jobs in the HPC that require more than some minutes, we cannot use the interactive session. 
Instead, we need to submit a bash script where all the tasks to be carried out 
need to be specified. 

For this particular study, we need to run `BASEML` for each partition under each tree hypothesis. 
For this kind of repeated tasks for which only the directory where they need 
to run changes as it uses different input files, we need to submit to the HPC what we call `job arrays`. 

```
=========================================================================
QUICK EXAMPLE OF A JOB ARRAY WITH 6 TREE HYPOTHESIS (EXCLUDING MAIN TREE T2) 
=========================================================================

In our case, we have 6 directories (one for each tree hypothesis) and, within
these 6 directories, 4 directories (one for each partition) where we have the
tmp* files previously generated. These tmp* files are the input files for 
BASEML. 

We just want to run the same command (i.e., running BASEML) within each of the 
directories that contain the tmp* files. Job arrays are like for loops. If we 
were to do this in a for loop, we would code it like this:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run from `baseml_method1` directory save working directory in var 
wd=$( pwd )
# Iterate through 6 directories with tree hypothesis 
for i in *
do 
# Iterate through 4 directories with partitions 
for j in `seq 1 4`
do 
# Move to the directory where tmp* files are 
cd $i/p0$j 
# Run baseml 
baseml *ctl
# Go back to working directory 
cd $wd 
done 
done 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To run this from the Apocrita HPC and to avoid time restrictions, we decided 
to have 6 bash scripts (one for each directory for the tree hypothesis) that would 
iterate through the 4 directories where the tmp* files are. In short, it is as if we 
now only needed this bit of the previous `for` loop: 

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Working directory is the tree directory, i.e., `0[0-9]..*`
wd=$( pwd )
# Iterate through 4 directories with partitions 
for j in `seq 1 4`
do 
# Move to the directory where tmp* files are 
cd p0$j 
# Run baseml 
baseml *ctl
# Go back to working directory 
cd $wd 
done 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When writing a job array, you do not need to write a `for` loop nor 
specify the variable that will change through each iteration (e.g., in the 
examples above, this would be `j` when going through `seq 1 4`). 
The job array has some options that are specified at the beginning of the script: 

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#!/bin/bash
#$ -cwd                     # Run the code from the current directory
#$ -V                       # Export environment to job
#$ -j y                     # Merge the standard output and standard error
#$ -l h_rt=240:00:00        # Limit each task to 10 days
#$ -l h_vmem=5G             # Request 5GB RAM
#$ -t 1-4
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The important option here is the last one, `#$ -t 1-4`. This option, `#$ -t`, says
that the code that is written below needs to be run 4 times. The variable used
for this "internal" `for` loop changes from HPC to HPC. The Apocrita HPC uses the
variable `SGE_TASK_ID`.
This means that, if it was written as a `for` loop, we would be writinf the following:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Working directory is the tree directory, i.e., `0[0-9]..*`
wd=$( pwd )
# Iterate through 4 directories with partitions 
for SGE_TASK_ID in `seq 1 4`
do 
# Move to the directory where tmp* files are 
cd p0$SGE_TASK_ID 
# Run baseml 
baseml *ctl
# Go back to working directory 
cd $wd 
done 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Nevertheless, this is not needed. When writing the job array, we would only need something like 
the following: 

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#!/bin/bash
#$ -cwd                     # Run the code from the current directory
#$ -V                       # Export environment to job
#$ -j y                     # Merge the standard output and standard error
#$ -l h_rt=240:00:00        # Limit each task to 10 days
#$ -l h_vmem=5G             # Request 5GB RAM
#$ -t 1-4

# Working directory is the tree directory, i.e., `0[0-9]..*`
wd=$( write_here_path_to_tree_directory )

# Move to the directory where tmp* files are 
cd $wd/p0$SGE_TASK_ID 

# Run baseml 
baseml *ctl
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As shown above, the only thing that we really need to specify is the path 
to the directory of the tree hypothesis (remember we had one bash script for each 
tree hypothesis, so this is the only bit that will change in the seven 
bash scripts). E.g., `wd=$( /data/scratch/MAMMALS_Hessian/baseml_method1/02_atlantogenata_tarver2016 )`.
Then, when this script is submitted to the HPC, it will be run in four different 
computers in parallel and, in each computer, the variable `$SGE_TASK_ID` will take one different 
number within the range specified in the option `#$ -t`, i.e., 1, 2, 3, or 4. Consequently, 
in each computer, BASEML will be run at the same time in parallel in each of the four 
directories where the corresponding partition was saved. 

To keep the output log files generated by each of the 6 bash scripts separated, 
we created the next file architecture: 

>>
pipelines_method1/
           |- 01_firsttree/   
           |     | 	        # Note: There are 6 directories for each tree hypothesis 
           |     |          # Each directory follows the same file architecture, only 
           |     |          # the directory name changes 		 
           |     |- pipeline_tree01.sh
           |
           | .
           | .
           | .
           | 
           | 
           |- 07_seventhtree/   
                 | 	       # Note: There are 6 directories for each tree hypothesis 
                 |         # Each directory follows the same file architecture, only 
                 |         # the directory name changes 		 
                 |- pipeline_tree07.sh
>>

This allows the output log files to be saved under the corresponding directory 
according to the tree hypothesis used, i.e., 4 log files, one for each partition 
run by the job array under each tree hypothesis.

```

## 4. Run `BASEML` 

Once all the bash scripts with the job arrays were written, they were submitted 
to the HPC. Depending on the partition, `BASEML` took between 6 to 20 hours, approximately, 
to calculate the Hessian and the gradient. 

### Preparing `in.BV` files for `MCMCtree`
When the jobs running `BASEML` finish, the software generates an output file called `rst2` where the 
gradient and Hessian are writen. This means that, for each partition under each tree hypothesis, we have 
an `rst2` file. When we are to later run `MCMCtree` for the divergence times estimation step, we need 
this information to enable the option that allows the likelihood to be approximated: what the documentation 
refers to the `in.BV` file.

To generate the `in.BV` file required by `MCMCtree`, we will concatenate the `rst2` files generated for each 
partition in a unique file. It is **very important** that we concatenate them in the correct order, i.e., from 
`p01` to `p04`, so each chunk of gradient and Hessian corresponds to the correct alignment that was used by `BASEML`
when they were calculated. 

The code that we have used is the following: 

```sh 
# Run from `baseml_method1` dir 
wd=$( pwd )
for i in *
do 

# 1. Move to directory for tree hypothesis i
cd $i
num_tree=$( echo $i | sed 's/\_..*//' ) 
for j in `seq 1 4`
do 
# 2. Copy the content of rst2 files inside 
#    each p0$j directory and save in 
#    file that will be later used by MCMCtree 
#    as in.BV. Add new lines between chunk of 
#    gradient + Hessian from one partiton to another
echo Adding gradient and Hessian of part$j under tree hypothesis $i
cat p0$j/rst2 >> in.BV.$num_tree".p1-4"
printf "\n\n" >> in.BV.$num_tree".p1-4"
done 

# 3. Print spaces for screen and go back to main dir 
cd $wd
printf "\n\n"

done 
```

The code above generates within each directory for each tree hypothesis a file called 
`in.BV.0X.p1-4`, where `X` is the index number used to identify the tree hypotheses.
This file will be later used by `MCMCtree` in the next step!
