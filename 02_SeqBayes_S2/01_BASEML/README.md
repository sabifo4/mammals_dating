# BASEML
All the details given below refer to the contents zipped in
[this file](https://www.dropbox.com/s/ajzlkhdnqdclgrh/SeqBayesS2_BASEML.zip?dl=0).
You may download it so you can go through all the steps easier.

## 1. Prepare architecture and prepare baseml files
The initial main architecture of the directories in the Apocrita HPC looks like this:

```
<data_subset>/
          |- [1-6]/
               |- prepare_baseml/
			        |- mcmctree_preparebaseml
						 |- mcmctree_preparebaseml.ctl
```

The control file is edited for each data subset so the correct paths for the directories where 
the alignment and tree files were saved is set. This procedure avoids linking file or copying files,
thus saving memory and file space in the cluster.
Directories from 1 to 6 will contain the gradient and Hessian for the following alignments:   

   * `1`: concatenated alignment 5 partitions.   
   * `2`: mt_12cp partition.   
   * `3`: mt_3cp partition.   
   * `4`: mt_rna partition.   
   * `5`: nt_12cp partition.   
   * `6`: nt_3cp partition.   
   
## 2. Generate directory for 5-partitions alignment
Once `BASEML` has run for all the data subsets (you will find the corresponding control files in each 
directory in the zip file provided above), we want to concatenate the `rst2` (i.e., the output 
files by `BASEML` that contain the gradient and the Hessian information for each data subset) generated 
for each data subset so we can generate the alignment with 5 individual partitions (5 blocks), the one that 
we will use for the Bayesian dating analysis with `MCMCtree`.
This is the content you will find for each data subset under directory `7` for each data subset.

```sh
# Run from <name_data_subset> directory
for i in `seq 2 6`
do

printf "Appending rst2 from dir "$i"...\n"
cat $i/rst2 >> 7/rst2 
printf "\n\n" >> 7/rst2

done
```
