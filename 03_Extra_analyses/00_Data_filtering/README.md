# 1. Process data subsets 
As an extra analysis, we assessed the difference in divergence time estimates when using 
nuclear and mitocohondrial data alone for the 72 taxa that are present in the main tree 
(i.e., the tree used in the first step of the Bayesian sequential-subtree analysis).

For that purpose, we screened the alignments generated for each data subset used in the second 
step of the Bayesian analysis (alignments generated with the data available for 4.7K taxa)
and extracted the mitochondrial (1st+2nd CPs) and nuclear (1st+2nd CPs) 
sequences of the 72 taxa that were common with the main tree (72 taxa and ~15K genomic data).

For that purpose, we used the data structure we had generated to run `MCMCtree`, which 
you can download 
[here](). 
If you want to reproduce our analysis, please download the file which link is provided above 
and save it in this directory under the name `alignments`. Then, 
run the following code:

```sh
# Run from `alignments` 
# We will extract the mit 1st+2nd (2)
# and nuc 1st+2n (5).
# First one, then another

# MIT-12CP
home_dir=$( pwd )
cp ../scripts/taxa_72sp.txt ../scripts/tmp_taxa_72sp.txt
for i in */2
do
cd $i
printf "<< DIR "$i" >>\n"
perl ../../../scripts/Extract_72sp.pl *aln ../../../scripts/tmp_taxa_72sp.txt | tee ../../../logs/log_extract_mit12cp.txt
# Update the "tmp_taxa_72sp.txt" file in its 
# file destination with only those taxa left to explore in 
# next subtrees. This avoids having duplicated taxa as some 
# subtrees have taxa present in others as we needed them for 
# specific calibrations in subsequent MCMCtree analyses
mv tmp_taxa_72sp.txt ../../../scripts
cd $home_dir 
done
cd $home_dir
# Check that 72sp have been collected!
grep ^[0-9]* -o */2/out*/*aln | sed 's/..*\://' | awk '{s+=$1} END {print s}'
# Remove tmp file 
rm ../scripts/tmp_taxa_72sp.txt

# NUC-12CP
home_dir=$( pwd )
cp ../scripts/taxa_72sp.txt ../scripts/tmp_taxa_72sp.txt
for i in */5
do
cd $i
printf "<< DIR "$i" >>\n"
perl ../../../scripts/Extract_72sp.pl *aln ../../../scripts/tmp_taxa_72sp.txt | tee ../../../logs/log_extract_nuc12cp.txt
mv tmp_taxa_72sp.txt ../../../scripts
cd $home_dir 
done
cd $home_dir
# Check that 72sp have been collected!
grep ^[0-9]* -o */5/out*/*aln | sed 's/..*\://' | awk '{s+=$1} END {print s}'
# Remove tmp file 
rm ../scripts/tmp_taxa_72sp.txt
```

Now, we just need to concatenate the output FASTA files for the 72sp as we will need 
to realign these sequences later:

```sh
# Run from `alignments`
# MIT-12CP 
for i in */2/out*/*fasta
do
cat $i >> mit12cp_unaln72sp.fasta 
done 

# NUC-12CP 
for i in */5/out*/*fasta
do
cat $i >> nuc12cp_unaln72sp.fasta
done
```

# 2. Generate alignments 
The next step is to generate the alignments for each partition that has been extracted 
in the previous step. We will use `MAFFT` for this purpose. The next snippet 
was used to generate three different alignments (i.e., nuc12CP and mit12CP):

```sh
# Run from Apocrita HPC 
# Note that $aln_name and $out_name 
# are variables updated within the bash 
# script submitted to the HPC for each 
# of the FAST alignemnts for each 
# partition
mafft --auto $aln_name > $out_name
```

The resulting alignments in FASTA format for each partition were downloaded and saved in directory 
`00_partition_alignments`. The file architecture is the following: 

```
# There are two directories labelled from
# `1` to `2` and, inside each of those, you 
# can find the input FASTA file for MAFFT and  
# the corresponding output alignments in FAST 
# format
00_partition_alignments 
     |- [1-2]             
	     |- *fasta 
	     |- *out.fasta
```

Now, we need to prepare the FASTA files in the correct format 
so the script `FASTAtoPHYL.pl` can be used to convert them into 
PHYLIP format. For that purpose, we do the following:

```sh
## Run from the directory `00_partition_alignments`
for i in `seq 1 2`
do

cd $i 
fasta=$( echo *_out.fasta )
printf "Parsing "$fasta"...\n"
one_line_fasta.pl $fasta
grep '>' $fasta | sed 's/>//' > species.txt
mkdir out_mafft 
mv *fasta out_mafft
fa=$( ls *fa )
name=$( echo $fa | sed 's/\_out..*//' | sed 's/unaln//' )
mv $fa $name".fasta"
cd ../

done
```

Last, we convert them into PHYLIP format:

```
## Run from the directory `00_partition_alignments`
for i in `seq 1 2`
do 

cd $i
num=$( grep '>' *fasta | wc -l )
len=$( sed -n '2,2p' *fasta | sed 's/\r//' | sed 's/\n//' | wc --m )
perl ../../scripts/FASTAtoPHYL.pl *fasta $num $len 
cd ..

done
```

Before generating the tree files, we will make sure that the final alignments do not contain the `~`
for some taxa
(i.e., `tursiops_truncatus~`, `cavia_aperea~`, and `felis_catus~`). We can run the following loop
for that purpose: 
 
```
## Run from the directory `00_partition_alignments`
for i in `seq 1 2`
do
sed -i 's/\~//g' $i/*aln
done
```

[Here]()
you can download the `00_partition_alignments` directory with all the files (input/output) that have been 
used/generated in the steps detailed in this section.

# 3. Generate uncalibrated and calibrated trees 
Before running `MCMCtree`, we need to prepare the tree files so we can fix the tree topology with the
corresponding node calibrations as well as the uncalibrated tree to run `BASEML` and compute the Hessian
and the gradient that will be used by `MCMCtree` when using the approximate likelihood.

For the calibrated trees, we will use the fossil calibrations that were used to run the divergence-time
estimation analysis with the first data set (72 taxa and ~15K loci):

```sh
# Run from `trees` 
cp original/72sp_atlantogenata_tarver2016_MCMCtree_calib.tree 00_fossil_calibs/72sp_fossilcalibs.tree
```

Note that we will use the main tree topology that we used in the first step of the Bayesian
sequential-subtree dating analysis. 
Before we start, however, we need to change the following taxa names in the calibrated and uncalibrated 
trees (the latter will be used to run `BASEML`):

```sh
# The following files are used and now saved in `original` directory
for i in 0*/*tree
do
sed -i 's/carlito\_syrichta/tarsius\_syrichta/' $i
sed -i 's/cricetulus\_griseus\_chok1gshd/cricetulus\_griseus/' $i
sed -i 's/heterocephalus\_glaber\_female/heterocephalus\_glaber/' $i
sed -i 's/notamacropus\_eugenii/macropus\_eugenii/' $i
sed -i 's/saimiri\_boliviensis\_boliviensis/saimiri\_boliviensis/' $i 
done
```

[Here]()
you can find the tree files that will be used to run `MCMCtree` and `BASEML`.

Now, we are ready to run `BASEML` to compute the Hessian and the gradient needed to use the 
approximate likelihood calculation implemented in `MCMCtree` to speed up divergence times inference!

# 4. Run `BASEML` 
We will compute the Hessian for each of the individual partitions. For each partition, we first generate 
the input files:

```sh
# Run from each corresponding `00_partition_alignments/Hessian` directory.
# Kill the job as soon as input files are generated 
mcmctree mcmctree_preparebaseml.ctl | tee logfile.txt
```

Then, we transfer the input files `tmp000*` to the cluster 
and run `BASEML`. 

Please note that we prepared the `BASEML` files in the PC and ran `BASEML` in Apocrita HPC.
[Here]()
you can download the output files generated by `BASEML` for each partition alignment 
as well as the input files we used to run the software.


# 5. Run `MCMCtree`
Once the `in.BV` files are generated, we can run `MCMCtree`.

