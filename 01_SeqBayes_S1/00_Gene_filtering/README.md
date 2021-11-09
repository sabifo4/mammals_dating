# Filtering: dataset 1

## 1. Download and get a summary of mammal genes
In order to start the filtering step, you will need to download the 
`genes.zip` file. Please click [here](https://www.dropbox.com/s/voy9mvwn6amdw95/genes.zip?dl=0) to start the download.
Make sure you have enough space in your disk as its size is around 760 Mb.

Once the download is complete, please unzip it and save it in a directory 
called `genes`. You can manually unzip this by using software such as [`7-Zip`](https://www.7-zip.org/) or, 
if you prefer the command line, you can open a terminal within the `01_SeqBayes_S1` directory in the 
`mammals_dating` GitHub repository you have cloned and run the following:   

```sh
unzip genes.zip 
rm genes.zip
```

Make sure that the `genes` directory is saved under the `01_SeqBayes_S1/00_Gene_filtering/`
directory if you want to follow this tutorial with the commands detailed in each of the 
subsequent steps. This is because all relative paths assume that the main directory is `00_Gene_filtering` and 
the scripts used can be found within the directory `00_Gene_filtering/scripts`. Your file 
architecture should look like this:

```
mammals/
    |- 00_Data_collection/
    |- 01_SeqBayes_S1/
    |               |- 00_Gene_filtering/ 
    |               |                |- data_subsets/
    |               |                |- genes/        <-- LOCATION OF `GENES` DIRECTORY ONCE YOU UNZIP THE FILE
    |               |                |- out_data/
    |               |                |- out_logs/ 
    |               |                |- out_RData/
    |               |                |- scripts/      <-- WHERE WE HAVE SAVED OUR SCRIPTS
    |               |                |- README.md		   
    |               |- <other_dirs>
    |- <other_dirs>
```
**NOTE 1:** The directory `data_subsets` is needed for an independent analysis (i.e., Bayesian model selection analysis), and 
it is not relevant to this filtering analysis.   

**NOTE 2:** You will see that we provide you with the output files you should have at the end of this tutorial inside 
the directories `out_data`, `out_logs`, and `out_RData`. If you want to run this tutorial from scratch, you 
might want to save a copy of these output files if you want to later check you have been able to reproduce the tutorial 
and get the same results we got. Therefore, we recommend running the following commands from the terminal: 

```sh 
# Run from 01_SeqBayes_S1/00_Gene_filtering

# 1. Create a directory where you will keep the 
#    output files we generated for this study
#    (i.e., generate a backup)
mkdir copy_outfiles 
# 2. Move the directories where output files were 
#    generated to the directory previously created
mv out_* copy_outfiles/
# 3. Create directories for the output data you 
#    will generate when running this tutorial 
mkdir out_data out_logs out_Rdata
```

Once everything is ready, you can run the script [`00_Filtering_genes.sh`](scripts/00_Filtering_genes.sh) 
to obtain a comma separated file, `genes.csv`, which will be saved inside the `out_data` 
directory to which the following information is appended:   

   * Gene name   
   * Number of sequences included in the alignment   
   * Presence of _Mus musculus_ in the alignment   
   * Presence of _Homo sapiens_ in the alignment   
   * Number of codons the sequences in the alignment have   
   * Number of nucleotides the sequences in the alignment have   

If you are following the architecture provided in this GitHub repository
and the instructions provided in the last paragraph, you can 
just copy and paste the following command in your terminal. Otherwise, make the 
changes you need so you can run the next command accordingly: 

```sh
# Run from 01_SeqBayes_S1/00_Gene_filtering
./scripts/00_Filtering_genes.sh
```

## 2. Use summary statistics to filter genes 
The script [`00_Get_filtered_genes_in_dir.R`](scripts/00_Get_filtered_genes_in_dir.R)
finds the genes that meet the requirements needed to proceed with the Bayesian inference:   

   1. There must be at least 10 sequences in the alignments (18 genes out).   
   2. The alignments must contain at least 100 codons (323 genes out).   
   3. The presence of _Mus musculus_ and _Homo sapiens_ is required (2 genes out).   
   
In total, there are 335 unique genes that do not meet these three requirements,
thus removed in this first filtering step. The 15,569 filtered genes that 
pass the filtering step are saved in a directory called `filtered_genes`
(you can download the directory [here](https://www.dropbox.com/s/c0djydlrh2qipd8/filtered_genes.zip?dl=0) so you can check you get the same results we did).
A log file recording the order in which they have been parsed is also
output ([`log_01_copy_filtered_genes.txt`](out_logs/log_01_copy_filtered_genes.txt)).

Additionally, this script finds only those genes present in the 72 
mammal species that, at the same time, meet the three requirements stated 
above. This results in 648 genes, which are saved in a directory called `filtered_genes_all72sp`
(you can download the directory [here](https://www.dropbox.com/s/wsgfy58ufbjtvo5/filtered_genes_all72sp.zip?dl=0) so you can check you get the same results we did).
A log file is also generated ([`log_01_copy_filtered_all72sp_genes.txt`](out_logs/log_01_copy_filtered_all72sp_genes.txt)).

You can open the R script [`00_Get_filtered_genes_in_dir.R`](scripts/00_Get_filtered_genes_in_dir.R)
in RStudio to run all the steps detailed above. We have added extra information to this R script so it can be used 
as a tutorial.

## 3. Getting the alignments and trees in PAML format
The script [`01.1_Get_data_for_baseml.sh`](scripts/01.1_Get_data_for_baseml.sh)
is used to get the alignments and trees in the format required to be read by
`MCMCtree` and `BASEML` (`PAML v4.9h`, [Yang 2007](https://academic.oup.com/mbe/article/24/8/1586/1103731)).
It does the following:   
   1. Within a folder called `baseml`, it creates directories labelled
   as `1`, `2`, `n`; where `n` equals to the number of directories 
   with filtered genes. There, the corresponding gene and tree files will be saved.     
   2. It gets the gene name and number of sequences and saves them in a variable.   
   3. It gets the fasta alignments in one line format using the [`one_line_fasta.pl`](../../src/00_get_seq_next_to_header.pl) script available in the `src` directory.   
   4. It uses the script [`00_get_seq_next_to_header.pl`](../../src/00_get_seq_next_to_header.pl) to format the output from step 3 
   so the sequences are next to the species name tab separated (script available in the `src` directory).   
   5. It uses the script [`01_concatenate_genes.pl`](../../src/01_concatenate_genes.pl) to now generate an alignment in PHYLIP format.   
   6. It removes unnecessary files generated during steps 3-5 (script available in the `src` directory).   
   7. It gets the file containing the RAxML best-scoring ML tree and saves it in the corresponding `baseml/$num` 
   directory in PAML format.   
	   
In the end, each `baseml/n` directory has 4 files:   

   1. *aln: the alignment in PHYLIP format.   
   2. *fasta: the alignment in FASTA format.   
   3. *log.txt: a log file with the length of each sequence in the alignment.   
   4. *tree: the tree file in PAML format.   
   
Also, a log file [`log_02_baseml_getdataformat.txt`](out_logs/log_02_baseml_getdataformat.txt)
is generated once this script finishes. This script can be called as it follows:

```sh
# Run from 01_SeqBayes_S1/00_Gene_filtering
./scripts/01.1_Get_data_for_baseml.sh
```

You can download the zip file `baseml.zip` [here](https://www.dropbox.com/s/9ogr853wara72la/baseml.zip?dl=0)
so you can check that you have been able to reproduce these filtering steps. 

After that, the script [`partition_alignments.pl`](../../src/partition_alignments.pl)
(script available in the `src` directory) needs to be run for each of the gene 
alignments. The script [`01.1_Get_data_for_baseml_partitioned.sh`](scripts/01.1_Get_data_for_baseml_partitioned.sh)
can do this in a `for` loop. You can run it like this:

```sh
# Run from 01_SeqBayes_S1/00_Gene_filtering
./scripts/01.1_Get_data_for_baseml_partitioned.sh
```

This generates the `partitions12.3_*aln` and `partitions12_*aln` alignments within
the corresponding directories for each gene.

We will only be using the `partitions12_*aln` alignments as these are the ones which 
contain the 1st and 2nd codon positions (12CP) partition (partition that we are
interested in using for the subsequent analyses).

## 4. Get distance between mouse and human
The script [`01.2_Get_mouse_human_aln.sh`](scripts/01.2_Get_mouse_human_aln.sh)
extracts human and mouse sequences from `partitions12_*aln` and generates a file
readable by R called `*_mouse_human.aln` to later calculate the sequence distance.
A log file ([`log_03_get_mouse_human.txt`](out_logs/log_03_get_mouse_human.txt))
is also generated.

You can call this script as it follows: 

```sh
# Run from 01_SeqBayes_S1/00_Gene_filtering
./scripts/01.2_Get_mouse_human_aln.sh
```

Once the distance is computed in R in the next step, it will be 
used to order the genes from slow- to fast-evolving.

>>**NOTE**: In the end, once they are sorted out from slow- to fast-evolving, 
>>we will concatenate them and then divide them into 4 partitions.


## 5. Second filtering step 

### 5.1. Use relative branch lengths for second filtering step
The script [`01_Analysis_filtered_genes.R`](scripts/01_Analysis_filtered_genes.R)
is used to save the best-scoring ML gene trees of each 
alignment in a list. The trees are input in R using the `ape` package ([Paradis et al., 2004](https://academic.oup.com/bioinformatics/article/20/2/289/204981)) and are 
saved as class `phylo`. In order to avoid doing this every time the script 
is loaded in R, the RData file [`mammal.trees.RData`](out_RData/mammal.trees.RData) 
is created and saved in the `out_RData` directory.

Afterwards, the tree lengths are used to obtain the relative branch lengths, respectively 
for each alignment. If the largest branch length of any of the sequences of the alignment 
is larger than 60% of the total tree length, then the gene is removed. 

A matrix with this information (genes have not been removed yet) is created at this 
stage and saved as [`mammals_summary_matrix_filtstep1.RData`](out_RData/mammals_summary_matrix_filtstep1.RData). 
Also, a log file, [`log_04_R_genes_blength_longer_60pc.txt`](out_logs/log_04_R_genes_blenghth_longer_60pc.txt),
is created to keep track of the genes deleted.

After this step, 133 genes have at least one branch length that is 
60% larger than the total tree length, which are later removed from the gene pool 
and leaves 15,436 genes.

### 5.2. Calculate pairwise distance Human-Mouse for third filtering step 
Then, we extract the human and mouse sequences from the alignments with 1st and 2nd CPs 
(files generated in step 4). We read and convert them in R using the `ape::DNAbin`
function. In order to avoid doing this every time the script is loaded in R, the RData
file [`MHalns_DNAbin.RData`](out_RData/MHalns_DNAbin.RData)
is created and can then be easily loaded, if needed.
Note that the distance between the human and mouse sequences (H-M distance) is computed
for all the genes using the function `ape::dist.dna` for all the gene alignments and using
three implemented models:   
   1. `TN93`: We used this one as `BASEML` uses the TN93 formula to 
            calculate the pairwise distances for more complicated models than 
			JC69, K80, F81, and F84; e.g., HKY85.   
   2. `JC69`: We also calculated them using the JC69 model.   
   3. `raw`: This is a model implemented in the function just in case the 
             previous models return `Inf` or `NaN` values - which happens when 
			 the sequences are very different and consequently the evolutionary distances 
			 are labelled as "undefined" in this R function.   
			 
### 5.3. Gene ordering from slow- to fast-evolving
Afterwards, the resulting distances are used to
order the genes from slow- to fast-evolving. When we performed this filtering step, 
we noticed that there were 4 genes for which the distances were returned as `NaN` or were
larger than 0.75 for at least one of the models under which `ape::dist.dna` was ran.
Therefore, for the subsequent filtering, we decided to remove the following genes: `ENSG00000132185`, `ENSG00000204544`,
`ENSG00000120937`, and `ENSG00000236699`. Furthermore, when we plotted the tree length of
each gene alignment VS the corresponding largest branch length, we found and outlier 
(`ENSG00000176973`, see plot [here](out_RData/check_relblVStreelength.pdf)),
which was additionally removed (see plot after removing outlier [here](out_RData/check_relblVStreelength_nooutlier.pdf)).
We also log-transformed the x-axis (tree length) to ease its 
visualisation and plotted the resulting data (see the plot with the outlier [here](out_RData/check_logtreelengthVSlargestbl.pdf) 
and without the outlier [here](out_RData/check_logtree_relblVStreelength_nooutlier.pdf).

The resulting filtered and ordered genes from slow- to fast-evolving 
are saved in a directory called `filtered_genes_step2` (you can also download this directory 
[here](https://www.dropbox.com/s/dd31bddinnxvfak/filtered_genes_step2.zip?dl=0))
in phylip format (only alignment with 1st+2nd CP) together
with the corresponding tree file. In total, there are 15,431 filtered genes. We also 
generate a table with the summary statistics for these filtered genes,
[`filtering_sum_stats.csv`](out_RData/filtering_sum_stats.csv).

After that, those genes that have been filtered in this step are matched against
the genes that were found to be present in all 72 species (i.e., 645 genes)
in a previous step.
The resulting ordered genes from slow- to fast-evolving - and common in all 72 species -
are saved in a directory called `filtered_genes_step2_all72sp` (you can also download this 
directory [here](https://www.dropbox.com/s/4btyr36rwcjtcto/filtered_genes_step2_all72sp.zip?dl=0)).

### 5.4. Remove common nuclear genes that might be present in dataset 2
There is a chance that some nuclear genes present in this dataset 1 might overlap with the nuclear genes  
present in the dataset used in the second step of the sequential Bayesian dating analysis.
We need to remove them so there is no overlap between them.
For that purpose, we run the following commands: 

```sh 
## Run from `00_Gene_filtering`
./scripts/02.1_Remove_nucmit_genes_2ndstep.sh

# 1. Get dir names from log file and save in scripts/del_dirs.txt 
# 2. Run this for loop from `00_Gene_filtering`
del_dirs=$( cat scripts/del_dirs.txt | sed 's/\r//' )
for i in $del_dirs 
do 
name=`ls filtered_genes_step3/$i/*tree`
name_gene=$( echo $name | sed 's/..*\///' | sed 's/\.tree//' )
rm -r filtered_genes_step3/$i
printf "Dir "$i" with gene "$name_gene" has been removed\n"
printf "Dir "$i" with gene "$name_gene" has been removed\n" >> out_logs/log_11_removed_genes_done.txt
done 
```

A log file with the list of genes that are excluded can be found in 
[`log_10_removed_genes.txt`](out_logs/log_10_removed_genes.txt).

A total of 15,268 genes pass the filter (you can download the `filtered_genes_step3` directory [here](https://www.dropbox.com/s/677ip3hliybe9v1/filtered_genes_step3.zip?dl=0)). 
We also generated a table with the summary statistics for these filtered genes,
[`filtering_sum_stats_nonucgenes.csv`](out_RData/filtering_sum_stats_nonucgenes.csv).
 
Now, we do the same for the 645 genes that are shared across all 72 taxa:

```sh
## Run from `00_Gene_filtering`
./scripts/02.1_Remove_nucmit_genes_2ndstep_common72sp.sh

# 1. Get dir names from log file and save in scripts/del_dirs.txt 
# 2. Run this for loop from `00_Gene_filtering`
del_dirs_com72sp=$( cat scripts/del_dirs_com72sp.txt | sed 's/\r//' )
for i in $del_dirs_com72sp 
do 
name=`ls filtered_genes_step3_all72sp/$i/*tree`
name_gene=$( echo $name | sed 's/..*\///' | sed 's/\.tree//' )
rm -r filtered_genes_step3_all72sp/$i
printf "Dir "$i" with gene "$name_gene" has been removed\n"
printf "Dir "$i" with gene "$name_gene" has been removed\n" >> out_logs/log_14_removed_genes_common72sp_done.txt
done 
```

In total, 634 genes are present in all the 72 mammal taxa that are not present 
in the second dataset that will be used in the second step of the sequential 
Bayesian dating analysis (you can download the `filtered_genes_step3_all72sp` directory [here](https://www.dropbox.com/s/u9boyowl9djuvxx/filtered_genes_step3_all72sp.zip?dl=0)).

>>**NOTE: There are additional output files generated during the steps detailed above that might not have been extensively described here.**
>>**If you want to know more about how they were generated and their content, please follow the detailed comments**
>>**in the bash and R scripts mentioned above.**

## 6. Concatenate the genes in 4 partitions

## 6.1. Generate concatenated and partitioned alignments
Once the ordered genes have been filtered in the directory `filtered_genes_step3`, we need to concatenate 
the alignments with 1st+2nd CPs of each gene until we get 4 partitions ranging from slow- to 
fast-evolving genes.

We first need to run the script [`03_Concatenate_genes_separated_for_partition.sh`](scripts/03_Concatenate_genes_separated_for_partition.sh)
so we can get renamed individual files from slow- to fast-evolving (i.e., `geneX_ENSG..*.fasta`, where 
`x` is the number of the gene and the rest of the file name corresponds to the name 
of the gene), separated in four different folders to later run the partitioning 
pipeline. You can run this script as it follows: 

```sh
# Run from 01_SeqBayes_S1/00_Gene_filtering
./scripts/03_Concatenate_genes_separated_for_partition.sh
```

Now, we can use the `fasta-phylip-partitions` pipeline to generate the partitioned
alignments:

```sh
# Open a terminal within the 01_SeqBayes_S1/00_Gene_filtering/000_alignments
# directory and run the following code:

# Loop over all the "partX" directories, move there, and run the pipeline 
# within each directory to obtain individual concatenated alignments for  
# each of the partitions
for i in `seq 1 4` 
do 

cd part$i 

# 1. Create `species_names.txt`. In order to create a file listing the 72 species 
#    for which most of the molecular data can be found, we are going to use the 
#    alignment we generated for one of the genes within the `filtered_genes_step2_all72sp` 
#    directory. Therefore, we can run the following command which will extract 
#    the first column of this file - the one with the taxa names - and then output it 
#    in the main directory from which you should run this command (`000_alignments`) 
#    in a file called `species_names.txt`, which will be later used by the pipeline:
sed -n '2,73p' ../../filtered_genes_step2_all72sp/1/partitions12_*.aln | awk '{print $1}' > species_names.txt

# 2. Now, run the pipeline. If you wanted more information about how to install this pipeline 
#    and the arguments used,
#    please read the tutorial [here](https://github.com/sabifo4/fasta-phylip-partitions/blob/main/README.md)
#    NOTE: If you are running this code in a directory that you have synched to Dropbox or another 
#    cloud storage and you have issues, just move the folder out of the synched directory and run the 
#    command below again.
Run_tasks.sh . mammals partN

# 3. Move back to main directory for next iteration to proceed 
cd ..

done
```

Once the code above finishes to concatenate all the genes for each partition, we just need 
to concatenate the four partitions into a unique file that we will name `alignment_4parts.aln`: 

```sh
# Open a terminal within the 01_SeqBayes_S1/00_Gene_filtering/000_alignments
# directory and run the following code:
wd=$( pwd )
for i in `seq 1 4`
do

cd part$i/phylip_format/02_concatenated_alignments 
cat *aln >> $wd"/alignment_4parts.aln" 
cd $wd 

done  
```

And this is the end of the gene filtering step! Now we can proceed with 
the Hessian and gradient estimation for each of the partitions with `BASEML` so we
can calculate the approximate likelihood during the Bayesian dating analysis with `MCMCtree` ! 

If you have followed the steps mentioned above properly, you will be able to have 
both the concatenated and the partitioned alignments in the `000_alignments` directory. 

```
000_alignments/				  
          |- part[1-4]  # One directory for each partition (4 in total).
          |   |- phylip_format     # Files here are input/output files generated by
          |   |           |        # the pipeline fasta-phylip-partitions and the perl 
          |   |           |        # script count_missingdat_72sp.pl 
          |   |           |- 00_alignments_per_locus
          |   |           |       |- [0-9]*                     # One directory per locus
          |   |           |       |- log01_phylip_format.txt 
          |   |           |       |- log02_generate_tmp_aln.txt 
          |   |           |- 01_alignment_all_loci
          |   |           |       |- mammals_all_loci
          |   |           |- 02_concatenated_alignments
          |   |                 |- out_count_NA      # Output directory generated by count_missingdat_72sp.pl 
          |   |                 |- prepare_baseml    # Directory used to generate input files needed to run BASEML later
          |   |                 |- mammals_concat*   # Several output files by pipeline that start with this name
          |   |- geneX_ENSGX.fasta
          |   |- species_names.txt	   
          |     
          |- alignment_4parts.aln      # Concatenated alignment
          |- count_missingdat_72sp.pl  # Perl script to count missing taxa
          |- README.md
```

You can also download the zip file with the content of this directory [here](https://www.dropbox.com/s/mrvzzvd4o6qqyqk/000_alignments.zip?dl=0).
Make sure you have enough space as its size is ~2.6Gb.

### 6.1. Count missing data
You might have seen that there is a perl script called `count_missingdat_72sp.pl` 
saved in the `000_alignments` directory, as shown above.
You can also find it in this GitHub repository [here](scripts/count_missingdat_72sp.pl).
Note that, if you want to use it to reproduce our results, you will have to download the zip file with the
[`000_alignments`](https://www.dropbox.com/s/mrvzzvd4o6qqyqk/000_alignments.zip?dl=0)
as it uses the file architecture described above.
Alternatively, you can download the alignment for each partition and the concatenated 
alignment [here](https://www.dropbox.com/s/twp5uligudl46po/SeqBayesS1_mainT2hyp_alignments.zip?dl=0) and modify the code snippet we provide below accordingly so it can fit your file 
architecture (do this at your own risk, make sure you have the bioinformatics skills to do that).

If you have downloaded the zip file with the `000_alignments` directory, you will have the perl script
already there and you can run the following code to get the summary of the missing data for each 
partition:

```sh
# E.g. Run from 000_alignments
for i in `seq 1 4`
do 
cd part$i/phylip_format/02_concatenated_alignments
name_aln=`ls *aln`
printf "Parsing alignment part$i"/"$name_aln"... ...\n"
printf "Parsing alignment part$i"/"$name_aln"... ...\n" > log_count_missdat_part$i".txt"
perl ../../../count_missingdat_72sp.pl $name_aln >> log_count_missdat_part$i".txt"
# This directory is created by the PERL script, so 
# move log file there
mv log_count_missdat_part$i".txt" out_count_NA
cd ../../../
printf "\n\n"
done

# Get a summary!
for i in `seq 1 4`
do
printf "<< DIR part"$i" >>\n" >> part$i"_countNA.tsv"
sed -n '1,2p' part$i/phylip_format/02_concatenated_alignments/out_count_NA/*avgmissdata.txt >> part$i"_countNA.tsv"
sed -n '7,7p' part$i/phylip_format/02_concatenated_alignments/out_count_NA/*avgmissdata.txt >> part$i"_countNA.tsv"
sed -n '9,9p' part$i/phylip_format/02_concatenated_alignments/out_count_NA/*avgmissdata.txt >> part$i"_countNA.tsv"
printf "\n" >> part$i"_countNA.tsv"
done 
```

This will generate the output files that are saved in the `phylip_format/02_concatenated_alignments` directories 
within each `partX` directory inside the `000_alignments` directory. In addition, in this same subdirectory, 
you will see that there is a directory called `prepare_baseml`. There, you will be able 
to see a file called `out.txt`, which we have used to extract the information regarding the alignment 
length and the site patter counts for each partition (note that we explain in more detail how we ran `BASEML` 
to get this information in the section `2. Generate BASEML files` in the `BASEML` tutorial
[here](../01_BASEML/02_Hessian/README.md)).

We have summarised the content of each partition here:

| Data subset | No. taxa | No. orthologs | Alignment length (base pairs) | Site pattern counts | Missing data |
|-------------|----------|---------------|-------------------------------|---------------------|--------------|
| Partition 1 | 72       | 3,817         | 8,926,316                     | 3,613,711           | 60.17        |
| Partition 2 | 72       | 3,817         | 8,339,196                     | 2,941,508           | 50.78        |
| Partition 3 | 72       | 3,817         | 8,605,264                     | 2,416,624           | 48.49        |
| Partition 4 | 72       | 3,817         | 7,302,398                     | 1,521,894           | 45.92        |


**NOTE**: Now, you can move onto the next step to learn how to generate the files 
required by `BASEML` and `MCMCtree`. Nevertheless, if you want to know how we generated the 
data subsets for the Bayesian model selection analysis, you may want to read the 
next section!

# Extracting data subsets for Bayesian model selection analysis
In addition, we randomly subsampled 1, 10, 30, 50, 100, and 500 genes out of the 645 genes shared across all the 72 mammal taxa. 
These are saved in the corresponding directories named as `sampleXsp`, where `X` corresponds to the number of genes subsampled. 
Then, we increased the sampling to 1,000 and 10,000 by including all the 645 genes shared across 
all the 72 taxa and then adding the rest of needed genes by randomly sampling them from those saved in the
`filtered_genes_step2` directory (i.e., these other genes are not shared across all the 72 mammal taxa).
We have ensured that these last genes randomly added to the 1,000 and 10,000 data subsets are not 
repeated among the previously 645 genes that are shared across all the 72 mammal taxa.

After that, we randomly subsampled 10, 30, 100, 500, and 1,000 genes out of all the genes (not only those shared 
across all the 72 mammal taxa). This should provide us with a more randomised sample than the one described above.

**NOTE 1**: We used the genes filtered in step 2 even though there were chances we picked any of the 11
nuclear genes that are also found in the dataset that we use in the second step of the sequential
Bayesian dating analysis. This does not affect the results we get in this analysis because it is a complete
independent analysis from the clock-dating analysis.   

**NOTE 2**: The code for this purpose has been appended to the R script [`01_Analysis_filtered_genes.R`](scripts/01_Analysis_filtered_genes.R),
please check step 11 from line 717. If you run the code from this line, you will see that, in your working directory,
the following directories will be generated: 

```
00_Gene_filtering
          |- sample1sp 
          |- sample10sp 
          |- sample10sp_all 
          |- sample30sp 
          |- sample30sp_all 
          |- sample50sp 
          |- sample100sp 
          |- sample100sp_all 
          |- sample500sp 
          |- sample500sp_all 
          |- sample1000sp 
          |- sample1000sp_all 
          |- sample10000sp_all 
          |- <other_dirs>
```

After that, we decided to partition the alignments into 2 or 4 blocks from slow- to fast-evolving. We rewrote the 
bash script [`03_Concatenate_genes_separated_for_partition.sh`](scripts/03_Concatenate_genes_separated_for_partition.sh) 
so we could generate the number of partitions required. The code for this purpose can be found in the script 
[`04.1_Concatenate_genes_separated_for_partition_subsamples.sh`](scripts/04.1_Concatenate_genes_separated_for_partition_subsamples.sh)
and can be run as it follows:

```sh 
# Subsample 1sp
./scripts/04.1_Concatenate_genes_separated_for_partition_subsamples.sh 1 sample1sp 1

# Subsample 10sp
./scripts/04.1_Concatenate_genes_separated_for_partition_subsamples.sh 2 sample10sp 4
./scripts/04.1_Concatenate_genes_separated_for_partition_subsamples.sh 5 sample10sp 2

# Subsample 30sp
./scripts/04.1_Concatenate_genes_separated_for_partition_subsamples.sh 7 sample30sp 4
./scripts/04.1_Concatenate_genes_separated_for_partition_subsamples.sh 15 sample30sp 2

# Subsample 50sp
./scripts/04.1_Concatenate_genes_separated_for_partition_subsamples.sh 12 sample50sp 4
./scripts/04.1_Concatenate_genes_separated_for_partition_subsamples.sh 25 sample50sp 2

# Subsample 100sp
./scripts/04.1_Concatenate_genes_separated_for_partition_subsamples.sh 25 sample100sp 4
./scripts/04.1_Concatenate_genes_separated_for_partition_subsamples.sh 50 sample100sp 2

# Subsample 500sp
./scripts/04.1_Concatenate_genes_separated_for_partition_subsamples.sh 125 sample500sp 4
./scripts/04.1_Concatenate_genes_separated_for_partition_subsamples.sh 250 sample500sp 2

# Subsample 1000sp
./scripts/04.1_Concatenate_genes_separated_for_partition_subsamples.sh 250 sample1000sp 4
./scripts/04.1_Concatenate_genes_separated_for_partition_subsamples.sh 500 sample1000sp 2

# Subsample 10000sp
./scripts/04.1_Concatenate_genes_separated_for_partition_subsamples.sh 2500 sample10000sp 4
./scripts/04.1_Concatenate_genes_separated_for_partition_subsamples.sh 5000 sample10000sp 2

## Now do the same for those folders that contain the subsampled genes 
## from the whole dataset 

# Subsample 10sp
./scripts/04.1_Concatenate_genes_separated_for_partition_subsamples.sh 2 sample10sp_all 4
./scripts/04.1_Concatenate_genes_separated_for_partition_subsamples.sh 5 sample10sp_all 2

# Subsample 30sp
./scripts/04.1_Concatenate_genes_separated_for_partition_subsamples.sh 7 sample30sp_all 4
./scripts/04.1_Concatenate_genes_separated_for_partition_subsamples.sh 15 sample30sp_all 2

# Subsample 100sp
./scripts/04.1_Concatenate_genes_separated_for_partition_subsamples.sh 25 sample100sp_all 4
./scripts/04.1_Concatenate_genes_separated_for_partition_subsamples.sh 50 sample100sp_all 2

# Subsample 500sp
./scripts/04.1_Concatenate_genes_separated_for_partition_subsamples.sh 125 sample500sp_all 4
./scripts/04.1_Concatenate_genes_separated_for_partition_subsamples.sh 250 sample500sp_all 2

# Subsample 1000sp
./scripts/04.1_Concatenate_genes_separated_for_partition_subsamples.sh 250 sample1000sp_all 4
./scripts/04.1_Concatenate_genes_separated_for_partition_subsamples.sh 500 sample1000sp_all 2

```

Now, we can use another script to generate the alignments in PHYLIP format. 

>>**NOTE**: If you are running this code in a directory that you have synched to Dropbox or another 
>>cloud storage and you have issues, just move the folder out of the synched directory and run the 
>>code again.

```sh
# Subsample 1sp
./scripts/04.2_Generate_partitions_alignments_subsamples.sh sample1sp

# Subsample 10sp
./scripts/04.2_Generate_partitions_alignments_subsamples.sh sample10sp

# Subsample 30sp
./scripts/04.2_Generate_partitions_alignments_subsamples.sh sample30sp

# Subsample 50sp
./scripts/04.2_Generate_partitions_alignments_subsamples.sh sample50sp

# Subsample 100sp
./scripts/04.2_Generate_partitions_alignments_subsamples.sh sample100sp

# Subsample 500sp
./scripts/04.2_Generate_partitions_alignments_subsamples.sh sample500sp

# Subsample 1000sp
./scripts/04.2_Generate_partitions_alignments_subsamples.sh sample1000sp

## Now do the same for those folders that contain the subsampled genes 
## from the whole dataset 

# Subsample 10sp
./scripts/04.2_Generate_partitions_alignments_subsamples.sh sample10sp_all

# Subsample 30sp
./scripts/04.2_Generate_partitions_alignments_subsamples.sh sample30sp_all

# Subsample 100sp
./scripts/04.2_Generate_partitions_alignments_subsamples.sh sample100sp_all

# Subsample 500sp
./scripts/04.2_Generate_partitions_alignments_subsamples.sh sample500sp_all

# Subsample 1000sp
./scripts/04.2_Generate_partitions_alignments_subsamples.sh sample1000sp_all


```

Remember to generate also the concatenated files with 2 and 4 partitions!

You can download the zip file with the final directories output by R and the alignments generated 
[here](https://www.dropbox.com/s/9tjoavpk04ufzay/analyses_data_subsets.zip?dl=0).
The plots comparing the mean time estimates for each partition and each data subset,
together with the data and the scripts to plot this, can be found in
[`data_subsets`](data_subsets),
where you will have further details about how to download the results obtained with `MCMCtree` used 
to generate those plots.

----

# SUMMARY
If you have correctly followed this tutorial, you should see the following 
general file architecture under the `01_SeqBayes_S1` directory:

```
01_SeqBayes_S1/
     |- 00_Gene_filtering/ 
                  |- 000_alignments/				  
                  |        |- part1/  
                  |        |- part2/  
                  |        |- part3/  
                  |        |- part4/  
                  |        |- alignments_4parts.aln
                  |        
                  |- 00_alignments_subsamples/				  
                  |        |- sample1sp/  
                  |        |- ... 
                  |        |- sample10000sp			  
                  |        
                  |- baseml/      
                  |     |- X/                     * There are 15,569 directories here. The name of the directories,    
                  |       |- ENSG..*.aln          * "X" here, goes from 1 to 15,569. Inside each "X" directory,       
                  |       |- ENSG..*.fasta        * you can find the 7 files detailed here.
                  |       |- ENSG..*.log.txt      
                  |       |- ENSG..*.tree       
                  |       |- ENSG..*_mouse_human.aln     
                  |       |- partitions12.3_ENSG.*.aln   
                  |       |- partitions12_ENSG.*.aln       
                  |      
                  |- data_subsets/     
                  |     |- scripts/    * There is 1 R script
                  |     |- *pdf/       * There are 6 PDF files with plots
                  |     |- README.md/  * File explaining the steps to obtain the PDF files 
                  |      
                  |- filtered_genes/    
                  |     |- X/                           * There are 15,569 directories here. The name of the directories, "X" here,
                  |       |- ENSG..*.fas.sorted         * corresponds to the gene name (i.e., ENSG..**). Inside each      
                  |       |- ENSG..*.cds.aln.fasta      * "X" directory, you can find the 7 files detailed here.  
                  |       |- RAxML_bestTree.ENSG..*     
                  |       |- RAxML_bipartitions.ENSG..*      
                  |       |- RAxML_bipartitionsBranchLabels.ENSG..*    
                  |       |- RAxML_bootstrap.ENSG..*   
                  |       |- RAxMLinfo.ENSG..*   	             
                  |         	             
                  |- filtered_genes_all72sp 
                  |     |- X/                              * There are 648 directories here. The name of the directories, "X" here, 
                  |       |- ENSG..*.fas.sorted            * corresponds to the gene name (i.e., ENSG..**). Inside each      
                  |       |- ENSG..*.cds.aln.fasta         * "X" directory, you can find the 7 files detailed here.  
                  |       |- RAxML_bestTree.ENSG..*     
                  |       |- RAxML_bipartitions.ENSG..*      
                  |       |- RAxML_bipartitionsBranchLabels.ENSG..*    
                  |       |- RAxML_bootstrap.ENSG..*   
                  |       |- RAxMLinfo.ENSG..*   
                  |         	             
                  |- filtered_genes_step2/   
                  |     |- X/                           * There are 15,431 directories here. The name of the directories,    
                  |       |- ENSG..*.tree               * "X" here, goes from 1 to 15,431. Inside each "X"directory,     
                  |       |- partitions12_ENS..*.aln    * "X" directory, you can find the 2 files detailed here.    
                  |         					  
                  |- filtered_genes_step2_all72sp/    
                  |     |- X/                           * There are 645 directories here. The name of the directories,    
                  |       |- ENSG..*.tree               * "X" here, goes from 1 to 645. Inside each "X"directory,     
                  |       |- partitions12_ENS..*.aln    * "X" directory, you can find the 2 files detailed here.    
                  |         	              
                  |- filtered_genes_step3/   
                  |     |- X/                           * There are 15,268 directories here. The name of the directories,    
                  |       |- ENSG..*.tree               * "X" here, goes from 1 to 15,431. Inside each "X"directory,     
                  |       |- partitions12_ENS..*.aln    * "X" directory, you can find the 2 files detailed here.    
                  |         					  
                  |- filtered_genes_step3_all72sp/    
                  |     |- X/                           * There are 634 directories here. The name of the directories,    
                  |       |- ENSG..*.tree               * "X" here, goes from 1 to 645. Inside each "X"directory,     
                  |       |- partitions12_ENS..*.aln    * "X" directory, you can find the 2 files detailed here.  
                  |         	              
                  |- genes/  
                  |     |- X/                           * There are 15,904 directories here. The name of the directories, "X" here,
                  |       |- ENSG..*.fas.sorted         * corresponds to the gene name (i.e., ENSG..**). Inside each      
                  |       |- ENSG..*.cds.aln.fasta      * "X" directory, you can find the 7 files detailed here.  
                  |       |- RAxML_bestTree.ENSG..*     
                  |       |- RAxML_bipartitions.ENSG..*      
                  |       |- RAxML_bipartitionsBranchLabels.ENSG..*    
                  |       |- RAxML_bootstrap.ENSG..*   
                  |       |- RAxMLinfo.ENSG..*  
                  |         	              
                  |- out_data/
                  |       |- genes.csv				  
                  |         	              
                  |- out_logs/ 
                  |       |- log_*txt        * There are 15 log files.	
                  |         	              
                  |- out_RData/
                  |       |- check..*pdf     * There are 4 pdf files.	
                  |       |- *csv            * There are 4 csv files.	
                  |       |- ..*RData        * There are 5 RData files.	
                  |         		          
                  |- scripts/ 
                  |       |- 0..*sh   * There are 9 bash scripts.	
                  |       |- 0..*pl   * There is 1 perl script.	
                  |       |- 0..*R    * There are 2 R scripts.	
                  |       |- *txt     * There are 3 text files.
                  | 
                  |- README.md		   
```

If you need specific data for which we have not provided a link above to download, 
please send an e-mail to 
<a href="mailto:s.alvarez-carretero@ucl.ac.uk">s.alvarez-carretero@ucl.ac.uk</a> with your request!
