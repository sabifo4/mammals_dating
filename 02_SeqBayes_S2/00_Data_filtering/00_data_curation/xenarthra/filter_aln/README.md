# XENARTHRA - Filtering alignment
The initial files that you can find in this directory are the following:

```
xenarthra
   |- filter_aln
         |- lineage.txt 
         |- names.txt 
         |- parse_lineage.R 
         |- partitions.txt 
         |- README.md
         |- summary.html 
         |- taxonomy_check.csv 
```
		 
The output files during the filtering step that is described below as well as the input 
files needed have been zipped in a file as they are very large. You can download 
this file [here](https://www.dropbox.com/s/qi2lvst5r61xaty/SeqBayesS2_filtaln_xenarthra.zip?dl=0).
To start the filtering step, you should have the following files here, which 
you can obtain once you unzip the file provided above: 

```
xenarthra 
   |- filter_aln
         |- checked_aln                        
         |     |- alignment_nt3cp.phylip      # Alignment with 3CP 
         |     |- RAxML_bestTree.BS_ML_GTRCAT
         |     
         |- alignment.phylip # Alignment with 12CP
         |- lineage.txt 
         |- names.txt 
         |- parse_lineage.R 
         |- partitions.txt 
         |- README.md
         |- summary.html 
         |- taxonomy_check.csv 
```

Note that you will see more files within the zipped file provided above. These are output 
files that you will generate if you follow all the steps below. Feel free to keep them so you can 
make sure that you reproduce the same results that we show here. 

## 1. Taxonomic filtering 
Within the [`filter_aln`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/afrotheria/filter_aln)
directory, you will find different files as detailed above. 
Specifically, the R script [`parse_lineage.R`](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/afrotheria/filter_aln/parse_lineage.R)
was written to carry out a first taxonomic filtering. Note that this 
R script will run if you have the same file architecture in this GitHub repository (i.e., it uses 
a function within the R script [`Filter_lineages.R`](https://github.com/sabifo4/mammals_dating/blob/main/src/Filter_lineages.R)
in the `src` directory and the [`genes.txt`](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/genes.txt)
file. There are no messages printed, which means that there are no taxa 
to be removed.

**NOTE**
Log files will be found in this same directory with the name `log_taxaNOTIN_<level>.txt`,
being `level` the one that has been checked (i.e., family, genus, order, or subfamily). In addition, 
you can also generate the `levels.checked.RData`, which you can use to explore the taxonomical 
levels. 

## 2. Check alignment with 3CP partition
The alignment with the first and second codon positions (12CP), file `alignment.phylip`
(you will find it once you unzip the file which link is provided above),
had previously undergone the previous check in 2012 (nothing to remove).
Nevertheless, the alignment 
with the third codon positions (3CP), file `alignment_nt3cp.phylip` (you will find it once 
you unzip the file which link provided above, inside a directory called `checked_aln` directory),
this check had not been done yet. 

To check if both alignment files (`alignment.phylip` and `alignment_nt3cp.phylip`) had the
same taxa, we ran the following check:

```sh
# Run the next code from `xenarthra/filter_aln/checked_aln`

# 1. Check they have the same info
grep -o '[a-z].* ' alignment_nt3cp.phylip > names_3nt.txt
grep -o '[a-z].* ' ../alignment.phylip > ../names_aln.txt
diff names_3nt.txt ../names_aln.txt # OK! 
```

Now that both files with the 12CP-alignment (`alignment.phylip`) and 3CP-alignment (`alignment_nt3cp.phylip`) 
have the same taxa names and have had the same filtering procedure carried out,
we can proceed to concatenate the alignment with the third 
codon positions to the alignment with 12CPs using the R script
[`Concatenate_seqs_for_MCMCtree.R`](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/01_alignments/Concatenate_seqs_for_MCMCtree.R)
so we can have the concatenated alignment. 

Instructions to follow:  

   * Open the RScript [`Concatenate_seqs_for_MCMCtree.R`](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/01_alignments/Concatenate_seqs_for_MCMCtree.R) 
   in RStudio and change line 24 so it is `subt    <- "xenarthra"`. Now, we can run it from RStudio. This script
   will generate a concatenated alignment file with all partitions, as well as 
   one alignment file for each individual partitions, inside a new dir called `00_mammal_alns/xenarthra`
   inside [`00_Data_filtering/01_alignments/`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/01_alignments). 
   Log and RData files will be saved inside 
   [`00_Data_filtering/01_alignments/Rout`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/01_alignments/Rout/log_concatenation). 
      > NOTE 1: Updated `partitions.txt` file inside `00_mammal_alns/xenarthra` generated.   
	  > NOTE 2: Output file called `check_appends.txt` to check that 3nt partition has been appended
	    to the right line of the alignment can be found in `Rout/Rdata". 
		If "0" in second column, sth wrong! So far, everything seems OK!   
   * Paths have been automatically set according to current file architecture in `mammals` dir, 
   do not change paths in the Rscript! 

Note that the alignment files and the R objects are not provided in the GitHub repository as they 
are too large. We provide a link to the final files generated once the curation step finishes (see 
last section of this file).

## 3. Second filtering step
Check if any further species should be removed:

```{sh}
# Run from `filter_aln/checked_aln` directory 
grep -o '[a-z]*_[a-z]*_[a-z]*' RAxML_bestTree.BS_ML_GTRCAT | sed 's/\_/\t/g' | awk '$2==$3{print $1,$2,$3}' 
# Blank output, OK!
```

No further changes required!

The final alignments generated at the end of this step can be downloaded from 
[here](https://www.dropbox.com/s/j50gk080m1juzca/SeqBayesS2_Raln_xenarthra.zip?dl=0).
They should be saved here if the same file architecture as the one set in the R scripts 
is to be used: `00_Data_filtering/01_alignments/00_mammals_alns/xenarthra`.