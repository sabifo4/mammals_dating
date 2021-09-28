# AFROTHERIA - Filtering alignment
The initial files that you can find in this directory are the following:

```
afrotheria
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
this file [here](https://www.dropbox.com/s/ymtr3guph5s6yzr/SeqBayesS2_filtaln_afrotheria.zip?dl=0).
To start the filtering step, you should have the following files here, which 
you can obtain once you unzip the file provided above: 

```
afrotheria 
   |- filter_aln
         |- checked_aln                        
         |     |- alignment_nt3cp.phylip      # Alignment with 3CP (found in `original_aln_tree`)
         |     |- RAxML_bestTree.BS_ML_GTRCAT # Unfiltered tree (found in `original_aln_tree`)
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
file. The messages printed out by this script are the following:

```
ORDER   Sirenia  has  3  taxa. Species  trichechus_manatus_latirostris  has  2  genes that are not shared
by any of the other taxa. It should be deleted!
```
>> ACTION: trichechus_manatus_latirostris was removed.

```
FAMILY   Trichechidae  has  2  taxa. Species  trichechus_manatus  has  35  genes that are not shared by any
of the other taxa. It should be deleted!
FAMILY   Trichechidae  has  2  taxa. Species  trichechus_manatus_latirostris  has  2  genes that are not
shared by any of the other taxa. It should be deleted!
```
>> ACTION: there are only 2 taxa for family Trichechidae. This means that subspecies
trichechus_manatus_latirostris needs to be removed as it has less genes than the sequence for 
at the species level. 

```
SUBFAMILY   Amblysominae  has  9  taxa. Species  amblysomus_iris  has  1  genes that are not shared by any
of the other taxa. It should be deleted!
```
>> ACTION: amblysomus_iris was removed.

```
GENUS   Oryzorictes  has  2  taxa. Species  oryzorictes_hova  has  2  genes that are not shared by any of
the other taxa. It should be deleted!
GENUS   Oryzorictes  has  2  taxa. Species  oryzorictes_talpoides  has  4  genes that are not shared by any
of the other taxa. It should be deleted!
```
>> ACTION: there are only 2 taxa for genus Oryzorictes. We keep the sequence with more genes,
oryzorictes_talpoides, and remove the other: oryzorictes_hova.

```
GENUS   Amblysomus  has  6  taxa. Species  amblysomus_iris  has  1  genes that are not shared by any of
the other taxa. It should be deleted!
```
>> ACTION: same finding as in subfamily: amblysomus_iris already removed.

```
GENUS   Trichechus  has  2  taxa. Species  trichechus_manatus  has  35  genes that are not shared by an
of the other taxa. It should be deleted!
GENUS   Trichechus  has  2  taxa. Species  trichechus_manatus_latirostris  has  2  genes that are not
shared by any of the other taxa. It should be deleted!
```
>> ACTION: same finding as in family: trichechus_manatus_latirostris already removed.


**NOTE**
Log files will be found in this same directory with the name `log_taxaNOTIN_<level>.txt`,
being `level` the one that has been checked (i.e., family, genus, order, or subfamily). In addition, 
you can also generate the `levels.checked.RData`, which you can use to explore the taxonomical 
levels. 

## 2. Check alignment with 3CP partition
The alignment with the first and second codon positions (12CP), file `alignment.phylip`
(you will find it once you unzip the file which link is provided above),
had previously undergone the previous filtering in 2012. Nevertheless, the alignment 
with the third codon positions (3CP), file `alignment_nt3cp.phylip` (you will find it once 
you unzip the file which link provided above, inside a directory called `checked_aln` directory), was not yet filtered. 

To check if both alignment files (`alignment.phylip` and `alignment_nt3cp.phylip`) had the
same taxa, we ran the following check:

```sh
# Run the next code from `afrotheria/filter_aln/checked_aln`

# 1. Check they have the same info
grep -o '[a-z].* ' alignment_nt3cp.phylip > names_3nt.txt
grep -o '[a-z].* ' ../alignment.phylip > ../names_aln.txt
diff names_3nt.txt ../names_aln.txt
## 22c22
##< elephantulus_rozeti
##---
##> petrosaltator_rozeti

## >NOTE: Apparently one species was renamed. Use python script to fix this in `names.txt`

# 2. Fix species names in `names.txt`. We used the RAxML file, the `names.txt` file, and 
#    the `lineages.txt` that had been generated when filtering the alignment 12CP.
#    The ML tree file is now saved inside `filter_aln/checked_aln`, from where the 
#    python script will be run. The `names.txt` and the `lineage.txt` files are in
#    the `filter_aln` dir, one up. Note that we keep a copy of original tree and
#    alignment files before filtering in dir `checked_aln/original_aln_tree`
mkdir original_aln_tree 
cp alignment_nt3cp.phylip original_aln_tree/
cp RAxML_bestTree.BS_ML_GTRCAT original_aln_tree/
cat > _species_name_change.txt << EOF
elephantulus_rozeti,petrosaltator_rozeti
EOF

mkdir screen_logs
python ../../../../../../src/fix_species_names.py _species_name_change.txt RAxML_bestTree.BS_ML_GTRCAT alignment_nt3cp.phylip ../names.txt ../lineage.txt > screen_logs/01_log1_fixspecies.txt

## > NOTE: The `names.txt` file was already fixed. Screen log shows that `elephantulus_rozeti`
##         is not found, so then `names.txt` (output file by python script) 
##         inside `checked_aln` does not need to be updated: 
##
##    60 names found
##    Error: elephantulus_rozeti not found!
##    60 names after 0 replacements

## 3. Move unncessary output files to 00_filt1 and rename output alignment
mkdir 00_filt1
mv _species_name_change.txt lineage.txt lineage.txt.1.bak 00_filt1/
mv names_3nt.txt 00_filt1/
rm alignment_nt3cp.phylip
mv alignment_nt3cp_out.phylip alignment_nt3cp.phylip

# 4. Check now if they have the same info 
grep -o '[a-z].* ' alignment_nt3cp.phylip > names_3nt.txt
diff names_3nt.txt ../names_aln.txt
## 22c22
## < petrosaltator_rozeti
## ---
## > petrosaltator_rozeti

# Apparently there are not 3 spaces in the `alignment_nt3cp.phylip`, so we can manually change this
# and check again 
sed -i 's/petrosaltator\_rozeti /petrosaltator\_rozeti   /g' alignment_nt3cp.phylip
grep -o '[a-z].* ' alignment_nt3cp.phylip > names_3nt.txt
diff names_3nt.txt ../names_aln.txt # OK! 
```

## 3. Generate alignment files
Now that both files with the 12CP-alignment (`alignment.phylip`) and 3CP-alignment (`alignment_nt3cp.phylip`) 
have the same taxa names and have had the same filtering procedure carried out,
we can proceed to concatenate the alignment with the third 
codon positions to the alignment with 12CPs using the R script
[`Concatenate_seqs_for_MCMCtree.R`](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/01_alignments/Concatenate_seqs_for_MCMCtree.R)
so we can have the concatenated alignment. 

Instructions to follow:   

   * Open the RScript [`Concatenate_seqs_for_MCMCtree.R`](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/01_alignments/Concatenate_seqs_for_MCMCtree.R) 
   in RStudio and change line 24 so it is `subt    <- "afrotheria"`. Now, we can run it from RStudio. This script
   will generate a concatenated alignment file with all partitions, as well as 
   one alignment file for each individual partitions, inside a new dir called `00_mammal_alns/afrotheria`
   inside [`00_Data_filtering/01_alignments/`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/01_alignments). 
   Log and RData files will be saved inside 
   [`00_Data_filtering/01_alignments/Rout`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/01_alignments/Rout/log_concatenation). 
      > NOTE 1: Updated `partitions.txt` file inside `00_mammal_alns/afrotheria` generated.   
	  > NOTE 2: Output file called `check_appends.txt` to check that 3nt partition has been appended
	    to the right line of the alignment can be found in `Rout/Rdata". 
		If "0" in second column, something wrong! So far, everything seems OK!   
   * Paths have been automatically set according to current file architecture in `mammals` dir, 
   do not change paths in the Rscript! 

Note that the alignment files and the R objects are not provided in the GitHub repository as they 
are too large. We provide a link to the final files generated once the curation step finishes (see 
last section of this file).

## 4. Second filtering step
Check if any further species should be removed:

```sh
# Run from `00_data_curation/afrotheria/filter_aln/checked_aln` directory 
grep -o '[a-z]*_[a-z]*_[a-z]*' RAxML_bestTree.BS_ML_GTRCAT | sed 's/\_/\t/g' | awk '$2==$3{print $1,$2,$3}'  > ../subsp_check.txt

# Count 
num=$( grep -o '[a-z]*_[a-z]*_[a-z]*' RAxML_bestTree.BS_ML_GTRCAT | sed 's/\_/\t/g' | awk '$2==$3{print $1,$2,$3}' | wc -l )
printf "There are $num subspecies with species\n"
## There are 2 subspecies with species
```

The output names are saved in a file called `subsp_check.txt` (inside `filter_aln`)
to further explore it:

```sh
# Run from `00_data_curation/afrotheria/filter_aln/checked_aln` directory 
input="../subsp_check.txt"
while IFS= read -r line 
do

subsp=$( echo $line | sed 's/ /\_/g' )
num_subsp=$( grep -ocP ''${subsp}'' ../../../genes.txt )
sp=$( echo $line | awk '{print $1,$2}' | sed 's/ /\_/g' )
num_sp=$( grep -ocP '\b'${sp}'\b' ../../../genes.txt )
printf "Num. sequences $sp : $num_sp\n"
printf "Num. sequences $subsp : $num_subsp\n"
printf "\n"

printf "Num. sequences $sp : $num_sp\n" >> ../sum_subspVSsp.txt
printf "Num. sequences $subsp : $num_subsp\n" >> ../sum_subspVSsp.txt
printf "\n" >> ../sum_subspVSsp.txt

done < $input
```

Output:

> **Remove subspecies, less genes**   
   Num. sequences elephas_maximus : 52      
   Num. sequences **elephas_maximus_maximus** : 2

> **Remove subspecies, less genes**   
   Num. sequences macroscelides_proboscideus : 40   
   Num. sequences **macroscelides_proboscideus_proboscideus** : 2
 
As this filtering was not applied before we generate the alignment, we are going to use the already
concatenated alignment with filtered `3nt` and 
`1st+2nd` CPs and remove these taxa. For that purpose, we are going to run the next commands: 

```sh
# Run from `filter_aln/checked_aln` directory 

# 1. Copy filtered alignment in `00_Data_filtering/01_alignments/00_mammal_alns/afrotheria/` 
#    called `afrotheria.aln` in the `filter_aln/checked_aln` directory. Remember that this is 
#    the concatenated filtered alignment that we previously generated. Once there, change the 
#    file name: 
cp ../../../../01_alignments/00_mammal_alns/afrotheria/afrotheria.aln .
mv afrotheria.aln alignment.phylip

# 2. Remove species 

# 2.1. Create input file for python script `prune_alignment.py`. This text file 
#    contains a list of the species that need to be removed 
cat > _species_remove.txt << EOF
elephas_maximus_maximus
macroscelides_proboscideus_proboscideus
EOF

# 2.2. Run four python scripts that will:
#    a) Convert the PHYLIP alignment into FASTA format 
#    b) Use the alignment in FASTA format and the previous text file `_species_remove.txt` 
#       to remove the lines containing the listed species in the text file.
#    c) Remove the gaps from the alignment. The alignment was done on AAs, so everything will 
#       come in lines of three and the sites should be homologous. Therefore, when one species is removed,
#       if there were any insertions in that removed species, so it will be in the whole column
#       Therefore, this gap column is removed and nothing should alter the alignment.
#    d) Convert fasta into phylip.
python ../../../../../../src/phy2fasta.py alignment.phylip > alignment.fasta
python ../../../../../../src/prune_alignment.py --remove _species_remove.txt alignment.fasta > alignment_pruned.fasta
python ../../../../../../src/remove_gap_columns.py alignment_pruned.fasta > alignment_pruned_nogaps.fasta
python ../../../../../../src/fasta2phy.py alignment_pruned_nogaps.fasta > alignment_pruned_nogaps.phylip

# 2.3. Now prune the tree.
#      Note that tool used before 2018 was called using 
#      this command: <path_to_binary>/nw_prune -f <aln_to_prune> <file_with_names> > outfile
#      Currently, option `-f` seems to be deprecated and it does not work. Therefore, we had to put 
#      all the names of taxa to be removed in a variable and pass it as the second argument instead 
#      of using a text file 
sp2rm=$( cat _species_remove.txt )
../../../../../../src/newick-utils-1.6/src/nw_prune RAxML_bestTree.BS_ML_GTRCAT $sp2rm > tmp.tree #COOL!

##  Move unnecessary output files to 01_filt2. Rename tree and alignment files.
mkdir 01_filt2 
mv *fasta alignment.phylip _species_remove.txt 01_filt2/
mv RAxML_bestTree.BS_ML_GTRCAT 01_filt2/RAxML_bestTree.BS_ML_GTRCAT_unpruned
mv tmp.tree RAxML_bestTree.BS_ML_GTRCAT
mv alignment_pruned_nogaps.phylip alignment.phylip

# 3. Remove duplicates 

# 3.1. Find duplicates
grep -oP '\w+_duplicate' RAxML_bestTree.BS_ML_GTRCAT # No duplicates!

# 4. Generate a list of names with the taxa in the updated filtered alignment 
grep -o '[a-z].* ' alignment.phylip > names_3nt_filt.txt
```

## 5. Data partitioning 
Now that we have the concatenated and filtered alignment ready, we need to generate the filtered
partitioned alignments by running the R script 
[`Partition_seqs_for_MCMCtree_after_filtering`](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/01_alignments/Partition_seqs_for_MCMCtree_after_filtering.R).
As we had already created partitioned alignments 
before the subspecies check had been applied, we will need to rearrange the files output for
this data subset in the
[`00_Data_filtering/01_alignments/`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/01_alignments)
subdirectories:

```sh
# Run from 00_Data_filtering/01_alignments 

# Update alignments dir 
cd 00_mammal_alns/afrotheria 
mkdir afrotheria_old 
mv *aln *txt afrotheria_old 

# Update Rout dir 
cd ../../Rout/log_concatenation/
mkdir afrotheria_old 
mv *afrotheria*txt afrotheria_old
cd ../Rdata/
mkdir afrotheria_old
mv *afrotheria*txt *afrotheria*RData afrotheria_old
```

Now, run the Rscript mentioned above following the next instructions:   

   * Open the RScript [`Partition_seqs_for_MCMCtree_after_filtering.R`](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/01_alignments/Partition_seqs_for_MCMCtree_after_filtering.R) 
   in RStudio and change line 24 so it is `subt    <- "afrotheria"` and uncomment line 30. Now, we can run it from RStudio. 
   This script will generate a concatenated alignment file with all partitions, as well as 
   one alignment file for each individual partitions, inside the directory `00_mammal_alns/afrotheria/`
   inside [`01_alignments`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/01_alignments).
   Log files and Rdata can be found [here](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/01_alignments/Rout/log_concatenation). 
      > NOTE: Updated `partitions.txt` file inside `00_mammal_alns/afrotheria/`.   
   * Paths have been automatically set according to current file architecture in `mammals` dir, 
   do not change paths in the Rscript! 
   
The final alignments generated at the end of this step can be downloaded from 
[here](https://www.dropbox.com/s/rdatnuzurpj3o8s/SeqBayesS2_Raln_afrotheria.zip?dl=0).
They should be saved here if the same file architecture as the one set in the R scripts 
is to be used: `00_Data_filtering/01_alignments/00_mammals_alns/afrotheria`.
