# MARSUPIALIA - Filtering alignment
The initial files that you can find in this directory are the following:

```
marsupialia
   |- filter_aln
         |- lineage.txt 
         |- marsupialia_taxonomy_check.csv 
         |- names.txt 
         |- parse_lineage.R 
         |- partitions.txt 
         |- README.md
         |- summary.html 
```
		 
The output files during the filtering step that is described below as well as the input 
files needed have been zipped in a file as they are very large. You can download 
this file [here](https://www.dropbox.com/s/slp7ufwk7759q7p/SeqBayesS2_filtaln_marsupialia.zip?dl=0).
To start the filtering step, you should have the following files here, which 
you can obtain once you unzip the file provided above: 

```
marsupialia 
   |- filter_aln
         |- checked_aln                        
         |     |- alignment_nt3cp.phylip      # Alignment with 3CP (found in `original_aln_tree`)
         |     |- RAxML_bestTree.BS_ML_GTRCAT # Unfiltered tree (found in `original_aln_tree`)
         |     
         |- alignment.phylip # Alignment with 12CP
         |- lineage.txt 
         |- marsupialia_taxonomy_check.csv 
         |- names.txt 
         |- parse_lineage.R 
         |- partitions.txt 
         |- README.md
         |- summary.html 
```

Note that you will see more files within the zipped file provided above. These are output 
files that you will generate if you follow all the steps below. Feel free to keep them so you can 
make sure that you reproduce the same results that we show here. 


## 1. Taxonomic filtering 
Within the [`filter_aln`](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/marsupialia/filter_aln)
directory, you will find different files as detailed above. 
Specifically, the R script [`parse_lineage.R`](02_SeqBayes_S2/00_Data_filtering/00_data_curation/marsupialia/filter_aln/parse_lineage.R)
was written to carry out a first taxonomic filtering. Note that this 
R script will run if you have the same file architecture in this GitHub repository (i.e., it uses 
a function within the R script [`Filter_lineages.R`](src/Filter_lineages.R)
in the `src` directory and the [`genes.txt`](02_SeqBayes_S2/00_Data_filtering/00_data_curation/genes.txt)
file. The messages printed out by this script are the following:

```
ORDER   Diprotodontia  has  114  taxa. Species  phalanger_maculatus_duplicate  has  0  genes that are not shared by any of the other taxa. 
It should be deleted!
```
>> ACTION: phalanger_maculatus_duplicate was removed.

```
ORDER   Dasyuromorphia  has  78  taxa. Species  antechinus_naso_duplicate  has  0  genes that are not shared by any
of the other taxa. It should be deleted!
```

>> ACTION: antechinus_naso_duplicate was removed.

```
FAMILY   Phalangeridae  has  18  taxa. Species  phalanger_maculatus_duplicate  has  0  genes that are not shared
by any of the other taxa. It should be deleted!
```
>> ACTION: same finding as in order:  phalanger_maculatus_duplicate already removed.

```
FAMILY   Dasyuridae  has  77  taxa. Species  antechinus_naso_duplicate  has  0  genes that are not shared by
any of the other taxa. It should be deleted!
```
>> ACTION: same finding as in order: antechinus_naso_duplicate already removed.

```
GENUS   Lasiorhinus  has  2  taxa. Species  lasiorhinus_krefftii  has  11  genes that are not shared by any of
the other taxa. It should be deleted!
GENUS   Lasiorhinus  has  2  taxa. Species  lasiorhinus_latifrons  has  3  genes that are not shared by any
of the other taxa. It should be deleted!
```
>> ACTION: this genus has only two taxa, we decided to keep both.   

```
GENUS   Phalanger  has  8  taxa. Species  phalanger_maculatus_duplicate  has  0  genes that are not shared by any
of the other taxa. It should be deleted!
```
>> ACTION: same finding as in order: phalanger_maculatus_duplicate already removed.   


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

The first manual checks suggest some species names have to be changed as well as some 
sequences are to be removed. The justifications can be found in the file 
[`marsupialia_taxonomy_check.csv`](02_SeqBayes_S2/00_Data_filtering/00_data_curation/euarchonta/filter_aln/marsupialia_taxonomy_check.csv). 
The filtering you will see in the commands below was already 
carried out for `alignment.phylip`, but has had to be done again for `alignment_nt3cp.phylip` so both 
alignment files are at the same filtering level. The code is the following:


```sh
# Run the next code from `marsupialia/filter_aln/checked_aln`

# 0. Names to be changed according to `marsupialia_taxonomy_check.csv`:
# 
# marmosa_murina_macrotarsus,marmosa_macrotarsus
# marmosa_murina_waterhousei,marmosa_waterhousei
# micoureus_paraguayanus,marmosa_micoureus_paraguayana
# micoureus_demerarae,marmosa_micoureus_demerarae
# micoureus_regina,marmosa_micoureus_regina
# micoureus_constantiae,marmosa_micoureus_constantiae
# antechinus_naso,murexia_naso
# phascomurexia_naso,murexia_naso
# phalanger_maculatus,spilocuscus_maculatus
# marmosa_robinsoni_simonsi,marmosa_simonsi
# marmosa_robinsoni_isthmica,marmosa_isthmica

## >NOTE
##  Problems observed: antechinus_naso and phascomurexia_naso are renamed the same:
##  `murexia_naso`. Also, `phalanger_maculatus` is renamed into 
##  `spilocuscus_maculatus`, which already exists.
##  >WHAT DO WE DO?
##           Remove any duplicates (those with fewest genes):
##
##           antechinus_naso,murexia_naso : remove this one (3 mit genes)
##           phascomurexia_naso,murexia_naso_duplicate : rename this one to murexia_naso_duplicate (all mit, 1 nuc gene)
##           phalanger_maculatus,spilocuscus_maculatus_duplicate : remove this one (1 gene)
##
## Now we can proceed with the filtering!

# 1. Fix species names in `names.txt`. We used the RAxML file, the `names.txt` file, and 
#    the `lineages.txt` that had been generated when filtering the alignment 12CP.
#    The ML tree file is now saved inside `filter_aln/checked_aln`, from where the 
#    python script will be run. The `names.txt` and the `lineage.txt` files are in
#    the `filter_aln` dir, one up. Note that we keep a copy of original tree and
#    alignment files before filtering in dir `checked_aln/original_aln_tree`
mkdir original_aln_tree 
cp alignment_nt3cp.phylip original_aln_tree/
cp RAxML_bestTree.BS_ML_GTRCAT original_aln_tree/
cat > _species_name_change.txt << EOF
marmosa_murina_macrotarsus,marmosa_macrotarsus
marmosa_murina_waterhousei,marmosa_waterhousei
micoureus_paraguayanus,marmosa_micoureus_paraguayana
micoureus_demerarae,marmosa_micoureus_demerarae
micoureus_regina,marmosa_micoureus_regina
micoureus_constantiae,marmosa_micoureus_constantiae
phascomurexia_naso,murexia_naso
marmosa_robinsoni_simonsi,marmosa_simonsi
marmosa_robinsoni_isthmica,marmosa_isthmica
EOF

mkdir screen_logs
python ../../../../../../src/fix_species_names.py _species_name_change.txt RAxML_bestTree.BS_ML_GTRCAT alignment_nt3cp.phylip ../names.txt ../lineage.txt > screen_logs/01_log1_fixspecies.txt

## Move unncessary output files to 00_filt1 and rename output alignment
mkdir 00_filt1
mv _species_name_change.txt 00_filt1/
rm alignment_nt3cp.phylip
mv alignment_nt3cp_out.phylip alignment_nt3cp.phylip

# 2. Remove species 

# 2.1. Create input file for python script `prune_alignment.py`. This text file 
#    contains a list of the species that need to be removed 
cat > _species_remove.txt << EOF
antechinus_naso
phalanger_maculatus
EOF

# 2.2. Run set of four python scripts
python ../../../../../../src/phy2fasta.py alignment_nt3cp.phylip > alignment_nt3cp.fasta
python ../../../../../../src/prune_alignment.py --remove _species_remove.txt alignment_nt3cp.fasta > alignment_nt3cp_pruned.fasta
python ../../../../../../src/remove_gap_columns.py alignment_nt3cp_pruned.fasta > alignment_nt3cp_pruned_nogaps.fasta
python ../../../../../../src/fasta2phy.py alignment_nt3cp_pruned_nogaps.fasta > alignment_nt3cp_pruned_nogaps.phylip

# 2.3. Now prune the tree
sp2rm=$( cat _species_remove.txt )
../../../../../../src/newick-utils-1.6/src/nw_prune RAxML_bestTree.BS_ML_GTRCAT $sp2rm > tmp.tree #COOL!

##  Move unnecessary output files to 01_filt2. Rename tree and alignment files.
mkdir 01_filt2 
mv *fasta alignment_nt3cp.phylip _species_remove.txt 01_filt2/
mv RAxML_bestTree.BS_ML_GTRCAT 01_filt2/RAxML_bestTree.BS_ML_GTRCAT_unpruned
mv tmp.tree RAxML_bestTree.BS_ML_GTRCAT
mv alignment_nt3cp_pruned_nogaps.phylip alignment_nt3cp.phylip

# 3. Are there duplicates?
grep -oP '\w+_duplicate' RAxML_bestTree.BS_ML_GTRCAT # No! 

# 4. I generated a list of names to check that they are the same 
grep -o '[a-z].* ' alignment_nt3cp.phylip > names_3nt.txt
grep -o '[a-z].* ' ../alignment.phylip > ../names_aln.txt
diff names_3nt.txt ../names_aln.txt # OK!
```

Now that both files with the 12CP-alignment (`alignment.phylip`) and 3CP-alignment (`alignment_nt3cp.phylip`) 
have the same taxa names and have had the same filtering procedure carried out,
we can proceed to concatenate the alignment with the third 
codon positions to the alignment with 12CPs using the R script
[`Concatenate_seqs_for_MCMCtree.R`](02_SeqBayes_S2/00_Data_filtering/01_alignments/Concatenate_seqs_for_MCMCtree.R)
so we can have the concatenated alignment. 

Instructions to follow:      

   * Open the RScript [`Concatenate_seqs_for_MCMCtree.R`](02_SeqBayes_S2/00_Data_filtering/01_alignments/Concatenate_seqs_for_MCMCtree.R) 
   in RStudio and change line 24 so it is `subt    <- "marsupialia"`. Now, we can run it from RStudio. This script
   will generate a concatenated alignment file with all partitions, as well as 
   one alignment file for each individual partitions, inside a new dir called `00_mammal_alns/marsupialia`
   inside [`00_Data_filtering/01_alignments/`](/02_SeqBayes_S2/00_Data_filtering/01_alignments). 
   Log and RData files will be saved inside 
   [`00_Data_filtering/01_alignments/Rout`](/02_SeqBayes_S2/00_Data_filtering/01_alignments/Rout/log_concatenation). 
      > NOTE 1: Updated `partitions.txt` file inside `00_mammal_alns/marsupialia` generated.   
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

```sh
# Run from `filter_aln/checked_aln` directory 
grep -o '[a-z]*_[a-z]*_[a-z]*' RAxML_bestTree.BS_ML_GTRCAT | sed 's/\_/\t/g' | awk '$2==$3{print $1,$2,$3}' | awk '$1!="rl"{print $1,$2,$3}' > ../subsp_check.txt

# Count 
num=$( grep -o '[a-z]*_[a-z]*_[a-z]*' RAxML_bestTree.BS_ML_GTRCAT | sed 's/\_/\t/g' | awk '$2==$3{print $1,$2,$3}' | awk '$1!="rl"{print $1,$2,$3}' | wc -l )
printf "There are $num subspecies with species\n"
## There are 1 subspecies with species
```
The output names have been saved in a file called `subsp_check.txt` to further explore it:

```sh
# Run from `00_data_curation/marsupialia/filter_aln/checked_aln` directory 
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

> Species had already been removed   
   ~~Num. sequences sminthopsis_leucopus : 3~~   
   Num. sequences sminthopsis_leucopus_leucopus : 3
   
No further changes required!

The final alignments generated at the end of this step can be downloaded from 
[here](https://www.dropbox.com/s/s3ujh36ei1lb34k/SeqBayesS2_Raln_marsupialia.zip?dl=0).
They should be saved here if the same file architecture as the one set in the R scripts 
is to be used: `00_Data_filtering/01_alignments/00_mammals_alns/marsupialia`.