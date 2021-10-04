# RODENTIA CTENOHYSTRICA - Filtering alignment
The initial files that you can find in this directory are the following:

```
rodentia_ctenohystrica
   |- filter_aln
         |- checked_aln                        
         |     |- taxonomical_check
         |           |- rodentia_ctenohystrica.png
         |           |- rodentia_ctenohystrica_FigTree
         |           |- RAxML_bestTree.BS_ML_GTRCAT
         |     
         |- extra_filtering   <-- You will read more about it at the end of this guideline                     
         |    
         |- lineage.txt 
         |- names.txt    
         |- parse_lineage.R 
         |- partitions.txt 
         |- README.md
         |- rodentia_ctenohystrica_taxonomy_check.csv
         |- summary.html 
```
		 
The output files during the filtering step that is described below as well as the input 
files needed have been zipped in a file as they are very large. You can download 
this file [here](https://www.dropbox.com/s/u3nhv7ra64wtv64/SeqBayesS2_filtaln_ctenohystrica.zip?dl=0).
To start the filtering step, you should have the following files arranged in the file 
architecture detailed above (you can obtain the files once you unzip the file
provided in the link above): 

```
rodentia_ctenohystrica 
   |- filter_aln
         |- checked_aln                        
         |     |- unfiltered_aln
         |     |    |- alignment_nt3cp.phylip         # Alignment with 3CP 
         |     |- RAxML_bestTree.BS_ML_GTRCAT         # Unfiltered tree (found in `original_tree`)
         |     |- taxonomical_check                   # Visual checks to evaluate 
         |          |- rodentia_ctenohystrica.png     # dubious taxa placement
         |          |- rodentia_ctenohystrica_FigTree
         |          |- RAxML_bestTree.BS_ML_GTRCAT    
         |     
         |- extra_filtering   <-- You will read more about it at the end of this guideline                     
         |          
         |- alignment.phylip # Alignment with 12CP
         |- lineage.txt 
         |- names.txt 
         |- parse_lineage.R 
         |- partitions.txt
         |- README.md
         |- rodentia_ctenohystrica_taxonomy_check.csv
         |- summary.html 
```

Note that you will see more files within the zipped file provided above. These are output 
files that you will generate if you follow all the steps below. Feel free to keep them so you can 
make sure that you reproduce the same results that we show here. 

## 1. Taxonomic filtering 
Within the `filter_aln`
directory, you will find different files as detailed above. 
Specifically, the R script [`parse_lineage.R`](parse_lineage.R)
was written to carry out a first taxonomic filtering. Note that this 
R script will run if you have the same file architecture in this GitHub repository (i.e., it uses 
a function within the R script [`Filter_lineages.R`](../../../../../src/Filter_lineages.R)
in the `src` directory and the [`genes.txt`](../../genes.txt)
file. There are several messages printed out by this script as all the rodents are 
parsed at the same time. You might have seen that, due to the large number of rodent species, 
we have had to divide them into different data subsets: "rodentia_squirrel",
"rodentia_ctenohystrica" (this one), "rodentia_subtree1", and "rodentia_subtree2". Below, you will find 
the messages printed out that refer to this data subset, "rodentia_ctenohystrica":

```
GENUS   Galea  has  4  taxa. Species  galea_monasteriensis  has  1  genes that are not
shared by any of the other taxa. It should be deleted!
```
>> ACTION: galea_monasteriensis removed.

```
GENUS   Microcavia  has  2  taxa. Species  microcavia_australis  has  2  genes that are not
shared by any of the other taxa. It should be deleted!
GENUS   Microcavia  has  2  taxa. Species  microcavia_niata  has  1  genes that are not
shared by any of the other taxa. It should be deleted!
```
>> ACTION: microcavia_australis & microcavia_niata merged, the label is "microcavia".

**NOTE**
Log files will be found in this same directory with the name `log_taxaNOTIN_<level>.txt`,
being `level` the one that has been checked (i.e., family, genus, order, or subfamily). In addition, 
you can also generate the `levels.checked.RData`, which you can use to explore the taxonomical 
levels. Remember that you will see more information about the other rodent taxa that are 
present in the other three data subsets.

# 2. First checks before applying filtering
First, we checked if there were any further species that should be removed:

```sh
# Run from `checked_aln` directory 
grep -o '[a-z]*_[a-z]*_[a-z]*' RAxML_bestTree.BS_ML_GTRCAT | sed 's/\_/\t/g' | awk '$2==$3{print $1,$2,$3}'  > subsp_check.txt

# Count 
num=$( grep -o '[a-z]*_[a-z]*_[a-z]*' RAxML_bestTree.BS_ML_GTRCAT | sed 's/\_/\t/g' | awk '$2==$3{print $1,$2,$3}'  | wc -l )
printf "There are $num subspecies with species\n"
## There are 3 subspecies with species
```
The output names have been saved in a file called `subsp_check.txt` to further explore it:

```sh
# Run from `00_data_curation/rodentia_ctenohystrica/filter_aln/checked_aln` directory 
input="subsp_check.txt"
while IFS= read -r line 
do

subsp=$( echo $line | sed 's/ /\_/g' )
num_subsp=$( grep -ocP ''${subsp}'' ../../../genes.txt )
sp=$( echo $line | awk '{print $1,$2}' | sed 's/ /\_/g' )
num_sp=$( grep -ocP '\b'${sp}'\b' ../../../genes.txt )
printf "Num. sequences $sp : $num_sp\n"
printf "Num. sequences $subsp : $num_subsp\n"
printf "\n"

printf "Num. sequences $sp : $num_sp\n" >> sum_subspVSsp.txt
printf "Num. sequences $subsp : $num_subsp\n" >> sum_subspVSsp.txt
printf "\n" >> sum_subspVSsp.txt

done < $input
```

Output:

> **REMOVE species, priority to CYTB**   
   Num. sequences **cryptomys_hottentotus** : 3   
   Num. sequences cryptomys_hottentotus_hottentotus : 3   
   ```
   cryptomys_hottentotus_hottentotus --> "CYB"  "RNR1"
   cryptomys_hottentotus             --> "RNR1"            "ENSG00000110799" "ENSG00000112964"
   # We give priority to CYB, so we keep subsp - NOTE: They also cluster together in the tree, so no issue.
   ```
   
> **REMOVE subspecies, less genes**   
   Num. sequences trinomys_setosus : 4   
   Num. sequences **trinomys_setosus_setosus** : 3   
   ```
   trinomys_setosus         --> "CYB" "ENSG00000110799" "ENSG00000112964"
   trinomys_setosus_setosus --> "CYB"  "RNR1"
   ```

> **REMOVE subspecies, less genes***   
   Num. sequences thrichomys_apereoides : 6   
   Num. sequences **thrichomys_apereoides_apereoides** : 2
   
   
In addition, we carried out other checks with regards to the species names and their 
placement in the phylogeny. No taxa were indicated there to be removed, only some 
taxa are to be renamed. You can find a list of these checks in the file 
`rodentia_ctenohystrica_taxonomy_check.csv` as well as below:   

```
Genus    Abrocoma        Maybe      cuscomys_ashaninka   --> I have not found a study that shows the position of the species of Abrocoma
Family   Echimyidae      Yes        carterodon_sulcidens --> Species seems to be in the wrong position
Genus    Echimyidae      Yes        trinomys_setosus     --> Species seems to be in the wrong position
Genus    Pattonomys      Yes        N/A                  --> Genus should be monophylotic
Genus    Fukomys         Yes        fukomys_damarensis   --> Species seems to be in the wrong position
Genus    Atherurus       Yes        N/A                  --> I did not find much information. I think it should be monophyletic
Species  Cavia tschudii  Yes        N/A                  --> Species seems to be in the wrong position
Species  Cavia aperea    Yes        N/A                  --> Species seems to be in the wrong position
```

Then, we opened the `RAxML_bestTree.BS_ML_GTRCAT` with `FigTree` to individually check the taxa with dubious 
plaecment. Those that are to be deleted are coloured in red, while those that they cluster
within the same clade with the species they belong to but the placement might be dubious have
been coloured in orange. The names for the latter will be changed so they have a `~` at the end. 
This means that, even though this is the placement found in the ML tree, other research have
found they cluster elsewhere. Therefore, we flagged them for future and further analyses that might be carried out,
but that are not part of the scope of this project. 
The output files (`png` and `FigTree` files) have been saved in the
[`checked_aln/taxonomical_check`](checked_aln/taxonomical_check)
directory. 

The taxa to be renamed and/or removed are the following:   

```
Taxa to rename:
hystrix_brachyurus,hystrix_brachyura
trinomys_setosus_denigratus,trinomys_setosus_setosus
cryptomys_anomalus,cryptomys_hottentotus_anomalus
cryptomys_holosericeus,cryptomys_hottentotus_holosericeus 
cuscomys_ashaninka,cuscomys_ashaninka~
carterodon_sulcidens,carterodon_sulcidens~
pattonomys_occasius,pattonomys_occasius~
pattonomys_semivillosus,pattonomys_semivillosus~
atherurus_africanus,atherurus_africanus~
atherurus_macrourus,atherurus_macrourus~
cavia_tschudii,cavia_tschudii~
cavia_aperea,cavia_aperea~


Taxa to remove:
trinomys_setosus
fukomys_damarensis
```

Note that one of these species has been renamed to `trinomys_setosus_setosus`.
This subspecies is already in the tree, but it has been flagged 
to be removed. Therefore, now we check which genes are in `trinomys_setosus_denigratus`
because we had previously generated the R object `levels.checked.RData`:

```R
# Open Rscript `parse_lineage.R` and run from lines 1-15.
# Then, uncomment line 67 and run it to load the R object 

load( "levels.checked.RData" )
levels.checked$genus$Trinomys$trinomys_setosus_denigratus
# [1] "CYB"  "RNR1"
```

According to this, we decide that it is best to remove this species too. Therefore, the
updated list of taxa to be removed is the following:

```
Taxa to rename:
hystrix_brachyurus,hystrix_brachyura
cryptomys_anomalus,cryptomys_hottentotus_anomalus
cryptomys_holosericeus,cryptomys_hottentotus_holosericeus 
cuscomys_ashaninka,cuscomys_ashaninka~
carterodon_sulcidens,carterodon_sulcidens~
pattonomys_occasius,pattonomys_occasius~
pattonomys_semivillosus,pattonomys_semivillosus~
atherurus_africanus,atherurus_africanus~
atherurus_macrourus,atherurus_macrourus~
cavia_tschudii,cavia_tschudii~
cavia_aperea,cavia_aperea~

Taxa to remove:
trinomys_setosus
fukomys_damarensis
trinomys_setosus_denigratus
```
   
   
**NOTE:** Data subsets (both alignments, 12CP and 3CP) for Afrotheria, Xenarthra, Euarchonta, and Marsupialia had already undergone 
these checks before 2018 (i.e., they had already been "cleaned", while the data subset for Lagomorpha 
had not yet).
Therefore, you do not see the last part of the filtering described in this section 
in the corresponding `README.md` files for these data subsets or the directory
[`taxonomical_check`](checked_aln/taxonomical_check)
You will see that the `filter_aln` directory for these four data subsets
contains a csv file with the details that were followed to filter the corresponding alignments (12CP and 3CP). 
Note that this csv file is equivalent to the excel sheet you find in this directory for 
data subset Rodentia ctenohystrica.

## 2. Check alignment with 3CP partition
Before we concatenate the alignment with the first and second codon positions (12CP), file `alignment.phylip`
(you will find it once you unzip the file which link is provided above),
and the alignment with the third codon
positions (3CP), file `alignment_nt3cp.phylip` (you will find it once 
you unzip the file which link provided above, inside a directory called `checked_aln` directory),
we ran the following code to make sure they had both undergone the same filtering steps and that were 
at the same "filtering stage":

```sh
# Run the next code from `rodentia_ctenohystrica/filter_aln/checked_aln/unfiltered_aln`

# Check they have the same info
grep -o '[a-z].* ' alignment_nt3cp.phylip > names_3nt.txt
grep -o '[a-z].* ' ../../alignment.phylip > ../../names_aln.txt
diff names_3nt.txt ../../names_aln.txt # OK! 
```

Now that we are sure that both files with the 12CP-alignment (`alignment.phylip`) and 3CP-alignment
(`alignment_nt3cp.phylip`) are at the same "fitering stage"
we can proceed to concatenate the alignment with the third 
codon positions to the alignment with 12CPs using the R script
[`Concatenate_seqs_for_MCMCtree_unfiltered.R`](../../../01_alignments/Concatenate_seqs_for_MCMCtree_unfiltered.R)
so we can have the concatenated alignment. 

Instructions to follow:   

   * Open the RScript [`Concatenate_seqs_for_MCMCtree_unfiltered.R`](../../../01_alignments/Concatenate_seqs_for_MCMCtree_unfiltered.R) 
   in RStudio and change line 25 so it is `subt    <- "rodentia_ctenohystrica"` and uncomment line 27.
   This script will generate a concatenated alignment file with all partitions, as well as 
   one alignment file for each individual partitions, inside a new dir called `00_mammal_alns/rodentia_ctenohystrica/unfiltered/`
   inside [`00_Data_filtering/01_alignments/`](../../../01_alignments). 
   Log and RData files will be saved inside
   [`00_Data_filtering/01_alignments/Rout`](../../../01_alignments/Rout/log_concatenation).   
      > NOTE 1: Updated `partitions.txt` file inside `00_mammal_alns/rodentia_ctenohystrica/unfiltered/` generated.   
	  > NOTE 2: Output file called `check_appends.txt` to check that 3nt partition has been appended
	    to the right line of the alignment can be found in `Rout/Rdata". 
		If "0" in second column, sth wrong! So far, everything seems OK!   
   * Paths have been automatically set according to current file architecture in `mammals` dir, 
   do not change paths in the Rscript! 
   
**NOTE: You can see that we have used a different R script in this step because, as mentioned above,**
**the data for laurasiatheria cetartiodactyla were not at the same stage as data subsets for Afrotheria, Xenarthra,**
**Euarchonta and Marsupialia.**

Now, we are going to use this unfiltered alignment phylip to apply the next filtering steps. 

## 3. Apply the filtering step
Run the following commands from `rodentia_ctenohystrica/filter_aln/checked_aln`

```sh
# Run the next code from `rodentia_ctenohystrica/checked_aln/checked_aln`

# 1. Fix species names in `names.txt`. We used the RAxML file, the `names.txt` file, and 
#    the `lineages.txt` that had been generated when filtering the alignment 12CP.
#    The ML tree file is now saved inside `filter_aln/checked_aln`, from where the 
#    python script will be run. The `names.txt` and the `lineage.txt` files are in
#    the `filter_aln` dir, one up. The alignment is the concatenated one that was generated in the 
#    previous step and can be found inside
#   `01_alignments/00_mammal_alns/rodentia_ctenohystrica/unfiltered/rodentia_ctenohystrica.aln`.
#
#    Note that we just copy this file inside `checked_aln` to proceed with the filtering. We also generate 
#    a copy of the original tree file 
cp ../../../../01_alignments/00_mammal_alns/rodentia_ctenohystrica/unfiltered/rodentia_ctenohystrica.aln .
mkdir original_tree 
cp RAxML_bestTree.BS_ML_GTRCAT original_tree/RAxML_bestTree.BS_ML_GTRCAT
cat > _species_name_change.txt << EOF
hystrix_brachyurus,hystrix_brachyura
cryptomys_anomalus,cryptomys_hottentotus_anomalus
cryptomys_holosericeus,cryptomys_hottentotus_holosericeus
cuscomys_ashaninka,cuscomys_ashaninka~
carterodon_sulcidens,carterodon_sulcidens~
pattonomys_occasius,pattonomys_occasius~
pattonomys_semivillosus,pattonomys_semivillosus~
atherurus_africanus,atherurus_africanus~
atherurus_macrourus,atherurus_macrourus~
cavia_tschudii,cavia_tschudii~
cavia_aperea,cavia_aperea~
EOF

mkdir screen_logs
python ../../../../../../src/fix_species_names.py _species_name_change.txt RAxML_bestTree.BS_ML_GTRCAT rodentia_ctenohystrica.aln ../names.txt ../lineage.txt > screen_logs/01_log1_fixspecies.txt

##>> NOTE: I needed to update the python script as it was not working. I just modifed `fix_alignment()` 
##>>       so the output file was different from the one being read

## Move unncessary output files to `00_filt1/` and rename output alignment
mkdir 00_filt1
mv _species_name_change.txt lineage.txt.1.bak lineage.txt 00_filt1/
rm rodentia_ctenohystrica.aln
mv rodentia_ctenohystrica_out.phylip alignment.phylip

# 2. Remove species 

# 2.1. Create input file for python script `prune_alignment.py`. This text file 
#      contains a list of the species that need to be removed 
cat > _species_remove.txt << EOF
cryptomys_hottentotus
trinomys_setosus_setosus
thrichomys_apereoides
trinomys_setosus_denigratus
trinomys_setosus
fukomys_damarensis
EOF

# 2.2. Run set of four python scripts to generate pruned alignment without gaps:
python ../../../../../../src/phy2fasta.py alignment.phylip > alignment.fasta
python ../../../../../../src/prune_alignment.py --remove _species_remove.txt alignment.fasta > alignment_pruned.fasta
python ../../../../../../src/remove_gap_columns.py alignment_pruned.fasta > alignment_pruned_nogaps.fasta
python ../../../../../../src/fasta2phy.py alignment_pruned_nogaps.fasta > alignment_pruned_nogaps.phylip

# 2.3. Now prune the tree.
#      All the names of taxa to be removed are saved in a variable and pass it as the second argument instead 
#      of using a text file to `nw_pruned`
sp2rm=$( cat _species_remove.txt )
../../../../../../src/newick-utils-1.6/src/nw_prune RAxML_bestTree.BS_ML_GTRCAT $sp2rm > tmp.tree #COOL!

##  Move unnecessary output files to `01_filt2/`. Rename tree and alignment files.
mkdir 01_filt2 
mv *fasta alignment.phylip _species_remove.txt 01_filt2/
mv RAxML_bestTree.BS_ML_GTRCAT 01_filt2/RAxML_bestTree.BS_ML_GTRCAT_unpruned
mv tmp.tree RAxML_bestTree.BS_ML_GTRCAT
mv alignment_pruned_nogaps.phylip alignment.phylip

# 3. Remove duplicates 

# 3.1. Find duplicates
grep -oP '\w+_duplicate' RAxML_bestTree.BS_ML_GTRCAT # No duplicates!

# 4. I generate a list of names to keep track of current species after filtering!
grep -o '[a-z].* ' alignment.phylip > names_filt.txt
```

## 4. Data partitioning 
Now that we have the concatenated and filtered alignment ready, we need to generate the filtered
partitioned alignments by running the R script 
[`Partition_seqs_for_MCMCtree_after_filtering.R`](../../../01_alignments/Partition_seqs_for_MCMCtree_after_filtering.R).

Instructions to follow: 

   * Open the RScript [`Partition_seqs_for_MCMCtree_after_filtering.R`](../../../01_alignments/Partition_seqs_for_MCMCtree_after_filtering.R) 
   in RStudio and change line 24 so it is `subt    <- "rodentia_ctenohystrica"` and uncomment line 25.
   Now, we can run it from RStudio. This script will generate a concatenated alignment file with all partitions, as well as 
   one alignment file for each individual partitions,
   inside [`01_alignments`](../../../01_alignments).
   Log files and Rdata can be found [here](../../../01_alignments/Rout/log_concatenation). 
      > NOTE: Updated `partitions.txt` file inside `00_mammal_alns/rodentia_ctenohystrica/` generated.   
   * Paths have been automatically set according to current file architecture in `mammals` dir, 
   do not change paths in the Rscript! 

The final alignments generated at the end of this step can be found in the directory 
`before_adding_new_taxa` when unzipping
[this file](https://www.dropbox.com/s/9cjmclt4q2h9p3n/SeqBayesS2_Raln_ctenohystrica_ALLFILTERS.zip?dl=0).
Please read the next section to understand how the final alignments (also present in the zip 
file which link is provided above in the main directory) were generated!

# EXTRA FILTERING -- ADDING TAXA TO THE ALIGNMENT
When we first carried out the analysis with the tree topology as in step 1, we realised that 
we should try and include a correct sequence for *F. damarensis* as this species was present in 
the 72-sp tree (firs step of sequential Bayesian dating approach, data 1). In addition, we wanted 
to include the calibration for the node "Rodentia" to avoid future issues when grafting the 
data subtrees onto the backbone tree at the end of the sequential Bayesian dating approach. 

The procedure followed was the following:

## 1. Generate data subsets
The directory
[`00_perl_parsing`](extra_filtering/00_perl_parsing)
has all the details about the scripts used to reformat the alignments for this data subset before we 
added included the data for the new taxa. You can download all the reformatted alignments generated as well 
as the data needed to run MAFFT for the re-alignment
[here](https://www.dropbox.com/s/fkwspa58qrvplyr/SeqBayesS2_filteraln2_ctenohystrica_00_perl_parsing.zip?dl=0).

## 2. Generate new alignments with new taxon 
The directory
[`01_mafft`](extra_filtering/01_mafft) 
contains all the data used to generate the alignment that includes the *Ictidomys_tridecemlineatus* sequence.
This is the species which sequence data was needed to include the calibration for the node "Rodentia".
Please access this directory using the link provided above to go through the steps followed. You can download the 
`aln` directory [here](https://www.dropbox.com/s/seqc0qbwuxk8ev2/SeqBayesS2_filteraln2_ctenohystrica_01_mafft.zip?dl=0),
which should be saved inside the 
[`01_mafft`](extra_filtering/01_mafft) 
directory mentioned above to reproduce our results. 

## 3. Tests with *Fukomys damarensis*
We carried out several tests to decide which *Fukomys damarensis* sequence was going to be included 
in this data subset for rodentia ctenohystrica. Please access 
[`02_tests_fukomys`](extra_filtering/02_tests_fukomys)
using the link provided to go through the steps followed. You can download the tests carried out 
[here](https://www.dropbox.com/s/dcgmy6syj4z87l6/SeqBayesS2_filteraln2_ctenohystrica_02_tests_fukomys.zip?dl=0). 
You should unzip and save the content inside the
[`02_tests_fukomys`](extra_filtering/02_tests_fukomys)
directory to reproduce the results. 

## 4. Generate new alignments with *Fukomys damarensis* 
The directory
[`03_mafft`](extra_filtering/03_mafft) 
contains all the data used to generate the alignment with both *Ictidomys_tridecemlineatus* and
*Fukomys damarensis* sequences. The information about how to proceed with the mitochondrial alignments can be found in the directory 
[`04_mafft_mitaln`](extra_filtering/04_mafft_mitaln). 
Please access the two directories using the links provided above to go through the steps followed to generate 
the nuclear and mitochondrial alignments. You can download the `aln` directory 
[here](https://www.dropbox.com/s/h73spqokux4v230/SeqBayesS2_filteraln2_ctenohystrica_03_mafft.zip?dl=0)
for the former and
[here](https://www.dropbox.com/s/cw3vfm0uqdctbyq/SeqBayesS2_filteraln2_ctenohystrica_04_mafft.zip?dl=0)
for the latter.
You should unzip and save the content inside 
[here](extra_filtering/03_mafft) 
and 
[here](extra_filtering/04_mafft_mitaln),
respectively, so you can reproduce our results. 

## 5. Final alignments
[Here](https://www.dropbox.com/s/9cjmclt4q2h9p3n/SeqBayesS2_Raln_ctenohystrica_ALLFILTERS.zip?dl=0)
you can download the final alignments for each individual partition as well as the alignment with the five partitions 
in separate blocks. The previous alignments have been moved to a directory called 
`before_adding_new_taxa`, also present in this zip file.
