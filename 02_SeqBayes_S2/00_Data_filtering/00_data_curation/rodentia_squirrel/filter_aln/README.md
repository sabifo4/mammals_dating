# RODENTIA SQUIRREL (SCIURIDAE AND RELATED) - Filtering alignment
The initial files that you can find in this directory are the following:

```
rodentia_squirrel
   |- filter_aln
         |- checked_aln                        
         |     |- taxonomical_check
         |           |- sciuridae_and_related.png
         |           |- sciuridae_and_related_FigTree
         |           |- RAxML_bestTree.BS_ML_GTRCAT
         |     
         |- lineage.txt 
         |- names.txt        # Updated `names.txt` during filtering steps described below
         |- names_orig.txt   # This is the original `names.txt`. It is modified
         |                   # in one of the filtering steps described below
         |- parse_lineage.R 
         |- partitions.txt 
         |- README.md
         |- rodentia_squirrel_taxonomy_check.xlsx
         |- summary.html 
```
		 
The output files during the filtering step that is described below as well as the input 
files needed have been zipped in a file as they are very large. You can download 
this file [here]().
To start the filtering step, you should have the following files arranged in the file 
architecture detailed above (you can obtain the files once you unzip the file
provided in the link above): 

```
rodentia_squirrel 
   |- filter_aln
         |- checked_aln                        
         |     |- unfiltered_aln
         |     |    |- alignment_nt3cp.phylip        # Alignment with 3CP 
         |     |- RAxML_bestTree.BS_ML_GTRCAT        # Unfiltered tree (found in `original_tree`)
         |     |- taxonomical_check                  # Visual checks to evaluate 
         |          |- sciuridae_and_related.png     # dubious taxa placement
         |          |- sciuridae_and_related_FigTree
         |          |- RAxML_bestTree.BS_ML_GTRCAT    
         |          
         |- alignment.phylip # Alignment with 12CP
         |- lineage.txt 
         |- names.txt 
         |- names_orig.txt 
         |- parse_lineage.R 
         |- partitions.txt
         |- README.md
         |- rodentia_squirrel_taxonomy_check.xlsx
         |- summary.html 
```

Note that you will see more files within the zipped file provided above. These are output 
files that you will generate if you follow all the steps below. Feel free to keep them so you can 
make sure that you reproduce the same results that we show here. 

## 1. Taxonomic filtering 
Within the [`filter_aln`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_squirrel/filter_aln)
directory, you will find different files as detailed above. 
Specifically, the R script [`parse_lineage.R`](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_squirrel/filter_aln/parse_lineage.R)
was written to carry out a first taxonomic filtering. Note that this 
R script will run if you have the same file architecture in this GitHub repository (i.e., it uses 
a function within the R script [`Filter_lineages.R`](https://github.com/sabifo4/mammals_dating/blob/main/src/Filter_lineages.R)
in the `src` directory and the [`genes.txt`](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/genes.txt)
file. There are several messages printed out by this script as all the rodents are 
parsed at the same time. You might have seen that, due to the large number of rodent species, 
we have had to divide them into different data subsets: "rodentia_squirrel" (this one),
"rodentia_ctenohystrica", "rodentia_subtree1", and "rodentia_subtree2". Below, you will find 
the messages printed out that refer to this data subset, "rodentia_squirrel", i.e., 
sciuridae and related species:   

```
GENUS   Funisciurus  has  2  taxa. Species  funisciurus_pyrropus  has  2  genes that are not 
shared by any of the other taxa. It should be deleted!
GENUS   Funisciurus  has  2  taxa. Species  funisciurus_carruthersi  has  1  genes that are 
not shared by any of the other taxa. It should be deleted!
```
>> ACTION: funisciurus_pyrropus & funisciurus_carruthersi merged, the label is "funisciurus".

```
GENUS   Callospermophilus  has  4  taxa. Species  callospermophilus_lateralis  has  3  genes
that are not shared by any of the other taxa. It should be deleted!
```
>> ACTION: callospermophilus_lateralis should have been removed, but it was not. This is because 
>> subspecies callospermophilus_lateralis_trepidus was removed instead. After furhter exploring the 
>> R object, this was decided because there are less genes at the subspecies level, hence 
>> the sequence at the species level is kept. 
```
callospermophilus_lateralis          --> "ND2" "RNR1" "ENSG00000265203"
callospermophilus_lateralis_trepidus --> "CYB"
callospermophilus_madrensi           --> "CYB"
callospermophilus_saturatus          --> "CYB"
```   

```
GENUS   Microsciurus  has  3  taxa. Species  microsciurus_flaviventer  has  1  genes that are not
shared by any of the other taxa. It should be deleted!
```
>> ACTION: microsciurus_flaviventer should have been removed, but it was not. This is because 
>> subspecies microsciurus_flaviventer_sabanillae was removed instead. After exploring the 
>> R object, the same genes in microsciurus_flaviventer_sabanillae are shared at the genus 
>> level in microsciurus_alfari. Therefore, it was decided to keep microsciurus_flaviventer 
>>instead. 

```
microsciurus_alfari                 --> "CYB"             "RNR1"            "ENSG00000265203"
microsciurus_flaviventer            --> "ENSG00000166349"
microsciurus_flaviventer_sabanillae --> "RNR1"            "ENSG00000265203"
```

**NOTE**
Log files will be found in this same directory with the name `log_taxaNOTIN_<level>.txt`,
being `level` the one that has been checked (i.e., family, genus, order, or subfamily). In addition, 
you can also generate the `levels.checked.RData`, which you can use to explore the taxonomical 
levels. Remember that you will see more information about the other rodent taxa that are 
present in the other three data subsets.

# 2. First checks before applying filtering
First, we checked if there were any further species that should be removed:

```{sh}
# Run from `filter_aln/checked_aln` directory 
grep -o '[a-z]*_[a-z]*_[a-z]*' RAxML_bestTree.BS_ML_GTRCAT | sed 's/\_/\t/g' | awk '$2==$3{print $1,$2,$3}' | awk '$1!="aeb"{print $1,$2,$3}' > subsp_check.txt

# Count 
num=$( grep -o '[a-z]*_[a-z]*_[a-z]*' RAxML_bestTree.BS_ML_GTRCAT | sed 's/\_/\t/g' | awk '$2==$3{print $1,$2,$3}'  | awk '$1!="aeb"{print $1,$2,$3}' | wc -l )
printf "There are $num subspecies with species\n"
## There are 9 subspecies with species
```
The output names have been saved in a file called `subsp_check.txt` to further explore it:

```{sh}
# Run from `00_data_curation/rodentia_squirrel/filter_aln/checked_aln` directory 
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

> **REMOVE subspecies, less genes**   
   Num. sequences petaurista_leucogenys : 2   
   Num. sequences **petaurista_leucogenys_leucogenys** : 2   
   ```
   petaurista_leucogenys_leucogenys --> "CYB"
   petaurista_leucogenys            --> "CYB"
   ```
   
> **REMOVE subspecies, less genes**   
   Num. sequences aplodontia_rufa : 15   
   Num. sequences **aplodontia_rufa_rufa** : 3

> No action   
   Num. sequences urocitellus_brunneus : 0   
   Num. sequences urocitellus_brunneus_brunneus : 2

> No action   
   Num. sequences urocitellus_townsendi : 0   
   Num. sequences urocitellus_townsendi_townsendi : 2

> **REMOVE subspecies, less genes**   
   Num. sequences urocitellus_columbianus : 2   
   Num. sequences **urocitellus_columbianus_columbianus** : 2   
   ```
   urocitellus_columbianus_columbianus --> "CYB"
   urocitellus_columbianus             --> "CYB"
   ```
> No action   
   Num. sequences spermophilus_elegans : 0   
   Num. sequences spermophilus_elegans_elegans : 2

> No action   
   Num. sequences ictidomys_mexicanus : 0   
   Num. sequences ictidomys_mexicanus_mexicanus : 2

> **REMOVE subspecies, less genes**   
   Num. sequences otospermophilus_beecheyi : 4   
   Num. sequences **otospermophilus_beecheyi_beecheyi** : 2

> **REMOVE subspecies, less genes**   
   Num. sequences tamias_panamintinus : 4   
   Num. sequences **tamias_panamintinus_panamintinus** : 2
   

All together, the taxa to remove are the following:

```
Taxa to remove:
dryomys_nitedula_ssp_aeb_2014
```

In addition, we carried out other checks with regards to the species names and their 
placement in the phylogeny. You can find a list of these checks in the file 
`rodentia_squirrel_taxonomy_check.xlsx` as well as below:

```
Genus       Spermophilus               spermophilus_dauricus                                       --> It seems to be in the wrong place. It should cluster with the rest of Spermophilus species. References show  support for: a:"((S. dauricus, S.xanthoprymmus), S. suslicus);" b: "(S. dauricus, ((S. citellus,S.xanthoprymnus), S. suslicus));" or c:"(S. dauricus, ((S. citellus,S.taurensis), S. xanthoprymmus));". So maybe at S. dauricus to this clade if decided to move taxon. However, reference d found one mitochondrila genome of S. dauricus clustering S. vulgaris, as we observe. So decide.
Genus       Paraxerus                  funisciurus                                                 --> In the referenced papers, all paraxerus cluster in same clade, here funiscirus clusters with some of the paraxerus in a praphyletic clade. We can keep it like this or try to place funiscirus as the outgroup of a monophyletic clade for paraxerus.
Genus       Hylopetes                  iomys_horsfieldi, petinomys_setosus, petaurillus_kinlochii  --> Ref. a shows how H. lepidus, H. spadiceus, and H. nigripes form a first cluster and then H. alboniger and H. phayrei form a second cluster before the monphyletic clade is formed. We see the same but H. alboniger and H. phayrei are clustered in a group with Petaurillus kinlochii, Petinomys setosus, and Iomys horsfieldi. Ref. b shows that P. setosus clusters with H. alboniger. Ref. c shows that they cluster separately, though, as well as P. kinlochi and I. horsfieldi are part of a sister clade, not within the same clade. Decide if we move these three taxa or not.
Genus       Sciurus                    microsciurus_flaviventer,microsciurus_alfari                --> References attached have found they cluster in different clades with other sciurus species. We can keep it as it is.  
Genus       Sciurus                    syntheosciurus_brochus                                      --> Reference shows that S. brochus is part of a clade with S. variegatoides, S. granatensis, S. aureogaster. Not the same topology as we find, but as there is not that much information, we can leave it as it is. 
Genus       Sciurus                    rheithrosciurus_macrotis                                    --> Ref. a found a supertree based on gene trees where R. macrotis is also clustering with S. list and S. vulgaris as in our topology (although ours has also Spermophilus dauricus). Ref. b seems to support a clade with Sciurus, Microsciurus, and Rheithrosciurus species, being the latter the outgroup for a clade with Sciurus and Microsciurus taxa. R
Species     Otospermophilus beecheyi   otospermophilus_atricapillus                                --> Ref shows it is sister taxa with O. beecheyi, although it seems O. beecheyi should be a monophyletic group. Note: its homotypic synonym is "Otospermophilus beecheyi atricapillus".
Species     Urocitellus mollis         urocitellus_townsendi_vigilis                               --> Not much information about this, only reference where different species of U. townsendii do not tend ot cluster together. If needed, try to cluster U. mollis with its  subsp.
Species     Sundasciurus tenuis        sundasciurus_jentinki, sundasciurus_tahan                   --> According to references, it seems all species of S. tenuis cluster together. Maybe try (((S. tenuis, S.tenuis subsp), S.tahan), S.jentinki);
Species     Dremomys pernyi            dremomys_pernyi,dremomys_pernyi_owstoni                     --> Ref. suggests that subsp and sp. cluster together.
Species     Callosciurus erythraeus    callosciurus_finlaysonii                                    --> Ref. suggests that callosciurus_finlaysonii clusters with C. caniceps and C. inornatus, not forming a paraphyletic group with C. erythraeus. Maybe try to put subsp. And sp. C. erythraeus together and move C. finlaysonii.
Species     Tamias minimus             tamias_minimus, tamias_minimus_scrutator                    --> Ref.  A shows same clade with same species. Maybe it would be good to cluster subsp. with sp.  Ref b. shows that T. minimums clusters with T. alpinus, so try to place subsp. and sp. There.
Species     Sciurus granatensis        sciurus_colliaei                                            --> Not much information about species S. colliaei. Either keep like that or try to put subsp. With sp. Of S. granatensis together.
```

Then, we opened the `RAxML_bestTree.BS_ML_GTRCAT` with `FigTree` to individually check the taxa with dubious 
plaecment. Those that are to be deleted are coloured in red, while those that they cluster
within the same clade with the species they belong to but the placement might be dubious have
been coloured in orange. The names for the latter will be changed so they have a `~` at the end. 
This means that, even though this is the placement found in the ML tree, other research have
found they cluster elsewhere. Therefore, we flagged them for future and further analyses that might be carried out,
but that are not part of the scope of this project. 
The output files (`png` and `FigTree` files) have been saved in the
[`checked_aln/taxonomical_check`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_squirrel/filter_aln/checked_aln/taxonomical_check)
directory. 

The taxa to be renamed and/or removed are the following:   

```
Taxa to rename:
spermophilus_elegans_elegans,urocitellus_elegans_elegans
spermophilus_spilosoma_marginatus,xerospermophilus_spilosoma_marginatus
funisciurus,funisciurus~
iomys_horsfieldi,iomys_horsfieldi~
petinomys_setosus,petinomys_setosus~
petaurillus_kinlochii,petaurillus_kinlochii~
microsciurus_flaviventer,microsciurus_flaviventer~
microsciurus_alfari,microsciurus_alfari~
syntheosciurus_brochus,syntheosciurus_brochus~
rheithrosciurus_macrotis,rheithrosciurus_macrotis~
otospermophilus_atricapillus,otospermophilus_atricapillus~
urocitellus_townsendi_vigilis,urocitellus_townsendi_vigilis~
sundasciurus_jentinki,sundasciurus_jentinki~
dremomys_pernyi_owstoni,dremomys_pernyi_owstoni~
callosciurus_finlaysonii,callosciurus_finlaysonii~
tamias_minimus_scrutator,tamias_minimus_scrutator~
sciurus_colliaei,sciurus_colliaei~

Taxa to remove:
spermophilus_dauricus
```

**NOTE:**Data subsets (both alignments, 12CP and 3CP) for Afrotheria, Xenarthra, Euarchonta, and Marsupialia had already undergone 
these checks before 2018 (i.e., they had already been "cleaned", while the data subset for Lagomorpha 
had not yet).
Therefore, you do not see the last part of the filtering described in this section 
in the corresponding `README.md` files for these data subsets or the directory
[`taxonomical_check`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_squirrel/filter_aln/checked_aln/taxonomical_check)
You will see that the `filter_aln` directory for these four data subsets
contains a csv file with the details that were followed to filter the corresponding alignments (12CP and 3CP). 
Note that this csv file is equivalent to the excel sheet you find in this directory for 
data subset Laurasiatheria cetartiodactyla.


## 2. Check alignment with 3CP partition
Before we concatenate the alignment with the first and second codon positions (12CP), file `alignment.phylip`
(you will find it once you unzip the file which link is provided above),
and the alignment with the third codon
positions (3CP), file `alignment_nt3cp.phylip` (you will find it once 
you unzip the file which link provided above, inside a directory called `checked_aln` directory),
we ran the following code to make sure they had both undergone the same filtering steps and that were 
at the same "filtering stage":

```{sh} 
# Run the next code from `rodentia_squirrel/filter_aln/checked_aln/unfiltered_aln`

# Check they have the same info
grep -o '[a-z].* ' alignment_nt3cp.phylip > names_3nt.txt
grep -o '[a-z].* ' ../../alignment.phylip > ../../names_aln.txt
diff names_3nt.txt ../../names_aln.txt # OK! 
```

Now that we are sure that both files with the 12CP-alignment (`alignment.phylip`) and 3CP-alignment
(`alignment_nt3cp.phylip`) are at the same "fitering stage"
we can proceed to concatenate the alignment with the third 
codon positions to the alignment with 12CPs using the R script
[`Concatenate_seqs_for_MCMCtree_unfiltered.R`](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/01_alignments/Concatenate_seqs_for_MCMCtree_unfiltered.R)
so we can have the concatenated alignment. 

Instructions to follow:   

   * Open the RScript [`Concatenate_seqs_for_MCMCtree_unfiltered.R`](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/01_alignments/Concatenate_seqs_for_MCMCtree_unfiltered.R) 
   in RStudio and change line 25 so it is `subt    <- "rodentia_squirrel"`, uncomment line 27, 
   and comment line 29. Now, we can run it from RStudio.
   This script will generate a concatenated alignment file with all partitions, as well as 
   one alignment file for each individual partitions, inside a new dir called `00_mammal_alns/rodentia_squirrel/unfiltered/`
   inside [`00_Data_filtering/01_alignments/`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/01_alignments). 
   Log and RData files will be saved inside
   [`00_Data_filtering/01_alignments/Rout`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/01_alignments/Rout/log_concatenation).   
      > NOTE 1: Updated `partitions.txt` file inside `00_mammal_alns/rodentia_squirrel/unfiltered/` generated.   
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
Run the following commands from `rodentia_squirrel/filter_aln/checked_aln`

```{sh}
# Run the next code from `rodentia_squirrel/filter_aln/checked_aln`

# 1. Fix species names in `names.txt`. We used the RAxML file, the `names.txt` file, and 
#    the `lineages.txt` that had been generated when filtering the alignment 12CP.
#    The ML tree file is now saved inside `filter_aln/checked_aln`, from where the 
#    python script will be run. The `names.txt` and the `lineage.txt` files are in
#    the `filter_aln` dir, one up. The alignment is the concatenated one that was generated in the 
#    previous step and can be found inside
#   `01_alignments/00_mammal_alns/rodentia_squirrel/unfiltered/rodentia_squirrel.aln`.
#    NOTE: I had to remove the names of the two taxa that had been combined to lead to 
#          "funisciurus" and "microsciurus_flaviventer" in the `names.txt` file, otherwise 
#          filtering was not working. I kept original file as `names_orig.txt`.
#
#    Note that we just copy this file inside `checked_aln` to proceed with the filtering. We also generate 
#    a copy of the original tree file 
cp ../../../../01_alignments/00_mammal_alns/rodentia_squirrel/unfiltered/rodentia_squirrel.aln .
mkdir original_tree 
cp RAxML_bestTree.BS_ML_GTRCAT original_tree/RAxML_bestTree.BS_ML_GTRCAT
cat > _species_name_change.txt << EOF
spermophilus_elegans_elegans,urocitellus_elegans_elegans
spermophilus_spilosoma_marginatus,xerospermophilus_spilosoma_marginatus
funisciurus,funisciurus~
iomys_horsfieldi,iomys_horsfieldi~
petinomys_setosus,petinomys_setosus~
petaurillus_kinlochii,petaurillus_kinlochii~
microsciurus_alfari,microsciurus_alfari~
syntheosciurus_brochus,syntheosciurus_brochus~
rheithrosciurus_macrotis,rheithrosciurus_macrotis~
otospermophilus_atricapillus,otospermophilus_atricapillus~
urocitellus_townsendi_vigilis,urocitellus_townsendi_vigilis~
sundasciurus_jentinki,sundasciurus_jentinki~
dremomys_pernyi_owstoni,dremomys_pernyi_owstoni~
callosciurus_finlaysonii,callosciurus_finlaysonii~
tamias_minimus_scrutator,tamias_minimus_scrutator~
sciurus_colliaei,sciurus_colliaei~
microsciurus_flaviventer,microsciurus_flaviventer~
EOF

mkdir screen_logs
python ../../../../../../src/fix_species_names.py _species_name_change.txt RAxML_bestTree.BS_ML_GTRCAT rodentia_squirrel.aln ../names.txt ../lineage.txt > screen_logs/01_log1_fixspecies.txt

##>> NOTE: I needed to update the python script as it was not working. I just modifed `fix_alignment()` 
##>>       so the output file was different from the one being read

## Move unncessary output files to `00_filt1/` and rename output alignment
mkdir 00_filt1
mv _species_name_change.txt 00_filt1/
rm rodentia_squirrel.aln
mv rodentia_squirrel_out.phylip alignment.phylip

# 2. Remove species 

# 2.1. Create input file for python script `prune_alignment.py`. This text file 
#      contains a list of the species that need to be removed 
cat > _species_remove.txt << EOF
dryomys_nitedula_ssp_aeb_2014
petaurista_leucogenys_leucogenys
aplodontia_rufa_rufa
urocitellus_columbianus_columbianus
otospermophilus_beecheyi_beecheyi
tamias_panamintinus_panamintinus
spermophilus_dauricus
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
[`Partition_seqs_for_MCMCtree_after_filtering.R`](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/01_alignments/Partition_seqs_for_MCMCtree_after_filtering.R).

Instructions to follow: 

   * Open the RScript [`Partition_seqs_for_MCMCtree_after_filtering.R`](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/01_alignments/Partition_seqs_for_MCMCtree_after_filtering.R) 
   in RStudio and change line 24 so it is `subt    <- "rodentia_squirrel"` and uncomment line 26.
   Now, we can run it from RStudio. This script will generate a concatenated alignment file with all partitions, as well as 
   one alignment file for each individual partitions,
   inside [`01_alignments`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/01_alignments).
   Log files and Rdata can be found [here](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/01_alignments/Rout/log_concatenation). 
      > NOTE: Updated `partitions.txt` file inside `00_mammal_alns/rodentia_squirrel/` generated.   
   * Paths have been automatically set according to current file architecture in `mammals` dir, 
   do not change paths in the Rscript! 

The final alignments generated at the end of this step can be downloaded from 
[here]().
They should be saved here if the same file architecture as the one set in the R scripts 
is to be used: `00_Data_filtering/01_alignments/00_mammals_alns/rodentia_squirrel`.