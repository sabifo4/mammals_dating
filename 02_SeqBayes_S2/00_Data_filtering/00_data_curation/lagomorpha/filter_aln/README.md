# LAGOMORPHA - Filtering alignment
The initial files that you can find in this directory are the following:

```
lagomorpha
   |- filter_aln
         |- checked_aln                        
         |     |- taxonomical_check
         |           |- lagomorpha.png
         |           |- lagomorpha_FigTree
         |           |- RAxML_bestTree.BS_ML_GTRCAT
         |     
         |- Lagomorpha_taxonomy_check.xlsx
         |- lineage.txt 
         |- names.txt 
         |- parse_lineage.R 
         |- partitions.txt 
         |- README.md
         |- summary.html 
```
		 
The output files during the filtering step that is described below as well as the input 
files needed have been zipped in a file as they are very large. You can download 
this file [here](https://www.dropbox.com/s/flc0tixotz286b4/SeqBayesS2_filtaln_lagmorpha.zip?dl=0).
To start the filtering step, you should have the following files arranged in the file 
architecture detailed above (you can obtain the files once you unzip the file
provided in the link above): 

```
lagomorpha 
   |- filter_aln
         |- checked_aln                        
         |     |- unfiltered_aln
         |     |    |- alignment_nt3cp.phylip      # Alignment with 3CP 
         |     |- RAxML_bestTree.BS_ML_GTRCAT      # Unfiltered tree (found in `original_tree`)
         |     |- taxonomical_check                # Visual checks to evaluate 
         |          |- lagomorpha.png              # dubious taxa placement
         |          |- lagomorpha_FigTree
         |          |- RAxML_bestTree.BS_ML_GTRCAT    
         |          
         |- alignment.phylip # Alignment with 12CP
         |- Lagomorpha_taxonomy_check.xlsx
         |- lineage.txt 
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
Within the [`filter_aln`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/lagomorpha/filter_aln)
directory, you will find different files as detailed above. 
Specifically, the R script [`parse_lineage.R`](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/lagomorpha/filter_aln/parse_lineage.R)
was written to carry out a first taxonomic filtering. Note that this 
R script will run if you have the same file architecture in this GitHub repository (i.e., it uses 
a function within the R script [`Filter_lineages.R`](https://github.com/sabifo4/mammals_dating/blob/main/src/Filter_lineages.R)
in the `src` directory and the [`genes.txt`](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/genes.txt)
file. The messages printed out by this script are the following:   

```
GENUS   Lepus  has  34  taxa. Species  lepus_starcki  has  1  genes that are not shared by any
of the other taxa. It should be deleted!
GENUS   Lepus  has  34  taxa. Species  lepus_crawshayi  has  3  genes that are not shared by
any of the other taxa. It should be deleted!
```
>> ACTION: both species were deleted. 

**NOTE**
Log files will be found in this same directory with the name `log_taxaNOTIN_<level>.txt`,
being `level` the one that has been checked (i.e., family, genus, order, or subfamily). In addition, 
you can also generate the `levels.checked.RData`, which you can use to explore the taxonomical 
levels. 

# 2. First checks before applying filtering
First, we checked if there were any further species that should be removed:

```sh
# Run from `filter_aln/checked_aln` directory 
grep -o '[a-z]*_[a-z]*_[a-z]*' RAxML_bestTree.BS_ML_GTRCAT | sed 's/\_/\t/g' | awk '$2==$3{print $1,$2,$3}' > subsp_check.txt

# Count 
num=$( grep -o '[a-z]*_[a-z]*_[a-z]*' RAxML_bestTree.BS_ML_GTRCAT | sed 's/\_/\t/g' | awk '$2==$3{print $1,$2,$3}' | wc -l)
printf "There are $num subspecies with species\n"
## There are 6 subspecies with species
```
The output names were saved in a file called `subsp_check.txt` to further explore it:

```sh
# Run from `00_data_curation/lagomorpha/filter_aln/checked_aln` directory 
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
   Num. sequences lepus_europaeus : 45   
   Num. sequences **lepus_europaeus_europaeus** : 1

> **REMOVE subspecies, less genes**   
   Num. sequences lepus_oiostolus : 8   
   Num. sequences **lepus_oiostolus_oiostolus**: 2

> **REMOVE subspecies, less genes**   
   Num. sequences oryctolagus_cuniculus : 177   
   Num. sequences **oryctolagus_cuniculus_cuniculus**: 11

> **REMOVE subspecies, less genes**   
   Num. sequences ochotona_dauurica : 5   
   Num. sequences **ochotona_dauurica_dauurica** : 4   
   ```
   # ochotona_dauurica          -->"CYB"             "ENSG00000121966" "ENSG00000174697" "ENSG00000213088"
   # ochotona_dauurica_dauurica -->"CYB" "ND4"
   ```

> **REMOVE species, less genes**   
   Num. sequences **ochotona_cansus** : 1   
   Num. sequences ochotona_cansus_cansus : 5   
   ```
   # ochotona_cansus        --> "ENSG00000174697"
   # ochotona_cansus_cansus -->"CYB"             "ND4"             "ENSG00000174697"
   ```
   
> **REMOVE subspecies, less genes**
   Num. sequences ochotona_thibetana : 4   
   Num. sequences **ochotona_thibetana_thibetana** : 2   
   ```
   # ochotona_thibetana           --> "CYB" "ND4"
   # ochotona_thibetana_thibetana -->  "CYB"
   ```

All together, the taxa to remove are the following:

```
Taxa to remove:
ochotona_thibetana_thibetana
lepus_oiostolus_oiostolus
ochotona_cansus
oryctolagus_cuniculus_cuniculus
lepus_europaeus_europaeus
ochotona_dauurica_dauurica
```

In addition, we carried out other checks with regards to the species names and their 
placement in the phylogeny. You can find a list of these checks in the file 
`Lagomorpha_taxonomy_check.xlsx` as well as below:

```
Lagomorpha	Genus	sylvilagus_bachmani            --> Seems it should cluster in same clade as S.palustris and S.aquaticus
Lagomorpha	Species	ochotona_thibetana_nangqenica  --> All subspecies of Thibetana seem to cluster together
Lagomorpha	Species	lepus_oiostolus_qinghaiensis   --> It seems it should cluster with species L. oiostolus
Lagomorpha	Species	lepus_insularis                --> Clade paraphyletic with L. insularis and L. californicus. Subsp. L. californicus xanti should cluster with clade with rest of L. californicus sp.
Lagomorpha	Species	ochotona_curzoniae_melanostoma --> It seems it should cluster with O. curzoniae, not forming a paraphyletic clade with O. thibetana.
Lagomorpha	Species	lepus_capensis_mediterraneus   --> Very far from Lepus capensis. Check if they should cluster together.
Lagomorpha	Species	ochotona_alpina_scorodumovi    --> ynonym: Ochotona scorodumovi. Reference does not use it as a subspecies, but as a species. Maybe leave as it is: "Lissovsky et al. (2007) also suggested the presence of at least one additional species, O. scorodumovi, which may be conspecific with O. mantchurica"
```

Then, we opened the `RAxML_bestTree.BS_ML_GTRCAT` with `FigTree` to individually check the taxa with dubious 
plaecment. Those that are to be deleted are coloured in red, while those that they cluster
within the same clade with the species they belong to but the placement might be dubious have
been coloured in orange. The names for the latter will be changed so they have a `~` at the end. 
This means that, even though this is the placement found in the ML tree, other research have
found they cluster elsewhere. Therefore, we flagged them for future and further analyses that might be carried out,
but that are not part of the scope of this project. 
The output files (`png` and `FigTree` files) have been saved in the
[`checked_aln/taxonomical_check`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/lagomorpha/filter_aln/checked_aln/taxonomical_check)
directory. 

The taxa to be renamed and/or removed are the following:

```
Taxa to rename:
ochotona_alpina_scorodumovi,ochotona_scorodumovi
ochotona_thibetana_nangqenica,ochotona_thibetana_nangqenica~
lepus_oiostolus_qinghaiensis,lepus_oiostolus_qinghaiensis~
lepus_insularis,lepus_insularis~
ochotona_curzoniae_melanostoma,ochotona_curzoniae_melanostoma~
lepus_capensis_mediterraneus,lepus_capensis_mediterraneus~

Taxa to remove:
sylvilagus_bachmani
```
   
   
**NOTE:** Data subsets (both alignments, 12CP and 3CP) for Afrotheria, Xenarthra, Euarchonta, and Marsupialia had already undergone 
these checks before 2018 (i.e., they had already been "cleaned", while the data subset for Lagomorpha 
had not yet).
Therefore, you do not see the last part of the filtering described in this section 
in the corresponding `README.md` files for these data subsets or the directory
[`taxonomical_check`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/lagomorpha/filter_aln/checked_aln/taxonomical_check)
You will see that the `filter_aln` directory for these four data subsets
contains a csv file with the details that were followed to filter the corresponding alignments (12CP and 3CP). 
Note that this csv file is equivalent to the excel sheet you find in this directory for 
data subset Lagomorpha.

## 2. Check alignment with 3CP partition
Before we concatenate the alignment with the first and second codon positions (12CP), file `alignment.phylip`
(you will find it once you unzip the file which link is provided above),
and the alignment with the third codon
positions (3CP), file `alignment_nt3cp.phylip` (you will find it once 
you unzip the file which link provided above, inside a directory called `checked_aln` directory),
we ran the following code to make sure they had both undergone the same filtering steps and that were 
at the same "filtering stage":

```sh
# Run the next code from `lagomorpha/filter_aln/checked_aln/unfiltered_aln`
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
   in RStudio and change line 25 so it is `subt    <- "lagomorpha"`, uncomment line 27, and comment line 29. Now, we
   can run it from RStudio. This script will generate a concatenated alignment file with all
   partitions, as well as one alignment file for each individual partitions, inside a new dir called `00_mammal_alns/lagomorpha/unfiltered/`
   inside [`00_Data_filtering/01_alignments/`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/01_alignments). 
   Log and RData files will be saved inside
   [`00_Data_filtering/01_alignments/Rout`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/01_alignments/Rout/log_concatenation).   
      > NOTE 1: Updated `partitions.txt` file inside `00_mammal_alns/lagomorpha/unfiltered/` generated.   
	  > NOTE 2: Output file called `check_appends.txt` to check that 3nt partition has been appended
	    to the right line of the alignment can be found in `Rout/Rdata". 
		If "0" in second column, sth wrong! So far, everything seems OK!   
   * Paths have been automatically set according to current file architecture in `mammals` dir, 
   do not change paths in the Rscript! 
   
**NOTE: You can see that we have used a different R script in this step because, as mentioned above,**
**the data for lagomorpha were not at the same stage as data subsets for Afrotheria, Xenarthra,**
**Euarchonta and Marsupialia.**

Now, we are going to use this unfiltered alignment phylip to apply the next filtering steps. 

## 3. Apply the filtering step
Run the following commands from `lagomorpha/filter_aln/checked_aln`

```sh
# Run the next code from `lagomorpha/filter_aln/checked_aln`

# 1. Fix species names in `names.txt`. We used the RAxML file, the `names.txt` file, and 
#    the `lineages.txt` that had been generated when filtering the alignment 12CP.
#    The ML tree file is now saved inside `filter_aln/checked_aln`, from where the 
#    python script will be run. The `names.txt` and the `lineage.txt` files are in
#    the `filter_aln` dir, one up. The alignment is the concatenated one that was generated in the 
#    previous step and can be found inside `01_alignments/00_mammal_alns/lagomorpha/unfiltered/lagomorpha.aln`.
#    Note that we just copy this file inside `checked_aln` to proceed with the filtering. We also generate 
#    a copy of the original tree file 
cp ../../../../01_alignments/00_mammal_alns/lagomorpha/unfiltered/lagomorpha.aln .
mkdir original_tree 
cp RAxML_bestTree.BS_ML_GTRCAT original_tree/RAxML_bestTree.BS_ML_GTRCAT
cat > _species_name_change.txt << EOF
ochotona_alpina_scorodumovi,ochotona_scorodumovi
ochotona_thibetana_nangqenica,ochotona_thibetana_nangqenica~
lepus_oiostolus_qinghaiensis,lepus_oiostolus_qinghaiensis~
lepus_insularis,lepus_insularis~
ochotona_curzoniae_melanostoma,ochotona_curzoniae_melanostoma~
lepus_capensis_mediterraneus,lepus_capensis_mediterraneus~
EOF

mkdir screen_logs
python ../../../../../../src/fix_species_names.py _species_name_change.txt RAxML_bestTree.BS_ML_GTRCAT lagomorpha.aln ../names.txt ../lineage.txt > screen_logs/01_log1_fixspecies.txt

##>> NOTE: I needed to update the python script as it was not working. I just modifed `fix_alignment()` 
##>>       so the output file was different from the one being read

## Move unncessary output files to `00_filt1/` and rename output alignment
mkdir 00_filt1
mv _species_name_change.txt lineage.txt lineage.txt.1.bak 00_filt1/
rm lagomorpha.aln
mv lagomorpha_out.phylip alignment.phylip

# 2. Remove species 

# 2.1. Create input file for python script `prune_alignment.py`. This text file 
#    contains a list of the species that need to be removed 
cat > _species_remove.txt << EOF
ochotona_thibetana_thibetana
lepus_oiostolus_oiostolus
ochotona_cansus
oryctolagus_cuniculus_cuniculus
lepus_europaeus_europaeus
ochotona_dauurica_dauurica
sylvilagus_bachmani
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
   in RStudio and change line 24 so it is `subt    <- "lagomorpha"` and uncomment line 26. Now, we can run it from RStudio. 
   This script will generate a concatenated alignment file with all partitions, as well as 
   one alignment file for each individual partitions,
   inside [`01_alignments`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/01_alignments).
   Log files and Rdata can be found [here](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/01_alignments/Rout/log_concatenation). 
      > NOTE: Updated `partitions.txt` file inside `00_mammal_alns/lagomorpha/`.   
   * Paths have been automatically set according to current file architecture in `mammals` dir, 
   do not change paths in the Rscript! 

The final alignments generated at the end of this step can be downloaded from 
[here](https://www.dropbox.com/s/ohvak6nwq9ldsv9/SeqBayesS2_Raln_lagomorpha.zip?dl=0).
They should be saved here if the same file architecture as the one set in the R scripts 
is to be used: `00_Data_filtering/01_alignments/00_mammals_alns/lagomorpha`.