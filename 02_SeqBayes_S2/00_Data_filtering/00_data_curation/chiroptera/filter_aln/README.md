# LAURASIATHERIA CHIROPTERA - Filtering alignment
The initial files that you can find in this directory are the following:

```
chiroptera
   |- filter_aln
         |- checked_aln                        
         |     |- taxonomical_check
         |           |- laurasiatheria_chiroptera.png
         |           |- laurasiatheria_chiroptera_FigTree
         |           |- RAxML_bestTree.BS_ML_GTRCAT
         |     
         |- extra_filtering   <-- You will read more about it at the end of this guideline                     
         |     
         |- laurasiatheria_chiroptera_taxonomy_check.csv
         |- lineage.txt 
         |- names.txt       
         |- parse_lineage.R 
         |- partitions.txt 
         |- README.md
         |- summary.html 
```
		 
The output files during the filtering step that is described below as well as the input 
files needed have been zipped in a file as they are very large. You can download 
this file [here](https://www.dropbox.com/s/keb98e01ekdhwih/SeqBayesS2_filtaln_chiroptera.zip?dl=0).
To start the filtering step, you should have the following files arranged in the file 
architecture detailed above (you can obtain the files once you unzip the file
provided in the link above): 

```
chiroptera 
   |- filter_aln
         |- checked_aln                        
         |     |- unfiltered_aln
         |     |    |- alignment_nt3cp.phylip   # Alignment with 3CP 
         |     |- RAxML_bestTree.BS_ML_GTRCAT   # Tree you can find in `original_tree` directory
         |     |- taxonomical_check                      # Visual checks to evaluate 
         |          |- laurasiatheria_chiroptera.png     # dubious taxa placement
         |          |- laurasiatheria_chiroptera_FigTree
         |          |- RAxML_bestTree.BS_ML_GTRCAT    
         |     
         |- extra_filtering   <-- You will read more about it at the end of this guideline                     
         |        
         |- alignment.phylip # Alignment with 12CP
         |- laurasiatheria_chiroptera_taxonomy_check.csv
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
Within the `filter_aln`
directory, you will find different files as detailed above. 
Specifically, the R script [`parse_lineage.R`](parse_lineage.R)
was written to carry out a first taxonomic filtering. Note that this 
R script will run if you have the same file architecture in this GitHub repository (i.e., it uses 
a function within the R script [`Filter_lineages.R`](../../../../../src/Filter_lineages.R)
in the `src` directory and the [`genes.txt`](../../genes.txt)
file. The messages printed out by this script are the following:   

```
GENUS   Saccolaimus  has  2  taxa. Species  saccolaimus_flaviventris  has  1  genes that are not shared by any
of the other taxa. It should be deleted!
GENUS   Saccolaimus  has  2  taxa. Species  saccolaimus_saccolaimus  has  1  genes that are not shared by any
of the other taxa. It should be deleted!
```
>> ACTION: saccolaimus_flaviventris and saccolaimus_saccolaimus were merged. Label is "saccolaimus".

```
GENUS   Choeroniscus  has  2  taxa. Species  choeroniscus_godmani  has  1  genes that are not shared by any of
the other taxa. It should be deleted!
GENUS   Choeroniscus  has  2  taxa. Species  choeroniscus_minor  has  6  genes that are not shared by any of
the other taxa. It should be deleted!
```
>> ACTION: choeroniscus_godmani and choeroniscus_minor were mereged. Label is "choeroniscus".

```
GENUS   Nyctophilus  has  3  taxa. Species  nyctophilus_arnhemensis  has  1  genes that are not shared by any
of the other taxa. It should be deleted!
```
>> ACTION: nyctophilus_arnhemensis was removed.

```
GENUS   Cheiromeles  has  2  taxa. Species  cheiromeles_torquatus  has  2  genes that are not shared by any of
the other taxa. It should be deleted!
GENUS   Cheiromeles  has  2  taxa. Species  cheiromeles_torquatus  has  1  genes that are not shared by any of
the other taxa. It should be deleted!
```
>> ACTION: cheiromeles_torquatus and cheiromeles_torquatus were merged. Label is "cheiromeles".

```
GENUS   Mormopterus  has  12  taxa. Species  mormopterus_jugularis  has  2  genes that are not shared by any
of the other taxa. It should be deleted!
```
>> ACTION: mormopterus_jugularis was removed. 

```
GENUS   Cynomops  has  2  taxa. Species  cynomops_planirostris  has  1  genes that are not shared by any
of the other taxa. It should be deleted!
GENUS   Cynomops  has  2  taxa. Species  cynomops_abrasus  has  4  genes that are not shared by any of
the other taxa. It should be deleted!
```
>> ACTION: cynomops_planirostris cynomops_abrasus were merged. Label is "saccolaimus".

**NOTE**
Log files will be found in this same directory with the name `log_taxaNOTIN_<level>.txt`,
being `level` the one that has been checked (i.e., family, genus, order, or subfamily). In addition, 
you can also generate the `levels.checked.RData`, which you can use to explore the taxonomical 
levels.

# 2. First checks before applying filtering
First, we checked if there were any further species that should be removed:

```sh
# Run from `checked_aln` directory 
grep -o '[a-z]*_[a-z]*_[a-z]*' RAxML_bestTree.BS_ML_GTRCAT | sed 's/\_/\t/g' | awk '$2==$3{print $1,$2,$3}' | awk '$1!="khz"{print $1,$2,$3}' > checked_aln.txt

# Count 
num=$( grep -o '[a-z]*_[a-z]*_[a-z]*' RAxML_bestTree.BS_ML_GTRCAT | sed 's/\_/\t/g' | awk '$2==$3{print $1,$2,$3}' | awk '$1!="pdh"{print $1,$2,$3}' | wc -l)
printf "There are $num subspecies with species\n"
## There are 19 subspecies with species
```
The output names have been saved in a file called `checked_aln.txt` to further explore it:

```sh
# Run from `checked_aln` directory 
input="checked_aln.txt"
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
   Num. sequences miniopterus_schreibersii : 26   
   Num. sequences **miniopterus_schreibersii_schreibersii** : 2

> NO ACTION   
   Num. sequences myotis_nesopolus : 0   
   Num. sequences myotis_nesopolus_nesopolus : 2

> **REMOVE subspecies, less genes**   
Num. sequences myotis_daubentonii : 16   
Num. sequences **myotis_daubentonii_daubentonii** : 2

> **REMOVE species, less genes**   
   Num. sequences **myotis_blythii** : 3   
   Num. sequences myotis_blythii_blythii : 5   
   ```
   # myotis_blythii_blythii --> "CYB" "ND1" "ENSG00000175097"   
   # myotis_blythii         --> "ND1" "ENSG00000129538"   
   ```

> **REMOVE subspecies, less genes**   
   Num. sequences myotis_adversus : 8   
   Num. sequences **myotis_adversus_adversus** : 4   
   ```
   # myotis_adversus_adversus --> "CYB" "ND1"   
   # myotis_adversus          --> "CYB"  "ND1"  "ND2"  "RNR1" "RNR2"   
   ```

> **REMOVE subspecies, less genes**   
   Num. sequences myotis_horsfieldii : 8   
   Num. sequences **myotis_horsfieldii_horsfieldii** : 3

> NO ACTION   
   Num. sequences eptesicus_isabellinus : 0   
   Num. sequences eptesicus_isabellinus_isabellinus : 1

> **REMOVE subspecies, less genes**   
   Num. sequences eptesicus_serotinus : 28   
   Num. sequences **eptesicus_serotinus_serotinus** : 1

> **REMOVE subspecies, less genes**   
   Num. sequences plecotus_austriacus : 6   
   Num. sequences **plecotus_austriacus_austriacus** : 2   
   ```
   # plecotus_austriacus_austriacus --> "ND1"   
   # plecotus_austriacus            --> "ND1" "RNR1" "RNR2" "ENSG00000152592" "ENSG00000175097"   
   ```

> **REMOVE subspecies, less genes**   
   Num. sequences lasiurus_cinereus : 7   
   Num. sequences **lasiurus_cinereus_cinereus** : 1

> **REMOVE subspecies, less genes**   
   Num. sequences artibeus_planirostris : 4   
   Num. sequences **artibeus_planirostris_planirostris** : 2   
   ```
   # artibeus_planirostris_planirostris --> "CYB"   
   # artibeus_planirostris              --> "CYB" "RNR2" "ENSG00000112964"   
   ```

> **REMOVE subspecies, less genes**   
   Num. sequences artibeus_jamaicensis : 42   
   Num. sequences **artibeus_jamaicensis_jamaicensis** : 2

> **REMOVE subspecies, less genes**   
   Num. sequences rhinolophus_darlingi : 4   
   Num. sequences **rhinolophus_darlingi_darlingi** : 2   
   ```
   # rhinolophus_darlingi_darlingi --> "CYB"   
   # rhinolophus_darlingi          --> "CO2" "CYB"   
   ```

> NO ACTION   
   Num. sequences hipposideros_alongensis : 0   
   Num. sequences hipposideros_alongensis_alongensis : 2

> NO ACTION   
   Num. sequences hipposideros_turpis : 0   
   Num. sequences hipposideros_turpis_turpis : 2

> NO ACTION   
   Num. sequences pteropus_seychellensis : 0   
   Num. sequences pteropus_seychellensis_seychellensis : 3

> NO ACTION   
   Num. sequences pteropus_pelewensis : 0   
   Num. sequences pteropus_pelewensis_pelewensis : 3


All together, the taxa to remove are the following:

```
Taxa to remove:
rhinolophus_rufus
```

In addition, we carried out other checks with regards to the species names and their 
placement in the phylogeny. You can find a list of these checks in the file 
`laurasiatheria_chiroptera_taxonomy_check.csv` as well as below:

```
Genus   Phyllostomus  Maybe  N/A                      --> Maybe it is worth forcing a monophyletic clade and compare likelihoods. Some authors find Phyllostomus as monophyletic (see Rojas et al 2016).
Genus   Chilonatalus  Maybe  nyctiellus_lepidus       --> This genus seems to be monophyletic; it's worth checking comparing likelihood after forcing a monophyletic clade.
Genus   Plecotus      Yes    plecotus_macrobullaris   --> Plecotus macrobullaris should group with the other species of Plecotus.
Genus   Dobsonia      Yes    aproteles_bulmerae       --> Dobsonia appears to be a monophyletic genus.
Genus   Myotis        Yes    myotis_moluccarum        --> This species myotis_moluccarum should group with the other Myotis species. It is in a very weird position.
Genus   Pteropus      Yes    acerodon                 --> Maybe Pteropus should form a monophyletic group. The species P. personatus is in a weird position, but there are only 2 sequences representing this taxa. Colgan and Costa 2012 recovered a monophyletic group, but they did not have any species of acerodon.
Genus   Vespertilio   Yes    eptesicus_serotinus      --> eptesicus_serotinus should group with other species of Eptesicus. Force a monophyletic group for Vespertilio and check whether eptesicus_serotinus grouped with the other Eptesicus species.
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

In addition, there are three taxa (`dermanura_incomitata`, `artibeus_glaucus_watsoni`, and
`artibeus_glaucus_rosenbergii`) that are to be renamed as `artibeus_watsoni`. In order to ease future filtering, 
we are going to see which of these 3 taxa has more genes, and thus keep it. The other two, will be removed.
We are going to use the object `levels.checked.RData` so it is easier to track them: 

```r
# Open `parse_lineage.R` in RStudio, run lines from 1-15 and then load the R object 
# that was previously output (i.e., run line 67 after uncommenting it).
# The commands used to track the genes for each of these three taxa are the following:

load( "levels.checked.RData" )
levels.checked$genus$Dermanura$dermanura_incomitata
# [1] "CYB"             "ENSG00000165240" "ENSG00000175097" "ENSG00000176697"
levels.checked$genus$Artibeus$artibeus_glaucus_watsoni
# [1] "CYB"             "RNR2"            "ENSG00000165240" "ENSG00000175097" "ENSG00000176697"
levels.checked$genus$Artibeus$artibeus_glaucus_rosenbergii
# [1] "CYB"  "ND4"  "RNR2"
```

According to the output shown above, we are to keep only `artibeus_glaucus_watsoni` as `artibeus_watsoni`,
thus removing `dermanura_incomitata` and `artibeus_glaucus_rosenbergii`. 

The taxa to be renamed and/or removed according to this filtering are the following:
```
Taxa to rename:
artibeus_hartii,enchisthenes_hartii
artibeus_glaucus_watsoni,artibeus_watsoni
dermanura_rava,artibeus_rava
pipistrellus_subflavus,perimyotis_subflavus
plecotus_rafinesquii,corynorhinus_rafinesquii
eptesicus_nasutus_batinensis,rhyneptesicus_nasutus_batinensis
eptesicus_nasutus_matschiei,rhyneptesicus_nasutus_matschiei
eptesicus_nasutus,rhyneptesicus_nasutus
histiotus_macrotus,eptesicus_macrotus
histiotus_magellanicus,eptesicus_magellanicus
eptesicus_sagittula,vespadelus_darlingtoni
eptesicus_regulus,vespadelus_regulus
eptesicus_vulturnus,vespadelus_vulturnus
micronycteris_nicefori,trinycteris_nicefori
vampyressa_bidens,vampyriscus_bidens
vampyressa_brocki,vampyriscus_brocki
vampyressa_nymphaea,vampyriscus_nymphaea
emballonura_tiavato,paremballonura_tiavato
emballonura_atrata,paremballonura_atrata
sturnira_lilium_paulsoni,sturnira_paulsoni
myotis_daubentonii_nathalinae,myotis_nathalinae
artibeus_glaucus_gnomus,artibeus_gnomus
platyrrhinus_helleri_incarum,platyrrhinus_incarum
hipposideros_turpis_pendelburyi,hipposideros_pendleburyi
eumops_glaucinus_floridanus,eumops_floridanus
lophostoma_silvicolum_occidentalis,lophostoma_occidentalis
scotonycteris_zenkeri_occidentalis,scotonycteris_occidentalis
platyrrhinus_lineatus_nigellus,platyrrhinus_nigellus
nyctiellus_lepidus,nyctiellus_lepidus~
aproteles_bulmerae,aproteles_bulmerae~
pteropus_personatus,pteropus_personatus~
rhyneptesicus_nasutus_batinensis,rhyneptesicus_nasutus_batinensis~
rhyneptesicus_nasutus,rhyneptesicus_nasutus~
rhyneptesicus_nasutus_matschiei,rhyneptesicus_nasutus_matschiei~

Taxa to remove:
myotis_moluccarum
eptesicus_serotinus
eptesicus_dimissus
dermanura_incomitata
artibeus_glaucus_rosenbergii
```

**NOTE:**Data subsets (both alignments, 12CP and 3CP) for Afrotheria, Xenarthra, Euarchonta, and Marsupialia had already undergone 
these checks before 2018 (i.e., they had already been "cleaned", while the data subset for Lagomorpha 
had not yet).
Therefore, you do not see the last part of the filtering described in this section 
in the corresponding `README.md` files for these data subsets or the directory
[`taxonomical_check`](checked_aln/taxonomical_check)
You will see that the `filter_aln` directory for these four data subsets
contains a csv file with the details that were followed to filter the corresponding alignments (12CP and 3CP). 
Note that this csv file is equivalent to the excel sheet you find in this directory for 
data subset Chiroptera.

## 2. Check alignment with 3CP partition
Before we concatenate the alignment with the first and second codon positions (12CP), file `alignment.phylip`
(you will find it once you unzip the file which link is provided above),
and the alignment with the third codon
positions (3CP), file `alignment_nt3cp.phylip` (you will find it once 
you unzip the file which link provided above, inside a directory called `checked_aln` directory),
we ran the following code to make sure they had both undergone the same filtering steps and that were 
at the same "filtering stage":

```sh
# Run the next code from `chiroptera/filter_aln/checked_aln/unfiltered_aln`

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
   in RStudio and change line 25 so it is `subt    <- "chiroptera"` and uncomment line 27. Now, we
   can run it from RStudio. This script will generate a concatenated alignment file with all
   partitions, as well as one alignment file for each individual partitions, inside a new dir called `00_mammal_alns/chiroptera/unfiltered/`
   inside [`00_Data_filtering/01_alignments/`](../../../01_alignments). 
   Log and RData files will be saved inside
   [`00_Data_filtering/01_alignments/Rout`](../../../01_alignments/Rout/log_concatenation).   
      > NOTE 1: Updated `partitions.txt` file inside `00_mammal_alns/chiroptera/unfiltered/` generated.   
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
Run the following commands from `chiroptera/filter_aln/checked_aln`

```sh
# Run the next code from `laurasiatheria_chiroptera/checked_aln/checked_aln`

# 1. Fix species names in `names.txt`. We used the RAxML file, the `names.txt` file, and 
#    the `lineages.txt` that had been generated when filtering the alignment 12CP.
#    The ML tree file is now saved inside `filter_aln/checked_aln`, from where the 
#    python script will be run. The `names.txt` and the `lineage.txt` files are in
#    the `filter_aln` dir, one up. The alignment is the concatenated one that was generated in the 
#    previous step and can be found inside
#    `01_alignments/00_mammal_alns/chiroptera/unfiltered/chiroptera.aln`.
#    Note that we just copy this file inside `checked_aln` to proceed with the filtering. We also generate 
#    a copy of the original tree file:
cp ../../../../01_alignments/00_mammal_alns/chiroptera/unfiltered/chiroptera.aln .
mkdir original_tree 
cp RAxML_bestTree.BS_ML_GTRCAT original_tree/RAxML_bestTree.BS_ML_GTRCAT
cat > _species_name_change.txt << EOF
artibeus_hartii,enchisthenes_hartii
artibeus_glaucus_watsoni,artibeus_watsoni
dermanura_rava,artibeus_rava
pipistrellus_subflavus,perimyotis_subflavus
plecotus_rafinesquii,corynorhinus_rafinesquii
eptesicus_nasutus_batinensis,rhyneptesicus_nasutus_batinensis~
eptesicus_nasutus_matschiei,rhyneptesicus_nasutus_matschiei~
eptesicus_nasutus,rhyneptesicus_nasutus~
histiotus_macrotus,eptesicus_macrotus
histiotus_magellanicus,eptesicus_magellanicus
eptesicus_sagittula,vespadelus_darlingtoni
eptesicus_regulus,vespadelus_regulus
eptesicus_vulturnus,vespadelus_vulturnus
micronycteris_nicefori,trinycteris_nicefori
vampyressa_bidens,vampyriscus_bidens
vampyressa_brocki,vampyriscus_brocki
vampyressa_nymphaea,vampyriscus_nymphaea
emballonura_tiavato,paremballonura_tiavato
emballonura_atrata,paremballonura_atrata
sturnira_lilium_paulsoni,sturnira_paulsoni
myotis_daubentonii_nathalinae,myotis_nathalinae
artibeus_glaucus_gnomus,artibeus_gnomus
platyrrhinus_helleri_incarum,platyrrhinus_incarum
hipposideros_turpis_pendelburyi,hipposideros_pendleburyi
eumops_glaucinus_floridanus,eumops_floridanus
lophostoma_silvicolum_occidentalis,lophostoma_occidentalis
scotonycteris_zenkeri_occidentalis,scotonycteris_occidentalis
platyrrhinus_lineatus_nigellus,platyrrhinus_nigellus
nyctiellus_lepidus,nyctiellus_lepidus~
aproteles_bulmerae,aproteles_bulmerae~
pteropus_personatus,pteropus_personatus~
EOF

mkdir screen_logs
python ../../../../../../src/fix_species_names.py _species_name_change.txt RAxML_bestTree.BS_ML_GTRCAT chiroptera.aln ../names.txt ../lineage.txt > screen_logs/01_log1_fixspecies.txt

##>> NOTE: I needed to update the python script as it was not working. I just modifed `fix_alignment()` 
##>>       so the output file was different from the one being read

## Move unncessary output files to `00_filt1/` and rename output alignment
mkdir 00_filt1
mv _species_name_change.txt lineage.txt lineage.txt.1.bak 00_filt1/
rm chiroptera.aln
mv chiroptera_out.phylip alignment.phylip

# 2. Remove species 

# 2.1. Create input file for python script `prune_alignment.py`. This text file 
#      contains a list of the species that need to be removed 
cat > _species_remove.txt << EOF
rhinolophus_rufus
dermanura_incomitata
artibeus_glaucus_rosenbergii
miniopterus_schreibersii_schreibersii
myotis_daubentonii_daubentonii
myotis_blythii
myotis_adversus_adversus
myotis_horsfieldii_horsfieldii
eptesicus_serotinus_serotinus
plecotus_austriacus_austriacus
lasiurus_cinereus_cinereus
artibeus_planirostris_planirostris
artibeus_jamaicensis_jamaicensis
rhinolophus_darlingi_darlingi
myotis_moluccarum
eptesicus_serotinus
eptesicus_dimissus
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
   in RStudio and change line 24 so it is `subt    <- "chiroptera"` and uncomment line 26. 
  Now, we can run it from RStudio. 
   This script will generate a concatenated alignment file with all partitions, as well as 
   one alignment file for each individual partitions,
   inside [`01_alignments`](../../../01_alignments).
   Log files and Rdata can be found [here](../../../01_alignments/Rout/log_concatenation). 
      > NOTE: Updated `partitions.txt` file inside `00_mammal_alns/laurasiatheria_chiroptera/`.   
   * Paths have been automatically set according to current file architecture in `mammals` dir, 
   do not change paths in the Rscript! 


## 5. Extra checks for dubious taxon: Rhynopithecus sp.
Due to the long branches that Rhynopithecus taxa have in the tree, we have decided to extract the part of the
tree that contains the 
clade they belong to and run RAxML again to see if their placement changes. The analysis will be carried out inside 
`chiroptera/filter_aln/checked_aln/02_check_rhynopithecus/`.

First, we use the `RAxML_bestTree.BS_ML_GTRCAT` tree file to extract that specific clade. We copy this file and
manually add a label to identify the clade we want to keep it

```sh
# Run the next code from `checked_aln`
cd 02_check_rhynopithecus
cp ../RAxML_bestTree.BS_ML_GTRCAT .
cp ../../../../../01_alignments/00_mammal_alns/chiroptera/chiroptera.aln .
cp ../../../../../01_alignments/00_mammal_alns/chiroptera/partitions.txt .
```

Now that we have the tree, we open it with `FigTree` to root it and we manually:   

   * Root the tree having Monotremata as the outgroup.   
   * Order the nodes in increasing order.   
   * Extract the clade going from `antrozous_dubiaquercus` to `laephotis_wintoni`.   
   * Remove branch lengths and save it as `check_rhynopithecus.tree`.   
   
We can now run the R script `extract_remove_sp.R` to generate the `_species_remove.txt`:

```sh
# Remove tree without filtering as it was a copy 
rm RAxML_bestTree.BS_ML_GTRCAT

# If you want a list of the species in new pruned tree, uncomment and run the following, but we do 
# not use it in the R script
#    sed 's/(//g' check_rhynopithecus.tree | sed 's/)//g' | sed 's/;//'  | sed 's/,/\n/g' > present_species.txt 

# Run the R script from RStudio `extract_remove_sp.R`
```

Now, we can remove the taxa not included in this new tree from the alignment: 

```sh
python ../../../../../../../src/phy2fasta.py chiroptera.aln > alignment.fasta
python ../../../../../../../src/prune_alignment.py --remove _species_remove.txt alignment.fasta > alignment_pruned.fasta
python ../../../../../../../src/remove_gap_columns.py alignment_pruned.fasta > alignment_pruned_nogaps.fasta
## NOTE: You might see some numbers printed on the screen, do not worry, this is OK.
python ../../../../../../../src/fasta2phy.py alignment_pruned_nogaps.fasta > alignment_pruned_nogaps.phylip

# Clean dir
mkdir PY_outfiles 
mv *fasta _species_remove.txt PY_outfiles/
rm chiroptera.aln
mv alignment_pruned_nogaps.phylip alignment.phylip
```

Now, we can run `RAxML` with the same arguments as used in previous subtrees:

```sh
raxmlHPC -m GTRGAMMA -p 12345 -n test_rhynopithecus -s alignment.phylip -t check_rhynopithecus.tree
# cluster 
raxmlHPC-PTHREADS-AVX -T 8 -m GTRGAMMA -p 12345 -n test_rhynopithecus -s alignment.phylip -t check_rhynopithecus.tree
```

As expected, they do not cluster where they should, and the new placement still leads to quite long branches.
Therefore, we will proceed to delete these three taxa too:

```sh
# 1.0. Move previous files to `03_filt3` and copy the already filtered alignment 
# from 01_alignments/00_mammal_alns/chiroptera/
# Run this code from `checked_aln`
mkdir 03_filt3 
mv alignment.phylip *txt 03_filt3 
cp RAxML_bestTree.BS_ML_GTRCAT 03_filt3
cp ../../../../01_alignments/00_mammal_alns/chiroptera/chiroptera.aln . 

# 1.1. Create input file for python script `prune_alignment.py`. This text file 
#      contains a list of the species that need to be removed 
cat > _species_remove.txt << EOF
rhyneptesicus_nasutus~
rhyneptesicus_nasutus_batinensis~
rhyneptesicus_nasutus_matschiei~
EOF

# 2.2. Run set of four python scripts to generate pruned alignment without gaps:
python ../../../../../../src/phy2fasta.py chiroptera.aln > alignment.fasta
python ../../../../../../src/prune_alignment.py --remove _species_remove.txt alignment.fasta > alignment_pruned.fasta
python ../../../../../../src/remove_gap_columns.py alignment_pruned.fasta > alignment_pruned_nogaps.fasta
python ../../../../../../src/fasta2phy.py alignment_pruned_nogaps.fasta > alignment_pruned_nogaps.phylip

# 2.3. Now prune the tree.
#      All the names of taxa to be removed are saved in a variable and pass it as the second argument instead 
#      of using a text file to `nw_pruned`
sp2rm=$( cat _species_remove.txt )
../../../../../../src/newick-utils-1.6/src/nw_prune RAxML_bestTree.BS_ML_GTRCAT $sp2rm > tmp.tree #COOL!

##  Move unnecessary output files to `04_filt4/`. Rename tree and alignment files.
mkdir 04_filt4 
mv *fasta _species_remove.txt 04_filt4/
mv RAxML_bestTree.BS_ML_GTRCAT 04_filt4/RAxML_bestTree.BS_ML_GTRCAT_unpruned
mv tmp.tree RAxML_bestTree.BS_ML_GTRCAT
mv alignment_pruned_nogaps.phylip alignment.phylip
rm chiroptera.aln

# 3. Remove duplicates 

# 3.1. Find duplicates
grep -oP '\w+_duplicate' RAxML_bestTree.BS_ML_GTRCAT # No duplicates!

# 4. I generate a list of names to keep track of current species after filtering!
grep -o '[a-z].* ' alignment.phylip > names_filt.txt
```

## 6. Data partitioning 
Before partitioning the data, we need to save the already filtered files in the `00_mammal_alns/chiroptera/`
and label them as old:

```{sh}
# Run from 00_Data_filtering/01_alignments/00_mammal_alns/chiroptera/
mkdir chiroptera_old
mv *aln *txt chiroptera_old/	

cd ../../Rout/log_concatenation/
mkdir chiroptera_old
mv *chiroptera*txt chiroptera_old

cd ../Rdata 
mkdir chiroptera_old
mv *chiroptera.RData chiroptera_old/
```

Now that we have the concatenated alignment ready, we need to generate the filtered partitioned alignments by 
running the R script
[`Partition_seqs_for_MCMCtree_after_filtering.R`](../../../01_alignments/Partition_seqs_for_MCMCtree_after_filtering.R).

Instructions to follow:   

   * Open the RScript [`Partition_seqs_for_MCMCtree_after_filtering.R`](../../../01_alignments/Partition_seqs_for_MCMCtree_after_filtering.R) 
   in RStudio and change line 24 so it is `subt    <- "chiroptera"` and uncomment lines 26 and 51. 
   Now, we can run it from RStudio. 
   This script will generate a concatenated alignment file with all partitions, as well as 
   one alignment file for each individual partitions,
   inside [`01_alignments`](../../../01_alignments).
   Log files and Rdata can be found [here](../../../01_alignments/Rout/log_concatenation). 
      > NOTE: Updated `partitions.txt` file inside `00_mammal_alns/chiroptera/`.   
   * Paths have been automatically set according to current file architecture in `mammals` dir, 
   do not change paths in the Rscript! 

The final alignments generated at the end of this step can be downloaded from 
[here](https://www.dropbox.com/s/2u0ad87tmtvq143/SeqBayesS2_Raln_chiroptera.zip?dl=0).
They should be saved here if the same file architecture as the one set in the R scripts 
is to be used: `00_Data_filtering/01_alignments/00_mammals_alns/chiroptera`.

# EXTRA FILTERING -- DATA SUBSETTING
When we first carried out the analysis with the tree topology as in step 1, we realised that 
it was too big to analyse with `MCMCtree` (i.e., too many taxa). Therefore, we decided to 
further subset this data and generate two data subsets: "chiroptera subtree 1" and 
"chiroptera subtree 2". 

The procedure followed was the following:

## 1. Explore partitioning the data set
The directory
[`00_R_parsing`](extra_filtering/00_R_parsing)
has the input/output files used to explore how to partition the data set into two data subsets. Please access this directory using the link 
provided to go through the steps followed. More details about this step can be found 
[here](extra_filtering),
in the first section `1. Obtain subtrees`.

## 2. Generate alignments 
[Here](extra_filtering), 
in section `2. Generating alignments`, 
you will read about all the steps we followed to extract the sequences of the taxa that need to be allocated to each data 
subset. The data can be downloaded 
[here](https://www.dropbox.com/s/coyqiu6d7owvtho/SeqBayesS2_filteraln2_chiroptera_01_perl_parsing.zip?dl=0)
if you want to check you have reproduced our results, which should be saved 
in the directory
[`01_perl_parsing`](extra_filtering/01_perl_parsing).
In this GitHub repository, due to limited space for big files, you will see only the Perl script we use 
to generate the alignments.

## 3. Add new taxa and generate alignments
In order to avoid issues when grafting the subtrees to the backbone tree, we decided to add 
extra taxa to each subtree (one species from data subset 1 is included in data subset 2, and 
viceversa). The steps followed to add taxa to the first subtree can be 
found
[here](extra_filtering/02_MAFFT_subt1),
while those followed to generate the one for the second subtree can be found 
[here](extra_filtering/02_MAFFT_subt2).

## 4. Final alignments 
[Here](https://www.dropbox.com/s/h5mqaqd8tabkww8/SeqBayesS2_Raln_chiroptera_subt1.zip?dl=0)
you can download the final alignments for subtree 1 
and
[here](https://www.dropbox.com/s/5zwds5t26r5l29f/SeqBayesS2_Raln_chiroptera_subt2.zip?dl=0)
for subtree 2. 
