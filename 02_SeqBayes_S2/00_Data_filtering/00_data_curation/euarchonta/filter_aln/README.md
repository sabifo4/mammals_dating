# EUARCHONTA - Filtering alignment
The initial files that you can find in this directory are the following:

```
euarchonta
   |- filter_aln
         |- euarchonta_taxonomy_check.csv 
         |- lineage.txt 
         |- names.txt 
         |- parse_lineage.R 
         |- partitions.txt 
         |- README.md
         |- summary.html 
```
		 
The output files during the filtering step that is described below as well as the input 
files needed have been zipped in a file as they are very large. You can download 
this file [here](https://www.dropbox.com/s/xbjm9qt1j194t70/SeqBayesS2_filtaln_euarchonta.zip?dl=0).
To start the filtering step, you should have the following files here, which 
you can obtain once you unzip the file provided above: 

```
euarchonta 
   |- filter_aln
         |- checked_aln                        
         |     |- alignment_nt3cp.phylip      # Alignment with 3CP (found in `original_aln_tree`)
         |     |- RAxML_bestTree.BS_ML_GTRCAT # Unfiltered tree (found in `original_aln_tree`)
         |     
         |- alignment.phylip # Alignment with 12CP
         |- euarchonta_taxonomy_check.csv 
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
FAMILY   Aotidae  has  12  taxa. Species  aotus_zonalis  has  3  genes that are not shared by any of the other taxa.
It should be deleted!
```
>> ACTION: aotus_zoalis was removed!

```
GENUS   Phaner  has  2  taxa. Species  phaner_furcifer  has  1  genes that are not shared by any of the other taxa.
It should be deleted!
GENUS   Phaner  has  2  taxa. Species  phaner_pallescens  has  1  genes that are not shared by any of the other taxa.
It should be deleted!
```
>> ACTION: these two species were merged into one, which tag name was `phaner` to represent the species.

```
GENUS   Aotus  has  12  taxa. Species  aotus_zonalis  has  3  genes that are not shared by any of the other taxa.
It should be deleted!
```
>> ACTION: same finding as in family: aotus_zonalis already removed.

```
GENUS   Miopithecus  has  3  taxa. Species  miopithecus_talapoin_talapoin  has  2  genes that are not shared
by any of the other taxa. It should be deleted!```
>> ACTION: miopithecus_talapoin_talapoin was removed.

```
GENUS   Presbytis  has  15  taxa. Species  presbytis_senex  has  1  genes that are not shared by any of the other taxa.
It should be deleted!
```
>> ACTION: presbytis_senex was removed.   


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
[`euarchonta_taxonomy_check.csv`](euarchonta_taxonomy_check.csv). 
The filtering you will see in the commands below was already 
carried out for `alignment.phylip`, but has had to be done again for `alignment_nt3cp.phylip` so both 
alignment files are at the same filtering level. The code is the following:

```sh
# Run the next code from `euarchonta/filter_aln/checked_aln`

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
callicebus_brunneus,plecturocebus_brunneus
callicebus_cupreus,plecturocebus_cupreus
callicebus_donacophilus,plecturocebus_donacophilus
callicebus_hoffmannsi,plecturocebus_hoffmannsi
callicebus_lugens,cheracebus_lugens
callicebus_moloch,plecturocebus_moloch
callicebus_purinus,cheracebus_purinus
callicebus_torquatus,cheracebus_torquatus
galago_granti,galagoides_granti
alouatta_seniculus_macconnelli,alouatta_macconnelli
sapajus_apella_macrocephalus,sapajus_macrocephalus
cercopithecus_albogularis_moloneyi,cercopithecus_mitis_moloneyi
cercopithecus_albogularis_labiatus,cercopithecus_mitis_labiatus
cercopithecus_albogularis_erythrarchus,cercopithecus_mitis_erythrarchus
cercopithecus_albogularis,cercopithecus_mitis_duplicate
cercopithecus_albogularis_monoides,cercopithecus_mitis_monoides
cercopithecus_albogularis_albotorquatus,cercopithecus_mitis_albotorquatus
cercopithecus_albogularis_francescae,cercopithecus_mitis_francescae
cercopithecus_albogularis_kolbi,cercopithecus_mitis_kolbi
cercopithecus_doggetti,cercopithecus_mitis_doggetti
cercopithecus_kandti,cercopithecus_mitis_kandti
propithecus_deckenii_coronatus,propithecus_coronatus
sapajus_nigritus_robustus,sapajus_robustus
loris_tardigradus_nordicus,loris_lydekkerianus_nordicus
eulemur_fulvus_albocollaris,eulemur_cinereiceps_duplicate
eulemur_fulvus_collaris,eulemur_collaris
EOF

mkdir screen_logs
python ../../../../../../src/fix_species_names.py _species_name_change.txt RAxML_bestTree.BS_ML_GTRCAT alignment_nt3cp.phylip ../names.txt ../lineage.txt > screen_logs/01_log1_fixspecies.txt

##>> NOTE: We had to update the python script as it was not working with this alignment.
##         Modification is noted down in the python script, which involves `fix_alignment()` 
##>>       so the output file was different from the one being read

## Move unncessary output files to 00_filt1 and rename output alignment
mkdir 00_filt1
mv _species_name_change.txt lineage.txt lineage.txt.1.bak 00_filt1/
rm alignment_nt3cp.phylip
mv alignment_nt3cp_out.phylip alignment_nt3cp.phylip

# 2. Remove species 

# 2.1. Create input file for python script `prune_alignment.py`. This text file 
#    contains a list of the species that need to be removed 
cat > _species_remove.txt << EOF
lophocebus_albigena
saimiri_sciureus
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
python ../../../../../../src/phy2fasta.py alignment_nt3cp.phylip > alignment_nt3cp.fasta
python ../../../../../../src/prune_alignment.py --remove _species_remove.txt alignment_nt3cp.fasta > alignment_nt3cp_pruned.fasta
python ../../../../../../src/remove_gap_columns.py alignment_nt3cp_pruned.fasta > alignment_nt3cp_pruned_nogaps.fasta
python ../../../../../../src/fasta2phy.py alignment_nt3cp_pruned_nogaps.fasta > alignment_nt3cp_pruned_nogaps.phylip

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
mv *fasta alignment_nt3cp.phylip _species_remove.txt 01_filt2/
mv RAxML_bestTree.BS_ML_GTRCAT 01_filt2/RAxML_bestTree.BS_ML_GTRCAT_unpruned
mv tmp.tree RAxML_bestTree.BS_ML_GTRCAT
mv alignment_nt3cp_pruned_nogaps.phylip alignment_nt3cp.phylip

# 3. Remove duplicates 

# 3.1. Find duplicates
grep -oP '\w+_duplicate' RAxML_bestTree.BS_ML_GTRCAT
## cercopithecus_mitis_duplicate
## eulemur_cinereiceps_duplicate

# 3.2. Create input file for python script `prune_alignment.py`. This text file 
#    contains a list of the species that need to be removed 
#    I havo also added "macaca_nemestrina_nemestrina" as it is the last one 
#    removed by Asif because there are only 2 genes
cat > _species_remove_duplicates.txt << EOF
cercopithecus_mitis_duplicate
eulemur_cinereiceps_duplicate
macaca_nemestrina_nemestrina
EOF

# 3.3. Run set of python scripts as in step 2.2
python ../../../../../../src/phy2fasta.py alignment_nt3cp.phylip > alignment_nt3cp.fasta
python ../../../../../../src/prune_alignment.py --remove _species_remove_duplicates.txt alignment_nt3cp.fasta > alignment_nt3cp_pruned.fasta
python ../../../../../../src/remove_gap_columns.py alignment_nt3cp_pruned.fasta > alignment_nt3cp_pruned_nogaps.fasta
python ../../../../../../src/fasta2phy.py alignment_nt3cp_pruned_nogaps.fasta > alignment_nt3cp_pruned_nogaps.phylip

# 3.4. Now prune the tree as in step 2.3
sp2rm=$( cat _species_remove_duplicates.txt )
../../../../../../src/newick-utils-1.6/src/nw_prune RAxML_bestTree.BS_ML_GTRCAT $sp2rm > tmp.tree #COOL!

##  Move unnecessary output files to 02_filt3. Rename tree and alignment files.
mkdir 02_filt3 
mv *fasta alignment_nt3cp.phylip _species_remove_duplicates.txt 02_filt3
mv RAxML_bestTree.BS_ML_GTRCAT 02_filt3/RAxML_bestTree.BS_ML_GTRCAT_unpruned
mv tmp.tree RAxML_bestTree.BS_ML_GTRCAT
mv alignment_nt3cp_pruned_nogaps.phylip alignment_nt3cp.phylip

# 4. Generate a list of names to check that they are the same 
grep -o '[a-z].* ' alignment_nt3cp.phylip > names_3nt.txt
grep -o '[a-z].* ' ../alignment.phylip > ../names_aln.txt
diff names_3nt.txt ../names_aln.txt # OK!
````

Now that both files with the 12CP-alignment (`alignment.phylip`) and 3CP-alignment (`alignment_nt3cp.phylip`) 
have the same taxa names and have had the same filtering procedure carried out,
we can proceed to concatenate the alignment with the third 
codon positions to the alignment with 12CPs using the R script
[`Concatenate_seqs_for_MCMCtree.R`](../../../01_alignments/Concatenate_seqs_for_MCMCtree.R)
so we can have the concatenated alignment. 

Instructions to follow:    

   * Open the RScript [`Concatenate_seqs_for_MCMCtree.R`](../../../01_alignments/Concatenate_seqs_for_MCMCtree.R),
   change line 24 so it is `subt    <- "euarchonta"`, and run it from RStudio. This script
   will generate a concatenated alignment file with all partitions, as well as 
   one alignment file for each individual partitions, inside a new dir called `00_mammal_alns/euarchonta`
   inside [`00_Data_filtering/01_alignments/`](../../../01_alignments). 
   Log and RData files will be saved inside 
   [`00_Data_filtering/01_alignments/Rout`](../../../01_alignments/Rout/log_concatenation). 
      > NOTE 1: Updated `partitions.txt` file inside `00_mammal_alns/euarchonta` generated.   
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
## There are 51 subspecies with species
```
The output names have been saved in a file called `subsp_check.txt` to further explore it:

```sh
# Run from `00_data_curation/euarchonta/filter_aln/checked_aln` directory 
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

> **Remove species, less genes**   
   Num. sequences **cercopithecus_nictitans** : 3   
   Num. sequences cercopithecus_nictitans_nictitans : 20

> **Remove subspecies, less genes**    
   Num. sequences cercopithecus_mitis : 33   
   Num. sequences **cercopithecus_mitis_mitis** : 20

> **Remove species, less genes**   
   Num. sequences **cercopithecus_preussi** : 2   
   Num. sequences cercopithecus_preussi_preussi : 20

> Species had already been removed   
   ~~Num. sequences cercopithecus_petaurista : 2~~   
   Num. sequences cercopithecus_petaurista_petaurista : 20

> **Remove subspecies, less genes**   
   Num. sequences cercopithecus_cephus : 22   
   Num. sequences **cercopithecus_cephus_cephus** : 12

> Species had already been removed   
   ~~Num. sequences cercopithecus_ascanius : 2~~   
   Num. sequences cercopithecus_ascanius_ascanius : 20

> Cannot be removed   
   Num. sequences cercopithecus_erythrotis : 0   
   Num. sequences cercopithecus_erythrotis_erythrotis : 20

> Cannot be removed   
   Num. sequences cercopithecus_pogonias : 0   
   Num. sequences cercopithecus_pogonias_pogonias : 20

> **Remove subspecies, less genes**   
   Num. sequences papio_ursinus : 27   
   Num. sequences **papio_ursinus_ursinus** : 2

> **Remove subspecies, less genes**  
   Num. sequences papio_cynocephalus : 28   
   Num. sequences **papio_cynocephalus_cynocephalus** : 2

> **Remove subspecies, less genes**   
   Num. sequences papio_hamadryas : 54   
   Num. sequences **papio_hamadryas_hamadryas** : 2

> Species had already been removed   
   Num. sequences macaca_nemestrina : 47   
   ~~Num. sequences macaca_nemestrina_nemestrina : 2~~

> **Remove subspecies, less genes**  
   Num. sequences macaca_fascicularis : 85   
   Num. sequences **macaca_fascicularis_fascicularis** : 1

> **Remove subspecies, less genes**   
   Num. sequences hylobates_agilis : 48   
   Num. sequences **hylobates_agilis_agilis** : 2

> **Remove subspecies, less genes**   
   Num. sequences hylobates_muelleri : 24   
   Num. sequences **hylobates_muelleri_muelleri** : 2

> **Remove subspecies, less genes**   
   Num. sequences hylobates_lar : 60   
   Num. sequences **hylobates_lar_lar** : 2

> **Remove subspecies, less genes**   
   Num. sequences nomascus_concolor : 29   
   Num. sequences **nomascus_concolor_concolor** : 2

> **Remove subspecies, less genes**   
   Num. sequences pongo_pygmaeus : 84   
   Num. sequences **pongo_pygmaeus_pygmaeus** : 3

> **Remove subspecies, less genes**   
   Num. sequences gorilla_gorilla : 186   
   Num. sequences **gorilla_gorilla_gorilla** : 28

> **Remove subspecies, less genes**   
   Num. sequences pan_troglodytes : 189   
   Num. sequences **pan_troglodytes_troglodytes** : 27

> **Remove subspecies, less genes**   
   Num. sequences galeopterus_variegatus : 41   
   Num. sequences **galeopterus_variegatus_variegatus** : 2

> **Remove subspecies, less genes**   
   Num. sequences loris_lydekkerianus : 26   
   Num. sequences **loris_lydekkerianus_lydekkerianus** : 2

> **Remove subspecies, less genes**   
   Num. sequences loris_tardigradus : 31   
   Num. sequences **loris_tardigradus_tardigradus** : 2

> **Remove subspecies, less genes**   
   Num. sequences galago_senegalensis : 33   
   Num. sequences **galago_senegalensis_senegalensis** : 2

> Species had already been removed   
   ~~Num. sequences galagoides_zanzibaricus : 3~~   
   Num. sequences galagoides_zanzibaricus_zanzibaricus : 2

> Species had already been removed   
   ~~Num. sequences otolemur_monteiri : 2~~   
   Num. sequences otolemur_monteiri_monteiri : 1   
   ```
   otolemur_monteiri_monteiri --> "RNR1"   
   otolemur_monteiri          --> "CYB"
   ```
   
> **Remove species, less genes**   
   Num. sequences **varecia_variegata** : 27   
   Num. sequences varecia_variegata_variegata : 31

> **Remove subspecies, less genes**   
   Num. sequences hapalemur_griseus : 29   
   Num. sequences **hapalemur_griseus_griseus** : 7

> **Remove species, less genes**   
   Num. sequences **eulemur_macaco** : 6   
   Num. sequences eulemur_macaco_macaco : 28

> **Remove species, less genes**   
   Num. sequences **eulemur_fulvus** : 6   
   Num. sequences eulemur_fulvus_fulvus : 30

> **Remove subspecies, priority to CYTB**   
   Num. sequences propithecus_deckenii : 2   
   Num. sequences **propithecus_deckenii_deckenii** : 2   
   ```
   propithecus_deckenii_deckenii --> "ND4"
   propithecus_deckenii          --> "CYB"
   ```

> **Remove subspecies, less genes**   
   Num. sequences propithecus_verreauxi : 40   
   Num. sequences **propithecus_verreauxi_verreauxi** : 5

> **Remove subspecies, less genes**   
   Num. sequences propithecus_diadema : 26   
   Num. sequences **propithecus_diadema_diadema** : 4

> **Remove subspecies, less genes**   
   Num. sequences cacajao_calvus : 31   
   Num. sequences **cacajao_calvus_calvus** : 2

> **Remove subspecies, less genes**   
   Num. sequences aotus_azarai : 31   
   Num. sequences **aotus_azarai_azarai** : 26

> Species had already been removed   
   ~~Num. sequences saimiri_sciureus : 62~~   
   Num. sequences saimiri_sciureus_sciureus : 22

> **Remove species, less genes**   
   Num. sequences **saimiri_oerstedii** : 7   
   Num. sequences saimiri_oerstedii_oerstedii : 20

> **Remove subspecies, less genes**   
   Num. sequences saimiri_boliviensis : 38   
   Num. sequences **saimiri_boliviensis_boliviensis** : 28

> Cannot be removed   
   Num. sequences saguinus_nigricollis : 0   
   Num. sequences saguinus_nigricollis_nigricollis : 2

> Species had already been removed   
   ~~Num. sequences saguinus_fuscicollis : 9~~   
   Num. sequences saguinus_fuscicollis_fuscicollis : 2   
   ```
   saguinus_fuscicollis_fuscicollis --> "CYB"
   saguinus_fuscicollis             --> "ENSG00000066279" "ENSG00000110243" "ENSG00000110244" "ENSG00000110799" "ENSG00000112964" "ENSG00000174951" "ENSG00000187554" "ENSG00000257138" "ENSG00000265203"
   ```
   
> **Remove subspecies, less genes**   
   Num. sequences alouatta_palliata : 8   
   Num. sequences **alouatta_palliata_palliata** : 2

> **Remove subspecies, less genes**   
   Num. sequences alouatta_seniculus : 31   
   Num. sequences **alouatta_seniculus_seniculus** : 3

> Species had already been removed   
   Num. sequences colobus_satanas : 0   
   Num. sequences colobus_satanas_satanas : 4

> **Remove subspecies, less genes**   
   Num. sequences piliocolobus_badius : 26   
   Num. sequences **piliocolobus_badius_badius** : 4

> **Remove subspecies, less genes**   
   Num. sequences trachypithecus_poliocephalus : 3   
   Num. sequences **trachypithecus_poliocephalus_poliocephalus** : 2   
   ```
   trachypithecus_poliocephalus_poliocephalus --> "ND4"
   trachypithecus_poliocephalus               --> "ND4" "ENSG00000077498"
   ```

> **Remove subspecies, less genes**   
   Num. sequences trachypithecus_francoisi : 56   
   Num. sequences **trachypithecus_francoisi_francoisi** : 2

> Cannot be removed   
   Num. sequences presbytis_hosei : 0   
   Num. sequences presbytis_hosei_hosei : 2

> Cannot be removed   
   Num. sequences presbytis_comata : 0   
   Num. sequences presbytis_comata_comata : 2

> **Remove subspecies, less genes**   
   Num. sequences presbytis_melalophos : 30   
   Num. sequences **presbytis_melalophos_melalophos** : 2

> Cannot be removed   
   Num. sequences presbytis_rubicunda : 0   
   Num. sequences presbytis_rubicunda_rubicunda : 2

> Cannot be removed   
   Num. sequences presbytis_potenziani : 0   
   Num. sequences presbytis_potenziani_potenziani : 2

> **Remove subspecies, less genes**   
   Num. sequences pygathrix_nemaeus : 41   
   Num. sequences **pygathrix_nemaeus_nemaeus** : 3

As this filtering was not applied before we generate the alignment, we are going to use the already
concatenated alignment with filtered `3nt` and 
`1st+2nd` CPs and remove these taxa. For that purpose, we are going to run the next commands: 

```sh
# Run from `filter_aln/checked_aln` directory 

# 1. Copy filtered alignment in `00_Data_filtering/01_alignments/00_mammal_alns/euarchonta/` 
#    called `euarchonta.aln` in the `filter_aln/checked_aln` directory. Remember that this is 
#    the concatenated filtered alignment that we previously generated. Once there, change the 
#    file name: 
cp ../../../../01_alignments/00_mammal_alns/euarchonta/euarchonta.aln .
mv euarchonta.aln alignment.phylip

# 2. Remove species 

# 2.1. Create input file for python script `prune_alignment.py`. This text file 
#    contains a list of the species that need to be removed 
cat > _species_remove.txt << EOF
cercopithecus_nictitans
cercopithecus_mitis_mitis
cercopithecus_preussi
cercopithecus_cephus_cephus
papio_ursinus_ursinus
papio_cynocephalus_cynocephalus
papio_hamadryas_hamadryas
macaca_fascicularis_fascicularis
hylobates_agilis_agilis
hylobates_muelleri_muelleri
hylobates_lar_lar
nomascus_concolor_concolor
pongo_pygmaeus_pygmaeus
gorilla_gorilla_gorilla
pan_troglodytes_troglodytes
galeopterus_variegatus_variegatus
loris_lydekkerianus_lydekkerianus
loris_tardigradus_tardigradus
varecia_variegata
hapalemur_griseus_griseus
eulemur_macaco
eulemur_fulvus
propithecus_deckenii_deckenii
propithecus_verreauxi_verreauxi
propithecus_diadema_diadema
cacajao_calvus_calvus
aotus_azarai_azarai
saimiri_oerstedii
saimiri_boliviensis_boliviensis
alouatta_palliata_palliata
alouatta_seniculus_seniculus
piliocolobus_badius_badius
trachypithecus_poliocephalus_poliocephalus
trachypithecus_francoisi_francoisi
presbytis_melalophos_melalophos
pygathrix_nemaeus_nemaeus
EOF

# 2.2. Run set of four python scripts that will
python ../../../../../../src/phy2fasta.py alignment.phylip > alignment.fasta
python ../../../../../../src/prune_alignment.py --remove _species_remove.txt alignment.fasta > alignment_pruned.fasta
python ../../../../../../src/remove_gap_columns.py alignment_pruned.fasta > alignment_pruned_nogaps.fasta
python ../../../../../../src/fasta2phy.py alignment_pruned_nogaps.fasta > alignment_pruned_nogaps.phylip

# 2.3. Now prune the tree.
sp2rm=$( cat _species_remove.txt )
../../../../../../src/newick-utils-1.6/src/nw_prune RAxML_bestTree.BS_ML_GTRCAT $sp2rm > tmp.tree #COOL!

##  Move unnecessary output files to 03_filt4. Rename tree and alignment files.
mkdir 03_filt4 
mv *fasta alignment.phylip _species_remove.txt 03_filt4/
mv RAxML_bestTree.BS_ML_GTRCAT 03_filt4/RAxML_bestTree.BS_ML_GTRCAT_unpruned
mv tmp.tree RAxML_bestTree.BS_ML_GTRCAT
mv alignment_pruned_nogaps.phylip alignment.phylip

# 3. Remove duplicates 

# 3.1. Find duplicates
grep -oP '\w+_duplicate' RAxML_bestTree.BS_ML_GTRCAT # No duplicates!

# 4. Generate a list of names with the taxa in the updated filtered alignment 
grep -o '[a-z].* ' alignment.phylip > names_3nt_filt.txt
```

## 4. Data partitioning 
Now that we have the concatenated and filtered alignment ready, we need to generate the filtered
partitioned alignments by running the R script 
[`Partition_seqs_for_MCMCtree_after_filtering.R`](../../../01_alignments/Partition_seqs_for_MCMCtree_after_filtering.R).
As we had already created partitioned alignments 
before the subspecies check had been applied, we will need to rearrange the files output for
this data subset in the
[`00_Data_filtering/01_alignments/`](../../../01_alignments)
subdirectories:

```sh
# Run from 00_Data_filtering/01_alignments 

# Update alignemnts dir 
cd 00_mammal_alns/euarchonta 
mkdir euarchonta_old 
mv *aln *txt euarchonta_old 

# Update Rout dir 
cd ../../Rout/log_concatenation/
mkdir euarchonta_old 
mv *euarchonta*txt euarchonta_old
cd ../Rdata/
mkdir euarchonta_old
mv *euarchonta*txt *euarchonta*RData euarchonta_old
```

Now, run the Rscript mentioned above following the next instructions:   

   * Open the RScript [`Partition_seqs_for_MCMCtree_after_filtering.R`](../../../01_alignments/Partition_seqs_for_MCMCtree_after_filtering.R) 
   in RStudio and change line 24 so it is `subt    <- "euarchonta"` and uncomment line 30. Now, we can run it from RStudio. 
   This script will generate a concatenated alignment file with all partitions, as well as 
   one alignment file for each individual partitions, inside the directory `00_mammal_alns/euarchonta/`
   inside [`01_alignments`](../../../01_alignments).
   Log files and Rdata can be found [here](../../../01_alignments/Rout/log_concatenation). 
      > NOTE: Updated `partitions.txt` file inside `00_mammal_alns/euarchonta/`.   
   * Paths have been automatically set according to current file architecture in `mammals` dir, 
   do not change paths in the Rscript! 
   
The final alignments generated at the end of this step can be downloaded from 
[here](https://www.dropbox.com/s/djlncitmsehsyap/SeqBayesS2_Raln_euarchonta.zip?dl=0).
They should be saved here if the same file architecture as the one set in the R scripts 
is to be used: `00_Data_filtering/01_alignments/00_mammals_alns/euarchonta`.
