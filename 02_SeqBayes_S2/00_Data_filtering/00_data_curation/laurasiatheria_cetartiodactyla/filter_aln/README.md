# LAURASIATHERIA CETARTIODACTYLA (ARTIODACTYLA) - Filtering alignment
The initial files that you can find in this directory are the following:

```
laurasiatheria_cetartiodactyla
   |- filter_aln
         |- checked_aln                        
         |     |- taxonomical_check
         |           |- laurasiatheria_cetartiodactyla.png
         |           |- laurasiatheria_cetartiodactyla_FigTree
         |           |- RAxML_bestTree.BS_ML_GTRCAT
         |     
         |- L.cetartiodactyla_taxonomy_check.xlsx
         |- lineage.txt 
         |- names.txt 
         |- parse_lineage.R 
         |- partitions.txt 
         |- README.md
         |- summary.html 
```
		 
The output files during the filtering step that is described below as well as the input 
files needed have been zipped in a file as they are very large. You can download 
this file [here](https://www.dropbox.com/s/4szw979qocelxf5/SeqBayesS2_filtaln_artiodactyla.zip?dl=0).
To start the filtering step, you should have the following files arranged in the file 
architecture detailed above (you can obtain the files once you unzip the file
provided in the link above): 

```
laurasiatheria_cetartiodactyla 
   |- filter_aln
         |- checked_aln                        
         |     |- unfiltered_aln
         |     |    |- alignment_nt3cp.phylip      # Alignment with 3CP 
         |     |- RAxML_bestTree.BS_ML_GTRCAT      # Unfiltered tree (found in `original_tree`)
         |     |- taxonomical_check                           # Visual checks to evaluate 
         |          |- laurasiatheria_cetartiodactyla.png     # dubious taxa placement
         |          |- laurasiatheria_cetartiodactyla_FigTree
         |          |- RAxML_bestTree.BS_ML_GTRCAT    
         |          
         |- alignment.phylip # Alignment with 12CP
         |- L.cetartiodactyla_taxonomy_check.xlsx
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
Within the [`filter_aln`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/laurasiatheria_cetartiodactyla/filter_aln)
directory, you will find different files as detailed above. 
Specifically, the R script [`parse_lineage.R`](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/laurasiatheria_cetartiodactyla/filter_aln/parse_lineage.R)
was written to carry out a first taxonomic filtering. Note that this 
R script will run if you have the same file architecture in this GitHub repository (i.e., it uses 
a function within the R script [`Filter_lineages.R`](https://github.com/sabifo4/mammals_dating/blob/main/src/Filter_lineages.R)
in the `src` directory and the [`genes.txt`](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/genes.txt)
file. The messages printed out by this script are the following:   

```
ORDER   Suina  has  23  taxa. Species  sus_scrofa_ussuricus  has  1  genes that are not shared by any of the other taxa. It should be deleted!```
>> ACTION: sus_scrofa_ussuricus removed.

```
FAMILY   Suidae  has  20  taxa. Species  sus_scrofa_ussuricus  has  1  genes that are not share by any of the other taxa. It should be deleted!
```
>> ACTION: Same finding as in order: sus_scrofa_ussuricus already removed.

```
SUBFAMILY   Alcelaphinae  has  9  taxa. Species  connochaetes_taurinus_taurinus  has  1  genes that are not share by any of the other taxa. It should be deleted!
SUBFAMILY   Alcelaphinae  has  9  taxa. Species  damaliscus_pygargus_phillipsi  has  1  genes that are not share by any of the other taxa. It should be deleted!
```
>> ACTION: Both species connochaetes_taurinus_taurinus & damaliscus_pygargus_phillipsi removed.

```
GENUS   Sus  has  13  taxa. Species  sus_scrofa_ussuricus  has  1  genes that are not share by any of the other taxa. It should be deleted!
```
>> ACTION: same finding as in order: sus_scrofa_ussuricus already removed.

```
GENUS   Syncerus  has  3  taxa. Species  syncerus_caffer  has  24  genes that are not share by any of the other taxa. It should be deleted!
GENUS   Connochaetes  has  3  taxa. Species  connochaetes_taurinus_taurinus  has  1  genes that are not share by any of the other taxa. It should be deleted!
```
>> ACTION: the other two species are subspecies of syncerus_caffer. They happen to share the same gene,
>> which is not found in syncerus_caffer -- we keep it. We remove the two subspecies.

```
GENUS   Damaliscus  has  3  taxa. Species  damaliscus_pygargus_phillipsi  has  1  genes that are not share by any of the other taxa. It should be deleted!
```
>> ACTION: Same finding as in subfamily: connochaetes_taurinus_taurinus & damaliscus_pygargus_phillipsi were removed.

**NOTE**
Log files will be found in this same directory with the name `log_taxaNOTIN_<level>.txt`,
being `level` the one that has been checked (i.e., family, genus, order, or subfamily). In addition, 
you can also generate the `levels.checked.RData`, which you can use to explore the taxonomical 
levels. 

# 2. First checks before applying filtering
First, we checked if there were any further species that should be removed:

```sh
# Run from `filter_aln/checked_aln` directory 
grep -o '[a-z]*_[a-z]*_[a-z]*' RAxML_bestTree.BS_ML_GTRCAT | sed 's/\_/\t/g' | awk '$2==$3{print $1,$2,$3}' | awk '$1!="pdh"{print $1,$2,$3}' > subsp_check.txt

# Count 
num=$( grep -o '[a-z]*_[a-z]*_[a-z]*' RAxML_bestTree.BS_ML_GTRCAT | sed 's/\_/\t/g' | awk '$2==$3{print $1,$2,$3}' | awk '$1!="pdh"{print $1,$2,$3}'| wc -l )
printf "There are $num subspecies with species\n"
## There are 26 subspecies with species
```

The output names have been saved in a file called `subsp_check.txt` to further explore it:

```sh
# Run from `00_data_curation/laurasiatheria_cetartiodactyla/filter_aln/checked_aln` directory 
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
   Num. sequences rangifer_tarandus : 29   
   Num. sequences **rangifer_tarandus_tarandus** : 3

> **REMOVE subspecies, less genes**   
   Num. sequences alces_alces : 27   
   Num. sequences **alces_alces_alces** : 1

> **REMOVE subspecies, less genes**   
   Num. sequences odocoileus_hemionus : 28   
   Num. sequences **odocoileus_hemionus_hemionus** : 4

> **REMOVE subspecies, less genes**   
   Num. sequences rucervus_eldii : 26   
   Num. sequences **rucervus_eldii_eldii** : 1

> **REMOVE subspecies, less genes**   
   Num. sequences axis_porcinus : 26   
   Num. sequences **axis_porcinus_porcinus** : 2

> **REMOVE subspecies, less genes**   
   Num. sequences cervus_elaphus : 40   
   Num. sequences **cervus_elaphus_elaphus** : 3

> No action   
   Num. sequences cervus_canadensis : 0   
   Num. sequences cervus_canadensis_canadensis : 3

> **REMOVE subspecies, less genes**   
   Num. sequences lama_guanicoe : 30   
   Num. sequences **lama_guanicoe_guanicoe** : 2

> **REMOVE subspecies, less genes**   
   Num. sequences sus_barbatus : 26   
   Num. sequences **sus_barbatus_barbatus** : 5

> **REMOVE subspecies, less genes**   
   Num. sequences sus_scrofa : 165   
   Num. sequences **sus_scrofa_scrofa** : 20

> **REMOVE subspecies, less genes**   
   Num. sequences inia_geoffrensis : 41   
   Num. sequences **inia_geoffrensis_geoffrensis** : 3

> **REMOVE subspecies, less genes**   
   Num. sequences neophocaena_asiaeorientalis : 29   
   Num. sequences **neophocaena_asiaeorientalis_asiaeorientalis** : 20   
   ```
   # neophocaena_asiaeorientalis_asiaeorientalis --> "ATP6" "CO1"  "CO2"  "CO3"  "CYB"  "ND1"  "ND2"  "ND4"  "ND5"  "RNR1" "RNR2"   
   # neophocaena_asiaeorientalis                 --> "ATP6"  "ATP8" "CO1" "CO2" "CO3" "CYB"  "ND1" "ND2"  "ND3" "ND4" "ND4L" "ND5"  "ENSG00000100644" "ENSG00000171885" "ENSG00000174697" "RNR1" "RNR2"              
   ```

> **REMOVE species, less genes**   
   Num. sequences **cephalorhynchus_hectori** : 3   
   Num. sequences cephalorhynchus_hectori_hectori : 18

> **REMOVE subspecies, less genes**   
   Num. sequences bubalus_bubalis : 62   
   Num. sequences **bubalus_bubalis_bubalis** : 1

> **REMOVE subspecies, less genes**   
   Num. sequences tragelaphus_scriptus : 26   
   Num. sequences **tragelaphus_scriptus_scriptus** : 2

> **REMOVE subspecies, less genes**   
   Num. sequences kobus_leche : 27   
   Num. sequences **kobus_leche_leche** : 2

> **REMOVE subspecies, less genes**   
   Num. sequences kobus_ellipsiprymnus : 26   
   Num. sequences **kobus_ellipsiprymnus_ellipsiprymnus** : 2

> **REMOVE subspecies, less genes**   
   Num. sequences gazella_gazella : 26   
   Num. sequences **gazella_gazella_gazella** : 1

> **REMOVE subspecies, less genes**   
   Num. sequences ovis_vignei : 27   
   Num. sequences **ovis_vignei_vignei** : 3

> **REMOVE subspecies, less genes**   
   Num. sequences ovis_ammon : 29   
   Num. sequences **ovis_ammon_ammon** : 2

> **REMOVE subspecies, less genes**   
   Num. sequences ovis_canadensis : 30   
   Num. sequences **ovis_canadensis_canadensis** : 4

> **REMOVE subspecies, less genes**   
   Num. sequences ovis_dalli : 15   
   Num. sequences **ovis_dalli_dalli** : 2

> **REMOVE subspecies, less genes**   
   Num. sequences rupicapra_rupicapra : 30   
   Num. sequences **rupicapra_rupicapra_rupicapra** : 1

> **REMOVE subspecies, less genes**   
   Num. sequences rupicapra_pyrenaica : 26   
   Num. sequences **rupicapra_pyrenaica_pyrenaica** : 20   
   ```
   # rupicapra_pyrenaica_pyrenaica --> "ATP6" "CO1"  "CO2"  "CO3"  "CYB"  "ND1"  "ND2"  "ND4"  "ND5"  "RNR1" "RNR2"   
   # rupicapra_pyrenaica           --> "ATP6" "ATP8" "CO1"  "CO2"  "CO3"  "CYB"  "ND1"  "ND2"  "ND3"  "ND4"  "ND4L" "ND5"  "RNR1" "RNR2"   
   ```

> **REMOVE subspecies, less genes**   
   Num. sequences pseudois_nayaur : 27   
   Num. sequences **pseudois_nayaur_nayaur** : 20   
   ```
   # pseudois_nayaur_nayaur --> "ATP6" "CO1"  "CO2"  "CO3"  "CYB"  "ND1"  "ND2"  "ND4"  "ND5"  "RNR1" "RNR2"   
   # pseudois_nayaur        --> "ATP6" "ATP8" "CO1" "CO2" "CO3" "CYB" "ND1" "ND2" "ND3" "ND4" "ND4L" "ND5" "ENSG00000125538" "RNR1" "RNR2"   
   ```

> **REMOVE subspecies, less genes**   
   Num. sequences capra_ibex : 27   
   Num. sequences **capra_ibex_ibex** : 6

All together, the taxa to remove are the following:

```
Taxa to remove:
alces_alces_alces
alces_alces_gigas
alces_alces_shirasi
axis_porcinus_porcinus
pseudonovibos_spiralis
tragelaphus_scriptus_scriptus
bos_javanicus_birmanicus
neophocaena_asiaeorientalis_asiaeorientalis
sus_barbatus_barbatus
bubalus_bubalis_bubalis
pseudois_nayaur_nayaur
gazella_gazella_gazella
```

In addition, we carried out other checks with regards to the species names and their 
placement in the phylogeny. You can find a list of these checks in the file 
`L.cetartiodactyla_taxonomy_check.xlsx` as well as below:

```
Artiodactyla	Genus	capra_pyrenaica                                                         --> C. pyrenaica should not cluster with R. pyrenaica, it should cluster with C. ibex. 
Artiodactyla	Genus	tursiops_australis,tursiops_truncatus,tursiops_aduncus                  --> (((((((S. longinostris, S. clymene), L. hoset), (D.delphis)), S. coeruleoalba), ((S.attenuata,S.frontalis), (T.aduncus, T.truncatus))) Sousa_chinensis), Sotalia_guianensis);
Artiodactyla	Genus	lama_glama_argentina                                                    --> Fix "lama_glama_argentina" within the clade where "lama_glama" is. 
Artiodactyla	Genus	mazama_americana, mazama_temama                                         --> Check: (((Od.therest, (Od. Sitkensi, M. pandora)), (M.temama, M.americana)), M. rufina);
Artiodactyla	Genus	pudu_puda, pudu_mephistophiles                                          --> Check: Move pudu_mephistophiles and place it (a) outgroup of pudu_puda or (b) sister taxa with pudu_puda
Artiodactyla	Genus	lagenorhynchus_acutus, lagenorhynchus_albirostris                       --> Check: The only difference is the species Orcinus orca, which actually seems to either be "((((rest), L. albirostris), O. orca), L. acutus);" or "((((REST), O. orca), L. albirostris), L. acutus);".
Artiodactyla	Genus	cervus_albirostris                                                      --> It should cluster with rest of Cervus. Check: Possible placements would be: "(((canadensis, nibbus),albirostris), elaphus);" or "(((canadensis,nibbus),elaphus),albirostris);"
Artiodactyla	Genus	nilgiritragus_hylocrius                                                 --> Check: Try placing it as outgroup of the rest of ovis, not cluster together with ovis.
Artiodactyla	Genus	hemitragus_jemlahicus                                                   --> Place as sister taxa with capra_sibirica
Artiodactyla	Genus	hylochoerus_meinertzhageni                                              --> Move and place it as an outgroup as we need the two species of Phacochoerus to cluster together
Artiodactyla	Genus	indopacetus_pacificus                                                   --> Place it as an outgroup of all Mesoplodon. I.e., (((Mesoplodon.sp), indopacetus_pacificus), Hyperoodon);
Artiodactyla	Genus	phocoenoides_dalli                                                      --> Check: Place as a sister taxon of Phocoena phocoena.
Artiodactyla	Species	tragelaphus_scriptus                                                    --> Check: Sister taxon of T. angasii.
Artiodactyla	Species	gazella_subgutturosa, gazella_marica                                    --> Check: Gazella subgutturosa should be in same clade as Gazella marica.
Artiodactyla	Species	bos_javanicus, bos_javanicus_lowi                                       --> Check: bos_javanicus should be sister taxon with gaurus.
Artiodactyla	Species	neophocaena_phocaenoides                                                --> Check: subspecies should cluster with this species.                                 
Artiodactyla	Species	cervus_timorensis                                                       --> Check: cervus_timorensis should cluster with cervus_timorensis_macassaricus
Artiodactyla	Species	cervus_canadensis_songaricus                                            --> Check: maybe force monophyletic with rest of C. canadensis.
Artiodactyla	Species	moschus_berezovskii_caobangis                                           --> Check: should cluster with M. berezovskii species.
Artiodactyla	Species	rupicapra_pyrenaica_ornata                                              --> Check: should cluster with R. pyrenaica.
Artiodactyla	Species	capricornis_milneedwardsii_maritimus                                    --> Check: should cluster with species capricornis_milneedwardsis.
Artiodactyla	Species	pseudois_schaeferi                                                      --> Check: pseudois_schaeferi should be outgroup of rest of pseudois_nayaur species.
Artiodactyla	Species	ovis_orientalis                                                         --> Check: All O. orientalis should be clustering in same clade.
Artiodactyla	Species	ovis_vignei_blanfordi, ovis_vignei_kermanensis, ovis_vignei_cycloceros  --> Check: Try to force all Ovis vignei together.
Artiodactyla	Species	rupicapra_rupicapra_cartusiana                                          --> Check: Should be clustering with R. rupicapra.
Artiodactyla	Species	All 6 species of Ovis orientalis and ovis_ammon_collium.                --> They are clustered in different clades. Check them all.
Artiodactyla	Species	capra_aegagrus_blythi                                                   --> Check: Force to be in smae clade with C. aegagrus sp.
Artiodactyla	Species	odocoileus_hemionus_sitkensis                                           --> Check:  (((Od.therest, (Od. Sitkensi, M. pandora)), (M.temama, M.americana)), M. rufina);
Artiodactyla	Species	lama_glama_argentina,lama_glama_chaku                                   --> Check: Bring lama_glama_chaku with rest of lama_glama sp.
```

Then, we opened the `RAxML_bestTree.BS_ML_GTRCAT` with `FigTree` to individually check the taxa with dubious 
plaecment. Those that are to be deleted are coloured in red, while those that they cluster
within the same clade with the species they belong to but the placement might be dubious have
been coloured in orange. The names for the latter will be changed so they have a `~` at the end. 
This means that, even though this is the placement found in the ML tree, other research have
found they cluster elsewhere. Therefore, we flagged them for future and further analyses that might be carried out,
but that are not part of the scope of this project. 
The output files (`png` and `FigTree` files) have been saved in the
[`checked_aln/taxonomical_check`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/laurasiatheria_cetartiodactyla/filter_aln/checked_aln/taxonomical_check)
directory. 

The taxa to be renamed and/or removed are the following:

```
Taxa to rename:
bison_bison,bos_bison
bison_bonasus,bos_bonasus
sus_salvanius,porcula_salvania
rusa_alfredi,cervus_alfredi
przewalskium_albirostris,cervus_albirostris
hemitragus_hylocrius,nilgiritragus_hylocrius
hemitragus_jayakari,arabitragus_jayakari
taurotragus_derbianus,tragelaphus_derbianus
gazella_subgutturosa_marica,gazella_marica
rusa_timorensis,cervus_timorensis
rusa_timorensis_macassaricus,cervus_timorensis_macassaricus
rusa_unicolor,cervus_unicolor
rusa_unicolor_cambojensis,cervus_unicolor_cambojensis
rusa_unicolor_swinhoei,cervus_unicolor_swinhoei
tursiops_australis,tursiops_australis~
tursiops_truncatus,tursiops_truncatus~
tursiops_aduncus,tursiops_aduncus~
lama_glama_argentina,lama_glama_argentina~
mazama_americana,mazama_americana~
mazama_temama,mazama_temama~
pudu_puda,pudu_puda~
pudu_mephistophiles,pudu_mephistophiles~
lagenorhynchus_acutus,lagenorhynchus_acutus~
lagenorhynchus_albirostris,lagenorhynchus_albirostris~
cervus_albirostris,cervus_albirostris~
hemitragus_jemlahicus,hemitragus_jemlahicus~
hylochoerus_meinertzhageni,hylochoerus_meinertzhageni~
indopacetus_pacificus,indopacetus_pacificus~
phocoenoides_dalli,phocoenoides_dalli~
tragelaphus_scriptus,tragelaphus_scriptus~
gazella_subgutturosa,gazella_subgutturosa~
gazella_marica,gazella_marica~
bos_javanicus,bos_javanicus~
neophocaena_phocaenoides,neophocaena_phocaenoides~
cervus_timorensis,cervus_timorensis~
cervus_canadensis_songaricus,cervus_canadensis_songaricus~
moschus_berezovskii_caobangis,moschus_berezovskii_caobangis~
rupicapra_pyrenaica_ornata,rupicapra_pyrenaica_ornata~
capricornis_milneedwardsii_maritimus,capricornis_milneedwardsii_maritimus~
pseudois_schaeferi,pseudois_schaeferi~
ovis_orientalis,ovis_orientalis~
ovis_vignei_blanfordi,ovis_vignei_blanfordi~
ovis_vignei_kermanensis,ovis_vignei_kermanensis~
ovis_vignei_cycloceros,ovis_vignei_cycloceros~
rupicapra_rupicapra_cartusiana,rupicapra_rupicapra_cartusiana~
ovis_orientalis_anatolica,ovis_orientalis_anatolica~
ovis_orientalis_isphahanica,ovis_orientalis_isphahanica~
ovis_ammon_collium,ovis_ammon_collium~
capra_aegagrus_blythi,capra_aegagrus_blythi~
odocoileus_hemionus_sitkensis,odocoileus_hemionus_sitkensis~
lama_glama_chaku,lama_glama_chaku~
nilgiritragus_hylocrius,nilgiritragus_hylocrius~

Taxa to remove:
capra_pyrenaica
```

**NOTE:**Data subsets (both alignments, 12CP and 3CP) for Afrotheria, Xenarthra, Euarchonta, and Marsupialia had already undergone 
these checks before 2018 (i.e., they had already been "cleaned", while the data subset for Lagomorpha 
had not yet).
Therefore, you do not see the last part of the filtering described in this section 
in the corresponding `README.md` files for these data subsets or the directory
[`taxonomical_check`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/laurasiatheria_cetartiodactyla/filter_aln/checked_aln/taxonomical_check)
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

```sh
# Run the next code from `laurasiatheria_cetartiodactyla/filter_aln/checked_aln/unfiltered_aln`
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
   in RStudio and change line 25 so it is `subt    <- "laurasiatheria_cetartiodactyla"`, uncomment line 27,
   and comment line 29. Now, we can run it from RStudio.
   This script will generate a concatenated alignment file with all partitions, as well as 
   one alignment file for each individual partitions, inside a new dir called `00_mammal_alns/laurasiatheria_cetartiodactyla/unfiltered/`
   inside [`00_Data_filtering/01_alignments/`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/01_alignments). 
   Log and RData files will be saved inside
   [`00_Data_filtering/01_alignments/Rout`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/01_alignments/Rout/log_concatenation).   
      > NOTE 1: Updated `partitions.txt` file inside `00_mammal_alns/laurasiatheria_cetartiodactyla/unfiltered/` generated.   
	  > NOTE 2: Output file called `check_appends.txt` to check that 3nt partition has been appended
	    to the right line of the alignment can be found in `Rout/Rdata". 
		If "0" in second column, sth wrong! So far, everything seems OK!   
	  > NOTE 3: In case function `load_files()` in line 180 does not work because of the path,
	  just uncomment lines 184-188 and it should work.   
   * Paths have been automatically set according to current file architecture in `mammals` dir, 
   do not change paths in the Rscript! 
   
**NOTE: You can see that we have used a different R script in this step because, as mentioned above,**
**the data for laurasiatheria cetartiodactyla were not at the same stage as data subsets for Afrotheria, Xenarthra,**
**Euarchonta and Marsupialia.**

Now, we are going to use this unfiltered alignment phylip to apply the next filtering steps. 

## 3. Apply the filtering step
Run the following commands from `laurasiatheria_cetartiodactyla/filter_aln/checked_aln`

```sh
# Run the next code from `laurasiatheria_cetartiodactyla/filter_aln/checked_aln`

# 1. Fix species names in `names.txt`. We used the RAxML file, the `names.txt` file, and 
#    the `lineages.txt` that had been generated when filtering the alignment 12CP.
#    The ML tree file is now saved inside `filter_aln/checked_aln`, from where the 
#    python script will be run. The `names.txt` and the `lineage.txt` files are in
#    the `filter_aln` dir, one up. The alignment is the concatenated one that was generated in the 
#    previous step and can be found inside
#   `01_alignments/00_mammal_alns/laurasiatheria_cetartiodactyla/unfiltered/laurasiatheria_cetartiodactyla.aln`.
#
#    Note that we just copy this file inside `checked_aln` to proceed with the filtering. We also generate 
#    a copy of the original tree file 
cp ../../../../01_alignments/00_mammal_alns/laurasiatheria_cetartiodactyla/unfiltered/laurasiatheria_cetartiodactyla.aln .
mkdir original_tree 
cp RAxML_bestTree.BS_ML_GTRCAT original_tree/RAxML_bestTree.BS_ML_GTRCAT
cat > _species_name_change.txt << EOF
bison_bison,bos_bison
bison_bonasus,bos_bonasus
sus_salvanius,porcula_salvania
rusa_alfredi,cervus_alfredi
przewalskium_albirostris,cervus_albirostris~
hemitragus_hylocrius,nilgiritragus_hylocrius~
hemitragus_jayakari,arabitragus_jayakari
taurotragus_derbianus,tragelaphus_derbianus
gazella_subgutturosa_marica,gazella_marica~
rusa_timorensis,cervus_timorensis~
rusa_timorensis_macassaricus,cervus_timorensis_macassaricus
rusa_unicolor,cervus_unicolor
rusa_unicolor_cambojensis,cervus_unicolor_cambojensis
rusa_unicolor_swinhoei,cervus_unicolor_swinhoei
tursiops_australis,tursiops_australis~
tursiops_truncatus,tursiops_truncatus~
tursiops_aduncus,tursiops_aduncus~
lama_glama_argentina,lama_glama_argentina~
mazama_americana,mazama_americana~
mazama_temama,mazama_temama~
pudu_puda,pudu_puda~
pudu_mephistophiles,pudu_mephistophiles~
lagenorhynchus_acutus,lagenorhynchus_acutus~
lagenorhynchus_albirostris,lagenorhynchus_albirostris~
hemitragus_jemlahicus,hemitragus_jemlahicus~
hylochoerus_meinertzhageni,hylochoerus_meinertzhageni~
indopacetus_pacificus,indopacetus_pacificus~
phocoenoides_dalli,phocoenoides_dalli~
tragelaphus_scriptus,tragelaphus_scriptus~
gazella_subgutturosa,gazella_subgutturosa~
bos_javanicus,bos_javanicus~
neophocaena_phocaenoides,neophocaena_phocaenoides~
cervus_canadensis_songaricus,cervus_canadensis_songaricus~
moschus_berezovskii_caobangis,moschus_berezovskii_caobangis~
rupicapra_pyrenaica_ornata,rupicapra_pyrenaica_ornata~
capricornis_milneedwardsii_maritimus,capricornis_milneedwardsii_maritimus~
pseudois_schaeferi,pseudois_schaeferi~
ovis_orientalis,ovis_orientalis~
ovis_vignei_blanfordi,ovis_vignei_blanfordi~
ovis_vignei_kermanensis,ovis_vignei_kermanensis~
ovis_vignei_cycloceros,ovis_vignei_cycloceros~
rupicapra_rupicapra_cartusiana,rupicapra_rupicapra_cartusiana~
ovis_orientalis_anatolica,ovis_orientalis_anatolica~
ovis_orientalis_isphahanica,ovis_orientalis_isphahanica~
ovis_ammon_collium,ovis_ammon_collium~
capra_aegagrus_blythi,capra_aegagrus_blythi~
odocoileus_hemionus_sitkensis,odocoileus_hemionus_sitkensis~
lama_glama_chaku,lama_glama_chaku~
EOF

mkdir screen_logs
python ../../../../../../src/fix_species_names.py _species_name_change.txt RAxML_bestTree.BS_ML_GTRCAT laurasiatheria_cetartiodactyla.aln ../names.txt ../lineage.txt > screen_logs/01_log1_fixspecies.txt

##>> NOTE: I needed to update the python script as it was not working. I just modifed `fix_alignment()` 
##>>       so the output file was different from the one being read

## Move unncessary output files to `00_filt1/` and rename output alignment
mkdir 00_filt1
mv _species_name_change.txt lineage.txt lineage.txt.1.bak 00_filt1/
rm laurasiatheria_cetartiodactyla.aln
mv laurasiatheria_cetartiodactyla_out.phylip alignment.phylip

# 2. Remove species 

# 2.1. Create input file for python script `prune_alignment.py`. This text file 
#    contains a list of the species that need to be removed 
cat > _species_remove.txt << EOF
alces_alces_alces
alces_alces_gigas
alces_alces_shirasi
axis_porcinus_porcinus
pseudonovibos_spiralis
tragelaphus_scriptus_scriptus
bos_javanicus_birmanicus
neophocaena_asiaeorientalis_asiaeorientalis
sus_barbatus_barbatus
bubalus_bubalis_bubalis
pseudois_nayaur_nayaur
gazella_gazella_gazella
rangifer_tarandus_tarandus
odocoileus_hemionus_hemionus
rucervus_eldii_eldii
cervus_elaphus_elaphus
lama_guanicoe_guanicoe
sus_scrofa_scrofa
inia_geoffrensis_geoffrensis
cephalorhynchus_hectori
kobus_leche_leche
kobus_ellipsiprymnus_ellipsiprymnus
ovis_vignei_vignei
ovis_ammon_ammon
ovis_canadensis_canadensis
ovis_dalli_dalli
rupicapra_rupicapra_rupicapra
rupicapra_pyrenaica_pyrenaica
capra_ibex_ibex
capra_pyrenaica
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
   in RStudio and change line 24 so it is `subt    <- "laurasiatheria_cetartiodactyla"` and uncomment line 26. 
   Now, we can run it from RStudio. This script will generate a concatenated alignment file with all partitions, as well as 
   one alignment file for each individual partitions,
   inside [`01_alignments`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/01_alignments).
   Log files and Rdata can be found [here](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/01_alignments/Rout/log_concatenation). 
      > NOTE: Updated `partitions.txt` file inside `00_mammal_alns/laurasiatheria_cetartiodactyla/` generated.   
   * Paths have been automatically set according to current file architecture in `mammals` dir, 
   do not change paths in the Rscript! 

The final alignments generated at the end of this step can be downloaded from 
[here](https://www.dropbox.com/s/94z6yura9rtaobb/SeqBayesS2_Raln_artiodactyla.zip?dl=0).
They should be saved here if the same file architecture as the one set in the R scripts 
is to be used: `00_Data_filtering/01_alignments/00_mammals_alns/laurasiatheria_cetartiodactyla`.