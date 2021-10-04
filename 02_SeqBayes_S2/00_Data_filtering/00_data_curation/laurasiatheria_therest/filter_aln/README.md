# LAURASIATHERIA THE REST - Filtering alignment
The initial files that you can find in this directory are the following:

```
laurasiatheria_therest
   |- filter_aln
         |- checked_aln                        
         |     |- taxonomical_check
         |           |- laurasiatheria_therest.png
         |           |- laurasiatheria_therest_FigTree
         |           |- RAxML_bestTree.BS_ML_GTRCAT
         |     
         |- extra_filtering   <-- You will read more about it at the end of this guideline                     
         |     
         |- L.threst_taxonomy_check.xlsx
         |- lineage.txt 
         |- names.txt       
         |- parse_lineage.R 
         |- partitions.txt 
         |- README.md
         |- summary.html 
```
		 
The output files during the filtering step that is described below as well as the input 
files needed have been zipped in a file as they are very large. You can download 
this file [here](https://www.dropbox.com/s/n0r1pq16ygu7ddv/SeqBayesS2_filtaln_laurasiatheria_therest.zip?dl=0).
To start the filtering step, you should have the following files arranged in the file 
architecture detailed above (you can obtain the files once you unzip the file
provided in the link above): 

```
laurasiatheria_therest 
   |- filter_aln
         |- checked_aln                        
         |     |- unfiltered_aln
         |     |    |- alignment_nt3cp.phylip   # Alignment with 3CP 
         |     |- RAxML_bestTree.BS_ML_GTRCAT   # Tree you can find in `original_tree` directory
         |     |- taxonomical_check                   # Visual checks to evaluate 
         |          |- laurasiatheria_therest.png     # dubious taxa placement
         |          |- laurasiatheria_therest_FigTree
         |          |- RAxML_bestTree.BS_ML_GTRCAT    
         |     
         |- extra_filtering   <-- You will read more about it at the end of this guideline                     
         |        
         |- alignment.phylip # Alignment with 12CP
         |- L.threst_taxonomy_check.xlsx
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
Within the [`filter_aln`](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/laurasiatheria_therest/filter_aln)
directory, you will find different files as detailed above. 
Specifically, the R script [`parse_lineage.R`](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/laurasiatheria_therest/filter_aln/parse_lineage.R)
was written to carry out a first taxonomic filtering. Note that this 
R script will run if you have the same file architecture in this GitHub repository (i.e., it uses 
a function within the R script [`Filter_lineages.R`](https://github.com/sabifo4/mammals_dating/blob/main/src/Filter_lineages.R)
in the `src` directory and the [`genes.txt`](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/genes.txt)
file. The messages printed out by this script are the following:   

```
GENUS   Conepatus  has  2  taxa. Species  conepatus_chinga  has  5  genes that are not shared
by any of the other taxa. It should be deleted!
GENUS   Conepatus  has  2  taxa. Species  conepatus_mesoleucus  has  3  genes that are not share
by any of the other taxa. It should be deleted!```
>> ACTION: conepatus_chinga and conepatus_mesoleucus were mereged. Label is "conepatus".

```
GENUS   Spilogale  has  2  taxa. Species  spilogale_gracilis_latifrons  has  1  genes that are not
shared by any of the other taxa. It should be deleted!
GENUS   Spilogale  has  2  taxa. Species  spilogale_putorius  has  19  genes that are not shared
by any of the other taxa. It should be deleted!```
>> ACTION: spilogale_putorius and spilogale_gracilis_latifrons were merged. Label is "spilogale".

```
GENUS   Mydaus  has  2  taxa. Species  mydaus_javanensis  has  6  genes that are not shared by any
of the other taxa. It should be deleted!
GENUS   Mydaus  has  2  taxa. Species  mydaus_marchei  has  3  genes that are not shared by any
of the other taxa. It should be deleted!```
>> ACTION: mydaus_javanensis and mydaus_marchei were merged. Label is "mydaus".

```
GENUS   Nyctereutes  has  3  taxa. Species  nyctereutes_procyonoides_procyonoides  has  2  genes
that are not shared by any of the other taxa. It should be deleted!```
>> ACTION: nyctereutes_procyonoides_procyonoides was removed.   

```
GENUS   Ursus  has  11  taxa. Species  ursus_arctos_isabellinus  has  1  genes that are not shared
by any of the other taxa. It should be deleted!```
>> ACTION: ursus_arctos_isabellinus was removed.

**NOTE**
Log files will be found in this same directory with the name `log_taxaNOTIN_<level>.txt`,
being `level` the one that has been checked (i.e., family, genus, order, or subfamily). In addition, 
you can also generate the `levels.checked.RData`, which you can use to explore the taxonomical 
levels. 

# 2. First checks before applying filtering
First, we checked if there were any further species that should be removed:

```sh
# Run from `checked_aln` directory 
grep -o '[a-z]*_[a-z]*_[a-z]*' RAxML_bestTree.BS_ML_GTRCAT | sed 's/\_/\t/g' | awk '$2==$3{print $1,$2,$3}' | awk '$1!="ft"{print $1,$2,$3}' > subsp_check.txt

# Count 
num=$( grep -o '[a-z]*_[a-z]*_[a-z]*' RAxML_bestTree.BS_ML_GTRCAT | sed 's/\_/\t/g' | awk '$2==$3{print $1,$2,$3}' | awk '$1!="pdh"{print $1,$2,$3}' | wc -l)
printf "There are $num subspecies with species\n"
## There are 29 subspecies with species
```
The output names have been saved in a file called `subsp_check.txt` to further explore it:

```sh
# Run from `checked_aln` directory 
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

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences felis_silvestris : 31   
   Num. sequences **felis_silvestris_silvestris** : 3

> **REMOVE SPECIES, we give priority to keep nuclear genes**  
   Num. sequences **urocyon_littoralis** : 5   
   Num. sequences urocyon_littoralis_littoralis : 20   
   ```
   # urocyon_littoralis_littoralis --> "ATP6" "CO1"  "CO2"  "CO3"  "CYB"  "ND1"  "ND2"  "ND4" "ND5"  "RNR1" "RNR2"
   # urocyon_littoralis            --> "ENSG00000012048" "ENSG00000112964" "ENSG00000166349" "ENSG00000176273" "ENSG00000176697"
   ```
   
> **REMOVE SUBSPECIES, less genes**   
   Num. sequences canis_lupus : 38   
   Num. sequences **canis_lupus_lupus** : 26

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences bassaricyon_gabbii : 12   
   Num. sequences **bassaricyon_gabbii_gabbii** : 2

> NO ACTION   
   Num. sequences bassaricyon_medius : 0   
   Num. sequences bassaricyon_medius_medius : 2

> NO ACTION   
   Num. sequences bassaricyon_neblina : 0   
   Num. sequences bassaricyon_neblina_neblina : 2

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences martes_melampus : 34   
   Num. sequences **martes_melampus_melampus** : 4

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences meles_meles : 43   
   Num. sequences **meles_meles_meles** : 3

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences ursus_thibetanus : 34   
   Num. sequences **ursus_thibetanus_thibetanus** : 26

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences odobenus_rosmarus : 39   
   Num. sequences **odobenus_rosmarus_rosmarus** : 27

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences pusa_hispida : 37   
   Num. sequences **pusa_hispida_hispida** : 6

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences phoca_vitulina : 46   
   Num. sequences **phoca_vitulina_vitulina** : 2

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences dicerorhinus_sumatrensis : 29   
   Num. sequences **dicerorhinus_sumatrensis_sumatrensis** : 1

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences rhinoceros_sondaicus : 26   
   Num. sequences **rhinoceros_sondaicus_sondaicus** : 1

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences ceratotherium_simum : 57   
   Num. sequences **ceratotherium_simum_simum** : 5

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences equus_hemionus : 26   
   Num. sequences **equus_hemionus_hemionus** : 2

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences crocidura_suaveolens : 3   
   Num. sequences **crocidura_suaveolens_suaveolens** : 2
   ```
   crocidura_suaveolens_suaveolens --> "RNR2" "ENSG00000012048"
   crocidura_suaveolens            --> "CYB"  "ENSG00000012048"
   ```

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences suncus_murinus : 31   
   Num. sequences **suncus_murinus_murinus** : 3

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences neomys_anomalus : 4   
   Num. sequences **neomys_anomalus_anomalus** : 2
   ```
   # neomys_anomalus_anomalus --> "CYB"
   # neomys_anomalus          --> "CYB" "RNR2" "ENSG00000012048"
   ```
   
> NO ACTION   
   Num. sequences sorex_dispar : 0   
   Num. sequences sorex_dispar_dispar : 2

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences sorex_hoyi : 2   
   Num. sequences **sorex_hoyi_hoyi** : 2   
   ```
   # sorex_hoyi_hoyi --> "CYB"
   # sorex_hoyi      --> "CYB"
   ```
   
> **REMOVE SUBSPECIES, less genes**   
   Num. sequences sorex_thompsoni : 2   
   Num. sequences **sorex_thompsoni_thompsoni** : 2
   ```
   # sorex_thompsoni_thompsoni --> "CYB"
   # sorex_thompsoni           --> "CYB"
   ```

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences sorex_palustris : 9   
   Num. sequences **sorex_palustris_palustris** : 2

> No action   
   Num. sequences sorex_bendirii : 0   
   Num. sequences sorex_bendirii_bendirii : 2

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences neofelis_nebulosa : 32   
   Num. sequences **neofelis_nebulosa_nebulosa** : 1

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences panthera_leo : 39   
   Num. sequences **panthera_leo_leo** : 20

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences panthera_tigris : 34   
   Num. sequences **panthera_tigris_tigris** : 20

All together, the taxa to remove are the following:

```
Taxa to remove:
sorex_palustris_palustris
meles_meles_meles
equus_hemionus_hemionus
crocidura_suaveolens_suaveolens
suncus_murinus_murinus
phoca_vitulina_vitulina
```

In addition, we carried out other checks with regards to the species names and their 
placement in the phylogeny. You can find a list of these checks in the file 
`L.therest_taxonomy_check.xlsx` as well as below:

```
Genus   Mustela               CHECK      neovison_vison                       --> It clusters with the rest of Mustela species as expected. Same clustering shown in Law et al. 2018, with the differnce that mustela_africana and mustela_felipei are sister taxa and then mustela_frenata comes. We can keep our topology as it is or try this one.
Genus   Helogale              CHECK      herpestes_naso                       --> In reference it also clusters with Atilax paludinosus. Nevertheless, it is clustered as an outgroup of the clade of Herspestes together with Gallerea,Bdeogale,Rhunchogale,Paracynictis,Ichenumia species. Try this topology too if needed. The two species of Helogale follow the clustering in the topology in reference.
Genus   Martes                CHECK      eira_barbara, pekania_pennanti       --> Sato et al. 2012 suggest elevating the subgenus Pekania for Pekania to genus to accommodate the fisher (Martes pennanti) as it clusters outside the rest of Martes and so to confirm Martes to be monophyletic. Match relationships as in our topology with Martes pennanti being renamed as Pekania pennanti.
                                                                                  Maybe try to separate Eira baraba and Pekania pennanti as in Koepfli et al. 2008, but the placement of these two taxa is a bit debated as Law et al. 2018 found the same clustering as we do.
Genus   Episoriculus          CHECK      soriculus_nigrescens,                --> There is not a consensus with the palcement of soriculus_nigrescens and episoriculus_fumidus. Nevertheless, most of them place episoriculus_fumidus as the "outgroup" of the rest of taxa shown here and then soriculus_nigrescens. We can leave it as it is now or try this approach.
                                         episoriculus_fumidus
Genus   Arctocephalus         YES        arctocephalus_tropicalis             --> There is only one gene, thus this might have affected the clustering. Force monophyletic.
Genus   Galictis              YES        lyncodon_patagonicus                 --> They cluster together in a paraphyletic group but lyncodon_patagonicus is not sister taxa with one of the galictis species. Fix galictis species to be sister taxa.
Genus   Canis                 YES        lupulella_adusta, canis_lupaster,    --> Change names to lupulella_adusta, lupulella_mesomelas, lupulella_mesomelas_elongae, canis_lupaster, canis_lupus_lycaon and topology needs to be updated to the following: 
                                         lupulella_mesomelas,                     ((((Lupulella adusta, Lupulella mesomelas), Lycaon pictus), Cuon alpinus), rest_canis); || (((Lycaon pictus, (Lupulella adusta, Lupulella mesomelas), Cuon alpinus), rest_canis);
                                         lupulella_mesomelas_elongae,                  a) ((lycaon_pictus, (lupulella_adusta, (lupullela_mesomelas, lupullela_mesomelas_elongae))), (cuon_alpinus, cuon_alpinus_lepturus));. Then cluster canis_rufus with canis_latrans
                                         canis_lupus_lycaon                            b) (((lupulella_adusta, (lupullela_mesomelas, lupullela_mesomelas_elongae)), lycaon_pictus), (cuon_alpinus, cuon_alpinus_lepturus));. Then cluster canis_rufus with canis_latrans                                                    
Species Hylomys suillus       CHECK      hylomys_parvus                       --> It could be forced monophyletic, but the same rearrangement has been observed in reference, so maybe we can leave it as it is.
Species Crocidura suaveolens  CHECK      crocidura_sibirica,                  --> Synonyms: C. caspica == C. russula caspica == C. suaveolens caspica | C. monacha == C. russula monacha | C. aleksandrisi == C. suaveolens aleksandrisi. Clustering with C. aleksandrisi and C. mimula appears in REF. No information about the others, so we can keep it like that or just put C. sibirica as the outgroup of them.
                                         crocidura_caspica
Species Equus burchellii      YES        equus_burchellii_cuninghamei         --> Force monophyletic and cluster these two sp. With the rest of equus burchellii as sister taxa of equus_grevyi.
                                         equus_burchellii_antiquorum
Species Felis silvestris      YES        felis_catus                          --> Force monophyletic. Cluster all felis_sylvestris sp together and have outgroup felis_catus
Species Equus hemionus        YES        equus_kiang                          --> Force Equus hemionus sp to cluster together in the clade where E. kiang is the outgroup
Species Crocidura orientalis  YES        crocidura_orientalis_lawuana         --> Change "(((C.beccarii, C.lepidura), (C.orientalis, C.brunnea)), C. orientalis_lawuana);" into "(((C. orientalis, C.orientalis.lawuana), C. brunnes), (C.becarii, C.lepidura));"
Species Panthera tigris       YES        panthera_pardus_saxicolor            --> Force monophyletic. Put panthera_pardus_saxicolor with the rest of panthera_pardus sp.
Species Suncus murinus        YES        suncus_montanus                      --> Force monophyletic suncus_murinus sp and then get suncus_montanus as their outgroup
Species Phoca vitulina        YES        phoca_largha                         --> Force monophyletic with phoca_largha being the outgroup for the rest of phoca_vitulina sp.
Species Sylvisorex granti     YES        sylvisorex_vulcanorum                --> Force monophyletic sylvisorex_granti with its subspecies and then put sylvisorex_vulcanorum.
Species Mogera insularis      YES        mogera_insularis_latouchei           --> Force monophyletic mogera_insularis with its subspecies mogera_insularis_latouchei and then put mogera_kanoana
```

Then, we opened the `RAxML_bestTree.BS_ML_GTRCAT` with `FigTree` to individually check the taxa with dubious 
plaecment. Those that are to be deleted are coloured in red, while those that they cluster
within the same clade with the species they belong to but the placement might be dubious have
been coloured in orange. The names for the latter will be changed so they have a `~` at the end. 
This means that, even though this is the placement found in the ML tree, other research have
found they cluster elsewhere. Therefore, we flagged them for future and further analyses that might be carried out,
but that are not part of the scope of this project. 
The output files (`png` and `FigTree` files) have been saved in the
[`checked_aln/taxonomical_check`](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/laurasiatheria_therest/filter_aln/checked_aln/taxonomical_check)
directory. 

The taxa to be renamed and/or removed are the following:

```
Taxa to rename:
martes_pennanti,pekania_pennanti
canis_adustus,lupulella_adusta
canis_mesomelas,lupulella_mesomelas
canis_mesomelas_elongae,lupulella_mesomelas_elongae
canis_lupus_lupaster,canis_lupaster
canis_lycaon,canis_lupus_lycaon
neovison_vison,neovison_vison~
herpestes_naso,herpestes_naso~
herpestes_ichneumon,herpestes_ichneumon~
eira_barbara,eira_barbara~
pekania_pennanti,pekania_pennanti~
soriculus_nigrescens,soriculus_nigrescens~
episoriculus_fumidus,episoriculus_fumidus~
arctocephalus_tropicalis,arctocephalus_tropicalis~
lyncodon_patagonicus,lyncodon_patagonicus~
lupulella_adusta,lupulella_adusta~
canis_lupaster,canis_lupaster~
lupulella_mesomelas,lupulella_mesomelas~
lupulella_mesomelas_elongae,lupulella_mesomelas_elongae~
canis_lupus_lycaon,canis_lupus_lycaon~
hylomys_parvus,hylomys_parvus~
crocidura_sibirica,crocidura_sibirica~
crocidura_caspica,crocidura_caspica~
equus_burchellii_cuninghamei,equus_burchellii_cuninghamei~
equus_burchellii_antiquorum,equus_burchellii_antiquorum~
felis_catus,felis_catus~
equus_kiang,equus_kiang~
crocidura_orientalis_lawuana,crocidura_orientalis_lawuana~
panthera_pardus_saxicolor,panthera_pardus_saxicolor~
suncus_montanus,suncus_montanus~
phoca_largha,phoca_largha~
sylvisorex_vulcanorum,sylvisorex_vulcanorum~
mogera_insularis_latouchei,mogera_insularis_latouchei~
```
   
   
**NOTE:** Data subsets (both alignments, 12CP and 3CP) for Afrotheria, Xenarthra, Euarchonta, and Marsupialia had already undergone 
these checks before 2018 (i.e., they had already been "cleaned", while the data subset for Lagomorpha 
had not yet).
Therefore, you do not see the last part of the filtering described in this section 
in the corresponding `README.md` files for these data subsets or the directory
[`taxonomical_check`](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/laurasiatheria_therest/filter_aln/checked_aln/taxonomical_check)
You will see that the `filter_aln` directory for these four data subsets
contains a csv file with the details that were followed to filter the corresponding alignments (12CP and 3CP). 
Note that this csv file is equivalent to the excel sheet you find in this directory for 
data subset Laurasiatheria the rest.

## 2. Check alignment with 3CP partition
Before we concatenate the alignment with the first and second codon positions (12CP), file `alignment.phylip`
(you will find it once you unzip the file which link is provided above),
and the alignment with the third codon
positions (3CP), file `alignment_nt3cp.phylip` (you will find it once 
you unzip the file which link provided above, inside a directory called `checked_aln` directory),
we ran the following code to make sure they had both undergone the same filtering steps and that were 
at the same "filtering stage":

```sh
# Run the next code from `laurasiatheria_therest/filter_aln/checked_aln/unfiltered_aln`
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
   in RStudio and change line 25 so it is `subt    <- "laurasiatheria_therest"` and uncomment line 27. Now, we
   can run it from RStudio. This script will generate a concatenated alignment file with all
   partitions, as well as one alignment file for each individual partitions, inside a new dir called `00_mammal_alns/laurasiatheria_therest/unfiltered/`
   inside [`00_Data_filtering/01_alignments/`](/02_SeqBayes_S2/00_Data_filtering/01_alignments). 
   Log and RData files will be saved inside
   [`00_Data_filtering/01_alignments/Rout`](/02_SeqBayes_S2/00_Data_filtering/01_alignments/Rout/log_concatenation).   
      > NOTE 1: Updated `partitions.txt` file inside `00_mammal_alns/laurasiatheria_therest/unfiltered/` generated.   
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
Run the following commands from `laurasiatheria_therest/filter_aln/checked_aln`

```sh
# Run the next code from `laurasiatheria_therest/filter_aln/checked_aln`

# 1. Fix species names in `names.txt`. We used the RAxML file, the `names.txt` file, and 
#    the `lineages.txt` that had been generated when filtering the alignment 12CP.
#    The ML tree file is now saved inside `filter_aln/checked_aln`, from where the 
#    python script will be run. The `names.txt` and the `lineage.txt` files are in
#    the `filter_aln` dir, one up. The alignment is the concatenated one that was generated in the 
#    previous step and can be found inside
#    `01_alignments/00_mammal_alns/laurasiatheria_therest/unfiltered/laurasiatheria_therest.aln`.
#    Note that we just copy this file inside `checked_aln` to proceed with the filtering. We also generate 
#    a copy of the original tree file:
cp ../../../../01_alignments/00_mammal_alns/laurasiatheria_therest/unfiltered/laurasiatheria_therest.aln .
mkdir original_tree 
cp RAxML_bestTree.BS_ML_GTRCAT original_tree/RAxML_bestTree.BS_ML_GTRCAT
cat > _species_name_change.txt << EOF
martes_pennanti,pekania_pennanti~
canis_adustus,lupulella_adusta~
canis_mesomelas,lupulella_mesomelas~
canis_mesomelas_elongae,lupulella_mesomelas_elongae~
canis_lupus_lupaster,canis_lupaster~
canis_lycaon,canis_lupus_lycaon~
neovison_vison,neovison_vison~
herpestes_naso,herpestes_naso~
herpestes_ichneumon,herpestes_ichneumon~
eira_barbara,eira_barbara~
soriculus_nigrescens,soriculus_nigrescens~
episoriculus_fumidus,episoriculus_fumidus~
arctocephalus_tropicalis,arctocephalus_tropicalis~
lyncodon_patagonicus,lyncodon_patagonicus~
hylomys_parvus,hylomys_parvus~
crocidura_sibirica,crocidura_sibirica~
crocidura_caspica,crocidura_caspica~
equus_burchellii_cuninghamei,equus_burchellii_cuninghamei~
equus_burchellii_antiquorum,equus_burchellii_antiquorum~
felis_catus,felis_catus~
equus_kiang,equus_kiang~
crocidura_orientalis_lawuana,crocidura_orientalis_lawuana~
panthera_pardus_saxicolor,panthera_pardus_saxicolor~
suncus_montanus,suncus_montanus~
phoca_largha,phoca_largha~
sylvisorex_vulcanorum,sylvisorex_vulcanorum~
mogera_insularis_latouchei,mogera_insularis_latouchei~
EOF

mkdir screen_logs
python ../../../../../../src/fix_species_names.py _species_name_change.txt RAxML_bestTree.BS_ML_GTRCAT laurasiatheria_therest.aln ../names.txt ../lineage.txt > screen_logs/01_log1_fixspecies.txt

##>> NOTE: I needed to update the python script as it was not working. I just modifed `fix_alignment()` 
##>>       so the output file was different from the one being read

## Move unncessary output files to `00_filt1/` and rename output alignment
mkdir 00_filt1
mv _species_name_change.txt lineage.txt lineage.txt.1.bak 00_filt1/
rm laurasiatheria_therest.aln
mv laurasiatheria_therest_out.phylip alignment.phylip

# 2. Remove species 

# 2.1. Create input file for python script `prune_alignment.py`. This text file 
#      contains a list of the species that need to be removed 
cat > _species_remove.txt << EOF
felis_silvestris_silvestris
canis_lupus_lupus
bassaricyon_gabbii_gabbii
martes_melampus_melampus
meles_meles_meles
ursus_thibetanus_thibetanus
odobenus_rosmarus_rosmarus
pusa_hispida_hispida
phoca_vitulina_vitulina
dicerorhinus_sumatrensis_sumatrensis
rhinoceros_sondaicus_sondaicus
ceratotherium_simum_simum
equus_hemionus_hemionus
crocidura_suaveolens_suaveolens
suncus_murinus_murinus
neomys_anomalus_anomalus
sorex_hoyi_hoyi
sorex_thompsoni_thompsoni
sorex_palustris_palustris
neofelis_nebulosa_nebulosa
panthera_leo_leo
panthera_tigris_tigris
urocyon_littoralis
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
../../../../../../src/newick-utils-1.6/src/nw_prune RAxML_bestTree.BS_ML_GTRCAT $sp2rm > tmp.tree
## COOL! Note that `urocyon_littoralis` is not found because it had already been removed 
## from the data set before starting this filtering step. Therefore, the unpruned tree did no 
## have this label either. 

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
   in RStudio and change line 24 so it is `subt    <- "laurasiatheria_therest"` and uncomment line 26. 
   Now, we can run it from RStudio. 
   This script will generate a concatenated alignment file with all partitions, as well as 
   one alignment file for each individual partitions,
   inside [`01_alignments`](/02_SeqBayes_S2/00_Data_filtering/01_alignments).
   Log files and Rdata can be found [here](/02_SeqBayes_S2/00_Data_filtering/01_alignments/Rout/log_concatenation). 
      > NOTE: Updated `partitions.txt` file inside `00_mammal_alns/laurasiatheria_therest/`.   
   * Paths have been automatically set according to current file architecture in `mammals` dir, 
   do not change paths in the Rscript! 

The final alignments generated at the end of this step can be found in the directory 
`before_adding_new_taxa` when downloading 
[this zip file](https://www.dropbox.com/s/f5ppv0hn168xnoq/SeqBayesS2_Raln_laurasiatheria_therest_ALLFILTERS.zip?dl=0).
Please continue reading how to generate the final alignments in the next section!

# EXTRA FILTERING -- ADDING TAXA TO THE ALIGNMENT
When we first carried out the analysis with the tree topology as in step 1, we realised that 
we should include the calibrations for nodes "Chiroptera", "Euungulata", and "Artiodactyla"
to avoid future issues when grafting the 
data subtrees onto the backbone tree at the end of the sequential Bayesian dating approach. 

**IMPORTANT NOTE**: this extra filtering step requires to have the first generated alignments with the R script
for the "Chiroptera" subtree, otherwise you will not be able to proceed!

The procedure followed was the following:

## 1. Generate data subsets
[Here](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/laurasiatheria_therest/filter_aln/extra_filtering),
in section `1. Generate fasta files to generate alignments with new taxa`, you will find 
all the details about how to reformat the alignments for this data subset before we 
included the data for the new taxa. You can download all the reformatted alignments generated as well 
as the data needed to run MAFFT for the re-alignment
[here](https://www.dropbox.com/s/6hvirfh37w8dw17/SeqBayesS2_filteraln2_ltherest_00_perl_parsing.zip?dl=0).
The content of this zip file should be saved inside directory
[`00_perl_parsing`](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/laurasiatheria_therest/filter_aln/extra_filtering/00_perl_parsing),
where you can find all the scripts used in this filtering step.

## 2. Generate new alignments with new taxa 
[Here](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/laurasiatheria_therest/filter_aln/extra_filtering),
in section `2. Getting alignment with new taxa`, you will find 
all the details to generate the final alignments. 
The directory
[`01_MAFFT`](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/laurasiatheria_therest/filter_aln/extra_filtering/01_MAFFT) 
contains the scripts used to generate the alignment that includes the sequences for 
*Pteropus vampyrus*, *Myotis lucifugus*, *Vicugna pacos*, and *Capra hircus*.
These are the species which sequence data were needed to include the calibrations for the nodes "Chiroptera", 
"Euungulata", and "Artiodactyla".
You can download the `maff_ltherest` directory
[here](https://www.dropbox.com/s/9kr87g1iddaqr8o/SeqBayesS2_filteraln2_ltherest_01_MAFFT.zip?dl=0),
which should be saved inside the 
[`01_MAFFT`](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/laurasiatheria_therest/filter_aln/extra_filtering/01_MAFFT) 
directory mentioned above when running the scripts and following the instructions detailed in the link 
provided above. 

## 3. Final alignments
[Here](https://www.dropbox.com/s/f5ppv0hn168xnoq/SeqBayesS2_Raln_laurasiatheria_therest_ALLFILTERS.zip?dl=0),
you can download the final alignments for this data set. The previous alignments have been moved to a directory called 
`before_adding_new_taxa`, which you can find inside this zip file too.
