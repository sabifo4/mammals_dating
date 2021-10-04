# RODENTIA THE REST - Filtering alignment
The initial files that you can find in this directory are the following:

```
rodentia_therest
   |- filter_aln
         |- checked_aln                        
         |     |- taxonomical_check
         |           |- rodentia_therest.png
         |           |- rodentia_therest_FigTree
         |           |- RAxML_bestTree.BS_ML_GTRCAT
         |     
         |- extra_filtering   <-- You will read more about it at the end of this guideline                     
         |     
         |- lineage.txt 
         |- names.txt       
         |- parse_lineage.R 
         |- partitions.txt 
         |- README.md
         |- rodentia_therest_taxonomy_check.xlsx
         |- summary.html 
```
		 
The output files during the filtering step that is described below as well as the input 
files needed have been zipped in a file as they are very large. You can download 
this file [here](https://www.dropbox.com/s/ksls0qu4xigxi80/SeqBayesS2_filtaln_rodentia_therest.zip?dl=0).
To start the filtering step, you should have the following files arranged in the file 
architecture detailed above (you can obtain the files once you unzip the file
provided in the link above): 

```
rodentia_therest 
   |- filter_aln
         |- checked_aln                        
         |     |- unfiltered_aln
         |     |    |- alignment_nt3cp.phylip   # Alignment with 3CP 
         |     |- RAxML_bestTree.BS_ML_GTRCAT   # Tree you can find in `original_tree` directory
         |     |- taxonomical_check             # Visual checks to evaluate 
         |          |- rodentia_therest.png     # dubious taxa placement
         |          |- rodentia_therest_FigTree
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
         |- rodentia_therest_taxonomy_check.xlsx
         |- summary.html 
```

Note that you will see more files within the zipped file provided above. These are output 
files that you will generate if you follow all the steps below. Feel free to keep them so you can 
make sure that you reproduce the same results that we show here. 

## 1. Taxonomic filtering 
Within the [`filter_aln`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_therest/filter_aln)
directory, you will find different files as detailed above. 
Specifically, the R script [`parse_lineage.R`](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_therest/filter_aln/parse_lineage.R)
was written to carry out a first taxonomic filtering. Note that this 
R script will run if you have the same file architecture in this GitHub repository (i.e., it uses 
a function within the R script [`Filter_lineages.R`](https://github.com/sabifo4/mammals_dating/blob/main/src/Filter_lineages.R)
in the `src` directory and the [`genes.txt`](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/genes.txt)
file. There are several messages printed out by this script as all the rodents are 
parsed at the same time. You might have seen that, due to the large number of rodent species, 
we have had to divide them into different data subsets: "rodentia_therest",
"rodentia_ctenohystrica", "rodentia_subtree1", and "rodentia_subtree2". This is "rodentia_therest", 
the first filtering step before obtaining the two data subset "rodentia_subtree1" and 
"rodentia_subtree".
Below, you will find the messages printed out that refer to this data subset, "rodentia_therest":   

```
FAMILY   Geomyidae  has  75  taxa. Species  geomys_personatus  has  1  genes that are not shared by any of the 
other taxa. It should be deleted!
GENUS   Geomys  has  31  taxa. Species  geomys_personatus  has  1  genes that are not shared by any of the
other taxa. It should be deleted!
```
>> ACTION: geomys_personatus was removed.

```
GENUS   Anomalurus  has  2  taxa. Species  anomalurus_beecrofti  has  14  genes that are not shared by any
of the other taxa. It should be deleted!
GENUS   Anomalurus  has  2  taxa. Species  anomalurus_sp_gp_2005  has  14  genes that are not shared by
any of the other taxa. It should be deleted!
```
>> ACTION: anomalurus_beecrofti & anomalurus_sp_gp_2005 were merged. The label is "anomalurus".

```
FAMILY   Calomyscidae  has  4  taxa. Species  calomyscus_mystax  has  3  genes that are not shared by any of
the other taxa. It should be deleted!
```
>> ACTION: calomyscus_mystax was removed. 

```
FAMILY   Cricetidae  has  661  taxa. Species  peromyscus_maniculatus_luteus  has  1  genes that are not shared
by any of the other taxa. It should be deleted!
```
>> ACTION: peromyscus_maniculatus_luteus was removed.

```
SUBFAMILY   Deomyinae  has  20  taxa. Species  acomys_russatus_russatus  has  1  genes that are not shared by
any of the other taxa. It should be deleted!
```
>> ACTION: acomys_russatus_russatus was removed.

```
SUBFAMILY   Calomyscinae  has  4  taxa. Species  calomyscus_mystax  has  3  genes that are not shared by any
of the other taxa. It should be deleted!
```
>> ACTION: calomyscus_mystax was removed. 

```
SUBFAMILY   Petromyscinae  has  2  taxa. Species  petromyscus_collinus  has  2  genes that are not shared by 
any of the other taxa. It should be deleted!
SUBFAMILY   Petromyscinae  has  2  taxa. Species  petromyscus_monticularis  has  3  genes that are not shared 
by any of the other taxa. It should be deleted!
```
>> ACTION: petromyscus_collinus and petromyscus_monticularis were merged. Label is "petromyscus".

```
SUBFAMILY   Neotominae  has  172  taxa. Species  peromyscus_maniculatus_luteus  has  1  genes that are not 
shared by any of the other taxa. It should be deleted!
```
>> ACTION: peromyscus_maniculatus_luteus was removed.

```
GENUS   Taterillus  has  4  taxa. Species  taterillus_emini  has  4  genes that are not shared by any of the
other taxa. It should be deleted!
```
>> ACTION: taterillus_emini was removed. 

```
GENUS   Acomys  has  15  taxa. Species  acomys_russatus_russatus  has  1  genes that are not shared by any of
the other taxa. It should be deleted!
```
>> ACTION: acomys_russatus_russatus was removed. 

```
GENUS   Pyromys  has  3  taxa. Species  mus_shortridgei  has  1  genes that are not shared by any of the other
taxa. It should be deleted!
```
>> ACTION: mus_shortridgei was removed. 

```
GENUS   Notomys  has  2  taxa. Species  notomys_alexis  has  6  genes that are not shared by any of the other
taxa. It should be deleted!
GENUS   Notomys  has  2  taxa. Species  notomys_fuscus  has  2  genes that are not shared by any of the other
taxa. It should be deleted!
```
>> ACTION: notomys_alexis & notomys_fuscus were merged. Label is "notomys". 

```
GENUS   Pseudomys  has  2  taxa. Species  pseudomys_australis  has  4  genes that are not shared by any of the
other taxa. It should be deleted!
GENUS   Pseudomys  has  2  taxa. Species  pseudomys_chapmani  has  14  genes that are not shared by any of the
other taxa. It should be deleted!
```
>> ACTION: pseudomys_australis & pseudomys_chapmani were merged. Label is "pseudomys". 

```
GENUS   Abeomelomys  has  2  taxa. Species  abeomelomys_sevia  has  3  genes that are not shared by any of the
other taxa. It should be deleted!
GENUS   Abeomelomys  has  2  taxa. Species  abeomelomys_sevia_tatei  has  1  genes that are not shared by any of
the other taxa. It should be deleted!
```
>> ACTION: abeomelomys_sevia & abeomelomys_sevia_tatei are flagged. We decided to remove 
>> "abeomelomys_sevia_tatei" because it had less genes than "abeomelomys_sevia".
>> CHECK USING THE R OBJECT:
>> ```abeomelomys_sevia       --> "CO2" "ENSG00000012048" "ENSG00000112964"```
>> ```abeomelomys_sevia_tatei --> "CYB"```

```
GENUS   Lemniscomys  has  7  taxa. Species  lemniscomys_barbarus_barbarus  has  1  genes that are not shared by
any of the other taxa. It should be deleted!
```
>> ACTION: lemniscomys_barbarus_barbarus was removed. 

```
GENUS   Calomyscus  has  4  taxa. Species  calomyscus_mystax  has  3  genes that are not shared by any of the
other taxa. It should be deleted!
```
>> ACTION: calomyscus_mystax was removed. 

```
GENUS   Petromyscus  has  2  taxa. Species  petromyscus_collinus  has  2  genes that are not shared by any of
the other taxa. It should be deleted!
GENUS   Petromyscus  has  2  taxa. Species  petromyscus_monticularis  has  3  genes that are not shared by any
of the other taxa. It should be deleted!
```
>> ACTION: petromyscus_collinus and petromyscus_monticularis were merged. Label is "petromyscus".

```
GENUS   Macrotarsomys  has  2  taxa. Species  macrotarsomys_bastardi  has  3  genes that are not shared by
any of the other taxa. It should be deleted!
GENUS   Macrotarsomys  has  2  taxa. Species  macrotarsomys_ingens  has  2  genes that are not shared by any
of the other taxa. It should be deleted!
```
>> ACTION: macrotarsomys_bastardi and macrotarsomys_ingens were merged. Label is "macrotarsomys".

```
GENUS   Peromyscus  has  100  taxa. Species  peromyscus_maniculatus_luteus  has  1  genes that are not shared by
any of the other taxa. It should be deleted!
```
>> ACTION: peromyscus_maniculatus_luteus was removed. 

```
GENUS   Alticola  has  15  taxa. Species  alticola_semicanus_alleni  has  1  genes that are not shared by any of the other taxa. It should be deleted!
GENUS   Alticola  has  15  taxa. Species  alticola_stracheyi  has  1  genes that are not shared by any of the other taxa. It should be deleted!
```
>> ACTION: alticola_semicanus_alleni & alticola_stracheyi were removed. 

```
GENUS   Rhagomys  has  2  taxa. Species  rhagomys_rufescens  has  1  genes that are not shared by any of the
other taxa. It should be deleted!
GENUS   Rhagomys  has  2  taxa. Species  rhagomys_longilingua  has  1  genes that are not shared by any of the
other taxa. It should be deleted!
```
>> ACTION: rhagomys_rufescens and rhagomys_longilingua were merged. Label is "rhagomys".

```
GENUS   Brucepattersonius  has  2  taxa. Species  brucepattersonius_soricinus  has  2  genes that are not shared
by any of the other taxa. It should be deleted!
GENUS   Brucepattersonius  has  2  taxa. Species  brucepattersonius_igniventris  has  2  genes that are not
shared by any of the other taxa. It should be deleted!
```
>> ACTION: brucepattersonius_soricinus and brucepattersonius_igniventris were merged. Label is "brucepattersonius".

```
GENUS   Bibimys  has  2  taxa. Species  bibimys_chacoensis  has  1  genes that are not shared by any of the
other taxa. It should be deleted!
GENUS   Bibimys  has  2  taxa. Species  bibimys_labiosus  has  1  genes that are not shared by any of the other
taxa. It should be deleted!
```
>> ACTION: bibimys_chacoensis and bibimys_labiosus were merged. Label is "bibimys".

```
GENUS   Rhipidomys  has  5  taxa. Species  rhipidomys_leucodactylus  has  1  genes that are not shared by any of
the other taxa. It should be deleted!
```
>> ACTION: rhipidomys_leucodactylus was removed.

**NOTE**   
Log files will be found in this same directory with the name `log_taxaNOTIN_<level>.txt`,
being `level` the one that has been checked (i.e., family, genus, order, or subfamily). In addition, 
you can also generate the `levels.checked.RData`, which you can use to explore the taxonomical 
levels. Remember that you will see more information about the other rodent taxa that are 
present in the other three data subsets.

# 2. First checks before applying filtering
First, we checked if there were any further species that should be removed:

```sh
# Run from `rodentia/checked_aln/` directory 
grep -o '[a-z]*_[a-z]*_[a-z]*' RAxML_bestTree.BS_ML_GTRCAT | sed 's/\_/\t/g' | awk '$2==$3{print $1,$2,$3}'  > subsp_check.txt

# Count 
num=$( grep -o '[a-z]*_[a-z]*_[a-z]*' RAxML_bestTree.BS_ML_GTRCAT | sed 's/\_/\t/g' | awk '$2==$3{print $1,$2,$3}'  | awk '$1!="aeb"{print $1,$2,$3}' | wc -l )
printf "There are $num subspecies with species\n"
## There are 39 subspecies with species
```
The output names have been saved in a file called `subsp_check.txt` to further explore it:

```sh
# Run from `00_data_curation/rodentia_therest/filter_aln/checked_aln` directory 
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
   Num. sequences oryzomys_couesi : 6   
   Num. sequences **oryzomys_couesi_couesi** : 2

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences oryzomys_palustris : 8   
   Num. sequences **oryzomys_palustris_palustris** : 2

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences abrothrix_olivaceus : 4   
   Num. sequences **abrothrix_olivaceus_olivaceus** : 2   
   ```
   abrothrix_olivaceus_olivaceus --> "CYB"
   abrothrix_olivaceus           --> "CYB" "ENSG00000152592" "ENSG00000265203"
   ```

> **REMOVE SUBSPECIES -- same amount of genes, we keep species**   
   Num. sequences phyllotis_xanthopygus : 3   
   Num. sequences **phyllotis_xanthopygus_xanthopygus** : 3   
   ```
   phyllotis_xanthopygus             --> "CYB" "ENSG00000265203"
   phyllotis_xanthopygus_xanthopygus --> "CYB" "ENSG00000166349"
   ```
   
> **REMOVE SUBSPECIES, less genes**   
   Num. sequences phyllotis_osilae : 5   
   Num. sequences **phyllotis_osilae_osilae** : 1

> **REMOVE SUBSPECIES -- same amount of genes, we keep species**   
   Num. sequences akodon_subfuscus : 2   
   Num. sequences **akodon_subfuscus_subfuscus** : 2   
   ```
   akodon_subfuscus_subfuscus --> "CYB"
   akodon_subfuscus           --> "CYB"
   ```
   
> **REMOVE SUBSPECIES -- same amount of genes, we keep species**   
   Num. sequences akodon_azarae : 3   
   Num. sequences **akodon_azarae_azarae** : 3   
   ```
   akodon_azarae_azarae --> "CYB" "ENSG00000152592"
   akodon_azarae        --> "CYB" "ENSG00000265203"
   ```

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences alticola_strelzowi : 3   
   Num. sequences **alticola_strelzowi_strelzowi** : 2   
   ```
   alticola_strelzowi_strelzowi --> "CYB"
   alticola_strelzowi           --> "CYB" "ENSG00000265203"
   ```

> **REMOVE SPECIES, less genes**   
   Num. sequences **microtus_fortis** : 6   
   Num. sequences microtus_fortis_fortis : 26   
   ```
   microtus_fortis_fortis --> "ATP6" "ATP8" "CO1"  "CO2"  "CO3"  "CYB"  "ND1"  "ND2" 
                              "ND3"  "ND4"  "ND4L" "ND5"  "RNR1" "RNR2"
   microtus_fortis        --> "CO1" "CYB" "ENSG00000112964" "ENSG00000265203"
   ```

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences cricetulus_migratorius : 8   
   Num. sequences **cricetulus_migratorius_migratorius** : 1

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences nyctomys_sumichrasti : 5   
   Num. sequences **nyctomys_sumichrasti_sumichrasti** : 2

> **REMOVE SUBSPECIES**   
   Num. sequences neotoma_albigula : 6   
   Num. sequences **neotoma_albigula_albigula** : 2

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences peromyscus_pectoralis : 3   
   Num. sequences **peromyscus_pectoralis_pectoralis** : 2   
   ```
   peromyscus_pectoralis_pectoralis --> "CYB"
   peromyscus_pectoralis            --> "CYB" "ENSG00000265203"
   ```

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences peromyscus_truei : 4   
   Num. sequences **peromyscus_truei_truei** : 2   
   ```
   peromyscus_truei_truei --> "CYB"
   peromyscus_truei       --> "CYB" "ENSG00000258839" "ENSG00000265203"
   ```
   
> **REMOVE SUBSPECIES, less genes**   
   Num. sequences peromyscus_gratus : 3   
   Num. sequences ** peromyscus_gratus_gratus**  : 2   
   ```
   peromyscus_gratus_gratus --> "CYB"
   peromyscus_gratus        --> "CYB" "ENSG00000265203"
   ```

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences peromyscus_difficilis : 4   
   Num. sequences ** peromyscus_difficilis_difficilis**  : 2   
   ```
   peromyscus_difficilis_difficilis --> "CYB"
   peromyscus_difficilis            --> "CYB" "ENSG00000152592" "ENSG00000258839"
   ```

> **REMOVE SUBSPECIES -- same amount of genes, we keep species**   
   Num. sequences peromyscus_nasutus : 2   
   Num. sequences ** peromyscus_nasutus_nasutus**  : 2   
   ```
   peromyscus_nasutus         --> "CYB" 
   peromyscus_nasutus_nasutus --> "CYB"
   ```

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences peromyscus_boylii : 8   
   Num. sequences **peromyscus_boylii_boylii** : 4

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences peromyscus_levipes : 4   
   Num. sequences **peromyscus_levipes_levipes** : 2   
   ```
   peromyscus_levipes_levipes --> "CYB"
   peromyscus_levipes         --> "CYB" "ENSG00000112964" "ENSG00000265203"
   ```

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences peromyscus_polionotus : 12   
   Num. sequences **peromyscus_polionotus_polionotus** : 3

> **NO ACTION, spalax_graecus had been already deleted!**   
   ~Num. sequences **spalax_graecus** : 2~   
   Num. sequences spalax_graecus_graecus : 4   
   ```
   spalax_graecus_graecus --> "ND1"  "RNR1" "RNR2"
   spalax_graecus         --> "CYB"   
   ```   
   
> **REMOVE SPECIES, less genes**   
   Num. sequences **sicista_subtilis** : 2   
   Num. sequences sicista_subtilis_subtilis : 3   
   ```
   sicista_subtilis_subtilis --> "CYB" "ENSG00000265203"
   sicista_subtilis          --> "CYB"
   ```

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences napaeozapus_insignis : 8   
   Num. sequences **napaeozapus_insignis_insignis** : 3

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences zapus_hudsonius : 7   
   Num. sequences **zapus_hudsonius_hudsonius** : 3

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences zapus_trinotatus : 5   
   Num. sequences **zapus_trinotatus_trinotatus** : 3   
   ```
   zapus_trinotatus_trinotatus --> "CYB" "ENSG00000012048"
   zapus_trinotatus            --> "CYB" "ENSG00000012048" "ENSG00000166349" "ENSG00000265203"
   ```

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences zapus_princeps : 6   
   Num. sequences **zapus_princeps_princeps** : 3

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences cratogeomys_castanops : 17   
   Num. sequences **cratogeomys_castanops_castanops** : 2

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences cratogeomys_goldmani : 4   
   Num. sequences **cratogeomys_goldmani_goldmani** : 2   
   ```
   cratogeomys_goldmani_goldmani  --> "CYB"
   cratogeomys_goldmani           --> "CO1" "CYB"
   ```

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences cratogeomys_fumosus : 4   
   Num. sequences **cratogeomys_fumosus_fumosus** : 2   
   ```
   cratogeomys_fumosus_fumosus  --> "CYB"
   cratogeomys_fumosus          --> "CO1" "CYB"
   ```

> **REMOVE SUBSPECIES -- same amount of genes, we keep species**   
   Num. sequences geomys_pinetis : 3   
   Num. sequences **geomys_pinetis_pinetis** : 3   
   ```
   geomys_pinetis_pinetis --> "CYB" "ENSG00000265203"
   geomys_pinetis         --> "CYB"  "RNR1"
   ```

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences geomys_breviceps : 13   
   Num. sequences **geomys_breviceps_breviceps** : 3   

> **NO ACTION, geomys_presonatus had been already deleted!**    
    ~Num. sequences **geomys_personatus** : 2~   
    Num. sequences geomys_personatus_personatus : 4   
   ```
   geomys_personatus_personatus --> "CYB" "RNR1" "ENSG00000265203"   
   geomys_personatus            --> "ND5"   
   ```   
	
> **REMOVE SUBSPECIES, less genes**   
   Num. sequences geomys_bursarius : 5   
   Num. sequences **geomys_bursarius_bursarius** : 2   
   ```
   geomys_bursarius_bursarius --> "CYB"
   geomys_bursarius           --> "CO2" "RNR1" "ENSG00000012048" "ENSG00000112964"
   ```

> NO ACTION   
   Num. sequences geomys_jugossicularis : 0   
   Num. sequences geomys_jugossicularis_jugossicularis : 4

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences castor_fiber : 28   
   Num. sequences **castor_fiber_fiber** : 2

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences rattus_leucopus : 29   
   Num. sequences **rattus_leucopus_leucopus** : 1

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences mus_musculus : 194   
   Num. sequences **mus_musculus_musculus** : 28

> **REMOVE SUBSPECIES, less genes**   
   Num. sequences otomys_denti : 3   
   Num. sequences **otomys_denti_denti** : 2   
   ```
   otomys_denti_denti --> "CYB"
   otomys_denti       --> "CYB"  "RNR1"
   ```
   
> **REMOVE SUBSPECIES, less genes**   
   Num. sequences otomys_typus : 3   
   Num. sequences **otomys_typus_typus** : 2   
   ```
   otomys_typus_typus --> "CYB"
   otomys_typus       --> "CYB"  "RNR1"
   ```

At this stage, there are also 3 duplicates. We check which genes are in `trinomys_setosus_denigratus`
using R object `levels.checked.RData` previously generated:

```R
# Open Rscript `parse_lineage.R` and run from lines 1-15.
# Then, uncomment line 67 and run it to load the R object 

load( "levels.checked.RData" )

# Check "habromys_delicatulus_duplicate"
levels.checked$genus$Peromyscus$peromyscus_sagax
# [1] "CYB"
levels.checked$genus$Habromys$habromys_delicatulus
# [1] "CYB"
# > Both of them have the same gene, so we keep habromys_delicatus 
# > and remove the duplicate 

# Check "neodon_sikimensis_duplicate"
levels.checked$genus$Microtus$microtus_sikimensis
# [1] "ENSG00000265203"
levels.checked$genus$Neodon$neodon_sikimensis
# [1] "CO1"             "CYB"             "ENSG00000112964"
# > We remove the duplicate as it only has 1 gene

# Check "rattus_andamanensis_duplicate"
levels.checked$genus$Rattus$rattus_tanezumi_sladeni
# [1] "CYB"
levels.checked$genus$Rattus$rattus_andamanensis
[# 1] "CYB" ="ENSG00000112964" "ENSG00000136869" "ENSG00000196664" "ENSG00000265203"
# > We remove the duplicate as it only has 1 gene
```

The updated list of taxa to be renamed and removed at this stage is the following:

```
Names to change:
eothenomys_smithii,myodes_smithii
akodon_latebricola,neomicroxus_latebricola
oxymycterus_iheringi,brucepattersonius_iheringi
bolomys_urichi,necromys_urichi
bolomys_amoenus,necromys_amoenus
bolomys_temchuki,necromys_temchuki
bolomys_lasiurus,necromys_lasiurus
muriculus_imberbis,mus_imberbis
sigmodontomys_aphrastus,tanyuromys_aphrastus
phyllotis_wolffsohni,tapecomys_wolffsohni
cricetulus_kamensis,urocricetus_kamensis
cricetulus_lama,urocricetus_lama
spalax_carmeli,nannospalax_carmeli
geomys_bursarius_lutescens,geomys_lutescens
geomys_bursarius_halli,geomys_jugossicularis_halli
perognathus_flavescens_apache,perognathus_apache
hylomyscus_denniae_vulcanorum,hylomyscus_vulcanorum
peromyscus_mexicanus_nicaraguae,peromyscus_nicaraguae
peromyscus_mexicanus_salvadorensis,peromyscus_salvadorensis
peromyscus_mexicanus_tropicalis,peromyscus_tropicalis
meriones_meridianus_penicilliger,meriones_penicilliger
otomys_orestes_zinki,otomys_zinki
perognathus_parvus_mollipilosus,perognathus_mollipilosus
thomomys_anitae,thomomys_bottae_anitae
chaetodipus_arenarius_siccus,chaetodipus_siccus
lasiopodomys_gregalis,microtus_gregalis

Taxa to remove:
rattus_leucopus_leucopus
peromyscus_difficilis_difficilis
oryzomys_couesi_couesi
cratogeomys_fumosus_fumosus
sicista_subtilis
phyllotis_osilae_osilae
cricetulus_migratorius_migratorius
phyllotis_xanthopygus_xanthopygus
abrothrix_olivaceus_olivaceus
apodemus_fulvipectus
rattus_tanezumi_sladeni
microtus_sikimensis
peromyscus_sagax
```

In addition, we carried out other checks with regards to the species names and their 
placement in the phylogeny. No taxa were indicated there to be removed, only some 
taxa are to be renamed. You can find a list of these checks in the file 
`rodentia_therest_taxonomy_check.csv` as well as below:   

```
Genus       Rheomys                  CHANGE PLACEMENT      rheomys_thomasi                                                   --> Ref. suggests Rhemoys species cluster together. They are sister taxa of ichithyomys stolzmanni and part of a paraphyletic clade with Sigmodon species.
Genus       Reithrodontomys          CHANGE PLACEMENT      peromyscus_polius                                                 --> Reference finds peromyscus clusters with the rest of peromyscus species as an outter group before Onychomys species closes this paraphyletic clade. Check if we want to move this taxa.
Genus       Handleyomys              CHANGE PLACEMENT      handleyomys_intectus                                              --> Seems to be in wrong placement. Move as the outter group of paraphyletic clade with Oecomys, Euryoryzomys, Transandinomys, Handleyomys, Hylaemys, and Nephelomys species. 
Genus       Myodes                   CHECK PLACEMENT       alticola_macrotis_vinogradovi, myodes_centralis                   --> Ref. A shows a similar clustering except for i) the relationships within subclades, ii) M. centralis and A. macrotis vinogradovi do not cluster in same clade but with M. rufocanus,
                                                                                                                                 iii) E. Smithii is M. smithii. Ref. B (which still uses synonym Craseomys instead of Myodes) shows M. smiithi (there as C. smiithi) clustering with M. andersoni (C. andersoni).
																																 Ref. c. shows "((M. centralis, A. macrotis), M. glareolus);". Check placement and decide if we move flagged outside taxa to match references.
Genus       Meriones                 CHECK PLACEMENT       meriones_tamariscinus                                             --> Reference found two inferred topologies: one as ours and the other having M. tamariscinus as the outer taxa of a monophyletic clade for Meriones. Decide what to do.
Genus       Orthogeomys              CHECK PLACEMENT       orthogeomys_grandis                                               --> Reference A found same clustering of all Orthogeomys sp. but the flagged one as a unique sister clade to another clade in which Pappogeomys, Cratogeomys, and O. grandis are.
                                                                                                                                 The difference is that we found O. grandis as an outgroup of these two subclades. Our topology: "((Orthogeomys sp., (Cratogeomys sp., Pappogeomys bulleri)), O. grandis);".
																																 Resolved phylogeny in ref: "(Orthogeomys sp., ((Cratogeomys sp., Pappogeomys bulleri), O. grandis));". See also recovered topology in ref. B with incongruences about its placement. Decide what we do with placement
Genus       Bolomys                  CHECK PLACEMENT       podoxymys_roraimae                                                --> Reference suggests that it is the outter group of the clade necromys, actually clustering with Thalpomys cerradensis.
Genus       Arvicanthis              CHECK PLACEMENT       arvicanthis_ansorgei                                              --> Reference found Lemniscomys species clustering with Arvicanthis species in a paraphyletic clade as we found. They did not include A. ansorgei, but it seems they should cluster together
                                                                                                                                 as it is shown in reference B with more Arvicanthis species. CHECK placement and if they should form monophyletic.
Genus       Allactaga                CHECK PLACEMENT       pygeretmus_pumilio                                                --> Reference shows that it tends to cluster with A. elater, but it changes depending on the gene used to infer topology. Decide if we change its placement.
Genus       Acomys                   CHECK PLACEMENT       acomys_subspinosus                                                --> Seems wrong placement. Reference places it with the rest of Acomys species, as the outter group of all of them. Check placement.
Genus       Grammomys                CHECK PLACEMENT       thallomys_loringi,thallomys_nigricauda,thallomys_paedulcus        --> Reference finds Grammomys form a monophyletic clade and Thallomy species cluster with Micaelamys namaquensis. Check placement.
Genus       Rattus                   CHECK PLACEMENT       diplothrix_legata,bandicota_savilei,nesokia_indica,               --> Reference finds the same clade of species but it seems it is placed as the outter clade for the rest of rattus species plus the clade with tarsomys, limnomys and rattus everetti. Check if we want to fix it like this.
                                                           bandicota_bengalensis,bandicota_indica
Genus       Synaptomys               CHECK PLACEMENT       myopus_schisticolor                                               --> Reference finds it clusters with Lemmus species, not with Synaptomys species. Decide if we rearrange taxa placement.
Genus       Parotomys                CHECK PLACEMENT       parotomys_brantsii,parotomys_littledalei                          --> Reference finds a paraphyletic clade with Myotomys and Parotomys species. They do not have parotomys_littledalei, but it seems to be sister species with P. brantsii according to MSW. So maybe just place them together in our topology.
Genus       Myomyscus                CHECK PLACEMENT       myomyscus_yemeni                                                  --> Reference finds they are sister taxa and cluster in a paraphyletic clade with Stenocephalemys species. Check placement and decide if we move myomyscus_yemeni as sister taxa of myomyscus_brockmani.
Genus       Neacomys                 CHECK PLACEMENT       neacomys_minutus                                                  --> Reference finds that Oreoryzomys balneator and Microryzomys  are sister taxa and part of a paraphyletic clade with Neacomys species. We can keep our topology as it is or try to move neacomys_minutus with the rest of neacomys species.
Genus       Aethomys                 CHECK PLACEMENT       dephomys_defua                                                    --> Reference suggests it clusters with Hybomys species and Stochomys longicaudatus (sister taxa of the latter). Check placement and decide if we move it here.
Genus       Dendromus                CHECK PLACEMENT       malacothrix_typica                                                --> Reference finds M. typica is the outgroup for Dendromus species. Check placement
Species     Peromyscus boylii        CHANGE PLACEMENT      peromyscus_boylii                                                 --> Move species with subspecies clade.
Species     Otomys tropicalis        CHANGE PLACEMENT      otomys_tropicalis                                                 --> Move specise with subspecies clade
Species     Myodes rufocanus         CHANGE PLACEMENT      myodes_rufocanus                                                  --> Reference shows myodes rufocanus clustering where the subsp. Is. Move species with subsp.
Species     Peromyscus crinitus      CHANGE PLACEMENT      peromyscus_crinitus                                               --> Not much about subspecies, but it seems subspecies clusters in clade according to reference and species should be moved there,
Species     Phyllotis osilae         CHANGE PLACEMENT      phyllotis_osilae_phaeus                                           --> Move subspecies to cluster with species
Species     Cricetulus migratorius   CHANGE PLACEMENT      cricetulus_migratorius                                            --> It seems to be a genetic differentiation between subspecies C. migratorius migratorious and C. phaeus, but as we are removing C. migratorius migratorius,  and they used to cluster together, 
                                                                                                                                 we should move species C. migratorius to cluster with subsp. C. migatorius phaeus as ref. B seems to support it cluster within this subclade.
Species     Dipodomys merriami       CHECK PLACEMENT       dipodomys_insularis                                               --> D. insularis seems to be sister taxa of D. melanurus. All D. merriami seem to form a monophyletic group. Ref. b finds sth similar as our topology, though. Check and decide if changed.
Species     Dipodomys simulans       CHECK PLACEMENT       dipodomys_simulans_peninsularis                                   --> Not much information. Keep as it is or try to cluster sp. With subsp.
Species     Neotoma ferruginea       CHECK PLACEMENT       neotoma_isthmica                                                  --> Keep as it is or cluster sp. N. ferruginea with subsp. so N. isthmica is outter taxa
Species     Oryzomys couesi          CHECK PLACEMENT       oryzomys_mexicanus                                                --> Try to cluster sp. With subsp.
Species     Cratogeomys fumosus      CHECK PLACEMENT       cratogeomys_fumosus_imparilis,cratogeomys_fumosus_angustirostris  --> Reference shows that all C. fumosus sequences cluster together. Keep as it is our topology or try to match subsp. With sp.
Species     Neotoma lepida           CHECK PLACEMENT       neotoma_lepida                                                    --> Not much information. Keep as it is or move neotoma lepida with its subsp to be sister taxa of N. bryanti as in ref.
Species     Peromyscus aztecus       CHECK PLACEMETN       peromyscus_aztecus_oaxacensis                                     --> Not much information as the ML tree inferred did not contain P. mexicanus. Maybe keep as it is.
```

Then, we opened the `RAxML_bestTree.BS_ML_GTRCAT` with `FigTree` to individually check the taxa with dubious 
plaecment. Those that are to be deleted are coloured in red, while those that they cluster
within the same clade with the species they belong to but the placement might be dubious have
been coloured in orange. The names for the latter will be changed so they have a `~` at the end. 
This means that, even though this is the placement found in the ML tree, other research have
found they cluster elsewhere. Therefore, we flagged them for future and further analyses that might be carried out,
but that are not part of the scope of this project. 
The output files (`png` and `FigTree` files) have been saved in the
[`checked_aln/taxonomical_check`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_therest/filter_aln/checked_aln/taxonomical_check)
directory. 

The taxa to be renamed and/or removed are the following:
```
Taxa to be renamed: 
rheomys_raptor,rheomys_raptor~
alticola_macrotis_vinogradovi,alticola_macrotis_vinogradovi~
myodes_centralis,myodes_centralis~
meriones_tamariscinus,meriones_tamariscinus~
orthogeomys_grandis,orthogeomys_grandis~
pygeretmus_pumilio,pygeretmus_pumilio~
thallomys_loringi,thallomys_loringi~
thallomys_nigricauda,thallomys_nigricauda~
thallomys_paedulcus,thallomys_paedulcus~
diplothrix_legata,diplothrix_legata~
bandicota_savilei,bandicota_savilei~
nesokia_indica,nesokia_indica~
bandicota_bengalensis,bandicota_bengalensis~
bandicota_indica,bandicota_indica~
parotomys_brantsii,parotomys_brantsii~
parotomys_littledalei,parotomys_littledalei~
neacomys_minutus,neacomys_minutus~
dephomys_defua,dephomys_defua~
otomys_tropicalis,otomys_tropicalis~
myodes_rufocanus,myodes_rufocanus~
peromyscus_crinitus,peromyscus_crinitus~
cricetulus_migratorius,cricetulus_migratorius~
cricetulus_migratorius_phaeus,cricetulus_migratorius_phaeus~
dipodomys_insularis,dipodomys_insularis~
dipodomys_simulans_peninsularis,dipodomys_simulans_peninsularis~
neotoma_isthmica,neotoma_isthmica~
oryzomys_mexicanus,oryzomys_mexicanus~
cratogeomys_fumosus_imparilis,cratogeomys_fumosus_imparilis~
cratogeomys_fumosus_angustirostris,cratogeomys_fumosus_angustirostris~
neotoma_lepida,neotoma_lepida~
peromyscus_aztecus_oaxacensis,peromyscus_aztecus_oaxacensis~

Taxa to be removed:
rheomys_thomasi
peromyscus_polius
handleyomys_intectus # Note: It is placed with Wilfredomys pictipes and Microakodontomys transitorius,
                     #       which seem to be missplaced too. W. pictipes seems to be mislabelled as 
					 #       there is a sequence labelled as Juliomys pictipes (homotypic synonym) and 
					 #       branch is very long here, hence these two taxa will be deleted too.
wilfredomys_pictipes
microakodontomys_transitorius
podoxymys_roraimae
acomys_subspinosus
myopus_schisticolor
myomyscus_yemeni
malacothrix_typica
peromyscus_boylii
phyllotis_osilae_phaeus
arvicanthis_ansorgei
```
   
   
**NOTE:** Data subsets (both alignments, 12CP and 3CP) for Afrotheria, Xenarthra, Euarchonta, and Marsupialia had already undergone 
these checks before 2018 (i.e., they had already been "cleaned", while the data subset for Lagomorpha 
had not yet).
Therefore, you do not see the last part of the filtering described in this section 
in the corresponding `README.md` files for these data subsets or the directory
[`taxonomical_check`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_therest/filter_aln/checked_aln/taxonomical_check)
You will see that the `filter_aln` directory for these four data subsets
contains a csv file with the details that were followed to filter the corresponding alignments (12CP and 3CP). 
Note that this csv file is equivalent to the excel sheet you find in this directory for 
data subset Rodentia the rest.

## 2. Check alignment with 3CP partition
Before we concatenate the alignment with the first and second codon positions (12CP), file `alignment.phylip`
(you will find it once you unzip the file which link is provided above),
and the alignment with the third codon
positions (3CP), file `alignment_nt3cp.phylip` (you will find it once 
you unzip the file which link provided above, inside a directory called `checked_aln` directory),
we ran the following code to make sure they had both undergone the same filtering steps and that were 
at the same "filtering stage":


```sh
# Run the next code from `rodentia_therest/filter_aln/checked_aln/unfiltered_aln`

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
   in RStudio and change line 25 so it is `subt    <- "rodentia_therest"` and uncomment line 27.
   This script will generate a concatenated alignment file with all partitions, as well as 
   one alignment file for each individual partitions, inside a new dir called `00_mammal_alns/rodentia_therest/unfiltered/`
   inside [`00_Data_filtering/01_alignments/`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/01_alignments). 
   Log and RData files will be saved inside
   [`00_Data_filtering/01_alignments/Rout`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/01_alignments/Rout/log_concatenation).   
      > NOTE 1: Updated `partitions.txt` file inside `00_mammal_alns/rodentia_therest/unfiltered/` generated.   
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
Run the following commands from `rodentia_therest/filter_aln/checked_aln`:

```sh
# Run the next code from `rodentia_therest/checked_aln/checked_aln`

# 1. Fix species names in `names.txt`. We used the RAxML file, the `names.txt` file, and 
#    the `lineages.txt` that had been generated when filtering the alignment 12CP.
#    The ML tree file is now saved inside `filter_aln/checked_aln`, from where the 
#    python script will be run. The `names.txt` and the `lineage.txt` files are in
#    the `filter_aln` dir, one up. The alignment is the concatenated one that was generated in the 
#    previous step and can be found inside
#    `01_alignments/00_mammal_alns/rodentia_therest/unfiltered/rodentia_therest.aln`.
#
#    Note that we just copy this file inside `checked_aln` to proceed with the filtering. We also generate 
#    a copy of the original tree file
cp ../../../../01_alignments/00_mammal_alns/rodentia_therest/unfiltered/rodentia_therest.aln .
mkdir original_tree 
cp RAxML_bestTree.BS_ML_GTRCAT original_tree/RAxML_bestTree.BS_ML_GTRCAT
cat > _species_name_change.txt << EOF
eothenomys_smithii,myodes_smithii
akodon_latebricola,neomicroxus_latebricola
oxymycterus_iheringi,brucepattersonius_iheringi
bolomys_urichi,necromys_urichi
bolomys_amoenus,necromys_amoenus
bolomys_temchuki,necromys_temchuki
bolomys_lasiurus,necromys_lasiurus
muriculus_imberbis,mus_imberbis
sigmodontomys_aphrastus,tanyuromys_aphrastus
phyllotis_wolffsohni,tapecomys_wolffsohni
cricetulus_kamensis,urocricetus_kamensis
cricetulus_lama,urocricetus_lama
spalax_carmeli,nannospalax_carmeli
geomys_bursarius_lutescens,geomys_lutescens
geomys_bursarius_halli,geomys_jugossicularis_halli
perognathus_flavescens_apache,perognathus_apache
hylomyscus_denniae_vulcanorum,hylomyscus_vulcanorum
peromyscus_mexicanus_nicaraguae,peromyscus_nicaraguae
peromyscus_mexicanus_salvadorensis,peromyscus_salvadorensis
peromyscus_mexicanus_tropicalis,peromyscus_tropicalis
meriones_meridianus_penicilliger,meriones_penicilliger
otomys_orestes_zinki,otomys_zinki
perognathus_parvus_mollipilosus,perognathus_mollipilosus
thomomys_anitae,thomomys_bottae_anitae
chaetodipus_arenarius_siccus,chaetodipus_siccus
lasiopodomys_gregalis,microtus_gregalis
rheomys_raptor,rheomys_raptor~
alticola_macrotis_vinogradovi,alticola_macrotis_vinogradovi~
myodes_centralis,myodes_centralis~
meriones_tamariscinus,meriones_tamariscinus~
orthogeomys_grandis,orthogeomys_grandis~
pygeretmus_pumilio,pygeretmus_pumilio~
thallomys_loringi,thallomys_loringi~
thallomys_nigricauda,thallomys_nigricauda~
thallomys_paedulcus,thallomys_paedulcus~
diplothrix_legata,diplothrix_legata~
bandicota_savilei,bandicota_savilei~
nesokia_indica,nesokia_indica~
bandicota_bengalensis,bandicota_bengalensis~
bandicota_indica,bandicota_indica~
parotomys_brantsii,parotomys_brantsii~
parotomys_littledalei,parotomys_littledalei~
neacomys_minutus,neacomys_minutus~
dephomys_defua,dephomys_defua~
otomys_tropicalis,otomys_tropicalis~
myodes_rufocanus,myodes_rufocanus~
peromyscus_crinitus,peromyscus_crinitus~
cricetulus_migratorius,cricetulus_migratorius~
cricetulus_migratorius_phaeus,cricetulus_migratorius_phaeus~
dipodomys_insularis,dipodomys_insularis~
dipodomys_simulans_peninsularis,dipodomys_simulans_peninsularis~
neotoma_isthmica,neotoma_isthmica~
oryzomys_mexicanus,oryzomys_mexicanus~
cratogeomys_fumosus_imparilis,cratogeomys_fumosus_imparilis~
cratogeomys_fumosus_angustirostris,cratogeomys_fumosus_angustirostris~
neotoma_lepida,neotoma_lepida~
peromyscus_aztecus_oaxacensis,peromyscus_aztecus_oaxacensis~
EOF

mkdir screen_logs
python ../../../../../../src/fix_species_names.py _species_name_change.txt RAxML_bestTree.BS_ML_GTRCAT rodentia_therest.aln ../names.txt ../lineage.txt > screen_logs/01_log1_fixspecies.txt

##>> NOTE: I needed to update the python script as it was not working. I just modifed `fix_alignment()` 
##>>       so the output file was different from the one being read

## Move unncessary output files to `00_filt1/` and rename output alignment
mkdir 00_filt1
mv _species_name_change.txt lineage.txt lineage.txt.1.bak 00_filt1/
rm rodentia_therest.aln
mv rodentia_therest_out.phylip alignment.phylip

# 2. Remove species 

# 2.1. Create input file for python script `prune_alignment.py`. This text file 
#      contains a list of the species that need to be removed 
cat > _species_remove.txt << EOF
rattus_leucopus_leucopus
peromyscus_difficilis_difficilis
oryzomys_couesi_couesi
cratogeomys_fumosus_fumosus
sicista_subtilis
phyllotis_osilae_osilae
cricetulus_migratorius_migratorius
phyllotis_xanthopygus_xanthopygus
abrothrix_olivaceus_olivaceus
apodemus_fulvipectus
rattus_tanezumi_sladeni
microtus_sikimensis
peromyscus_sagax
oryzomys_palustris_palustris
akodon_subfuscus_subfuscus
akodon_azarae_azarae
alticola_strelzowi_strelzowi
microtus_fortis
nyctomys_sumichrasti_sumichrasti
neotoma_albigula_albigula
peromyscus_pectoralis_pectoralis
peromyscus_truei_truei
peromyscus_gratus_gratus
peromyscus_nasutus_nasutus
peromyscus_boylii_boylii
peromyscus_levipes_levipes
peromyscus_polionotus_polionotus
napaeozapus_insignis_insignis
zapus_hudsonius_hudsonius
zapus_trinotatus_trinotatus
zapus_princeps_princeps
cratogeomys_castanops_castanops
cratogeomys_goldmani_goldmani
geomys_pinetis_pinetis
geomys_breviceps_breviceps
castor_fiber_fiber
mus_musculus_musculus
otomys_denti_denti
otomys_typus_typus
geomys_bursarius_bursarius
rheomys_thomasi
peromyscus_polius
handleyomys_intectus
wilfredomys_pictipes
microakodontomys_transitorius
podoxymys_roraimae
acomys_subspinosus
myopus_schisticolor
myomyscus_yemeni
malacothrix_typica
peromyscus_boylii
phyllotis_osilae_phaeus
arvicanthis_ansorgei
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
   in RStudio and change line 24 so it is `subt    <- "rodentia_therest"` and uncomment line 25.
   Now, we can run it from RStudio. This script will generate a concatenated alignment file with all partitions, as well as 
   one alignment file for each individual partitions,
   inside [`01_alignments`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/01_alignments).
   Log files and Rdata can be found [here](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/01_alignments/Rout/log_concatenation). 
      > NOTE: Updated `partitions.txt` file inside `00_mammal_alns/rodentia_therest/` generated.   
   * Paths have been automatically set according to current file architecture in `mammals` dir, 
   do not change paths in the Rscript! 

The final alignments generated at the end of this step can be downloaded from 
[here](https://www.dropbox.com/s/bo4shh56jrlswwi/SeqBayesS2_Raln_rodentia_therest.zip?dl=0). 
Nevertheless, these are not the final alignments as this data set was further subset 
into two data subsets. Please continue reading to learn more about this and how to 
obtain the final alignments!

# EXTRA FILTERING -- DATA SUBSETTING
When we first carried out the analysis with the tree topology as in step 1, we realised that 
it was too big to analyse with `MCMCtree` (i.e., too many taxa). Therefore, we decided to 
further subset this data and generate two data subsets: "rodentia subtree 1" and 
"rodentia subtree 2". 

The procedure followed was the following:

## 1. Explore partitioning the data set
The directory
[`00_R_parsing`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_therest/filter_aln/extra_filtering/00_R_parsing)
has the input/output files used to explore how to partition the data set into two data subsets. Please access this directory using the link 
provided to go through the steps followed. More details about this step can be found 
[here](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_therest/filter_aln/extra_filtering),
in the first section `1. Obtain subtrees`.

## 2. Generate alignments 
[Here](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_therest/filter_aln/extra_filtering), 
in section `2. Generating alignments`, 
you will read about all the steps we followed to extract the sequences of the taxa that need to be allocated to each data 
subset. The data can be downloaded 
[here](https://www.dropbox.com/s/h99ciuwwrvcbfep/SeqBayesS2_filteraln2_rodtherest_01_perl_parsing.zip?dl=0)
if you want to check you have reproduced our results, which should be saved 
in the directory
[`01_perl_parsing`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_therest/filter_aln/extra_filtering/01_perl_parsing).
In this GitHub repository, due to limited space for big files, you will see only the Perl script we use 
to generate the alignments.

## 3. Add new taxa and generate alignments
In order to avoid issues when grafting the subtrees to the backbone tree, we decided to add 
extra taxa to each subtree (two taxa from data subset 1 are included in data subset 2, and 
viceversa). The steps followed to add taxa to the first subtree can be 
found
[here](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_therest/filter_aln/extra_filtering/02_MAFFT_subt1),
while those followed to generate the one for the second subtree can be found 
[here](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_therest/filter_aln/extra_filtering/02_MAFFT_subt2).

## 4. Final alignments 
[Here](https://www.dropbox.com/s/5cxvn2fvqdevti8/SeqBayesS2_Raln_rod_subt1.zip?dl=0)
you can download the final alignments for subtree 1 
and
[here](https://www.dropbox.com/s/5cxvn2fvqdevti8/SeqBayesS2_Raln_rod_subt1.zip?dl=0)
for subtree 2. 