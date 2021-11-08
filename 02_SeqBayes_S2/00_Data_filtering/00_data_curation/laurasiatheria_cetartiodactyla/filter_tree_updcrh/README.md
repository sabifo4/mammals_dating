# Laurasiatheria cetartiodactyla - phylogeny

## 1. Get tree topology and add calibrations

### Setting maximum age for nodes that only had a minimum assigned following the fossil record
We note that there are several nodes that need to have the maximum age established. For that purpose, 
we found the ancestor of this node that was calibrated with an ST distribution and computed the 
2.5% quantile as it follows:

```
# First column contains all the nodes that do require a maximum age to be established 
# to use a soft-bound calibration
# The second column lists the ancestor of each of the nodes.
# The third column includes the ST distribution with which the ancestral nodes have 
# been calibrated -- this ST calibration is used to set the maximum age, see below.

SUINA|'B(0.347,x)'                     | ARTIOFABULA   | ST(0.553999686032752,0.00768136937401687,-0.964454809591589,7.34267207273677)

WHIPPOMORPHA|'B(0.524,x)'              | CETRUMINANTIA | ST(0.506969434346495,0.00916624957035315,-1.51323271125153,6.29387336904601)
MYSTICETI|'B(0.1599,x)'                | CETRUMINANTIA | ST(0.506969434346495,0.00916624957035315,-1.51323271125153,6.29387336904601)
ODONTOCETI|'B(0.2304,x)'               | CETRUMINANTIA | ST(0.506969434346495,0.00916624957035315,-1.51323271125153,6.29387336904601)
DELPHINIDA|'B(0.1599,x)'               | CETRUMINANTIA | ST(0.506969434346495,0.00916624957035315,-1.51323271125153,6.29387336904601)
PHOCOENIDAE-MONODONTIDAE|'B(0.076,x)'  | CETRUMINANTIA | ST(0.506969434346495,0.00916624957035315,-1.51323271125153,6.29387336904601)
HIPPOPOTAMIDAE|'B(0.0774,x)'           | CETRUMINANTIA | ST(0.506969434346495,0.00916624957035315,-1.51323271125153,6.29387336904601)
PHYSETEROIDEA|'B(0.1382,x)'            | CETRUMINANTIA | ST(0.506969434346495,0.00916624957035315,-1.51323271125153,6.29387336904601)
BALAENOPTERIDAE|'B(0.073,x)'           | CETRUMINANTIA | ST(0.506969434346495,0.00916624957035315,-1.51323271125153,6.29387336904601)
 
GIRAFFIDAE|'B(0.14,x)'                 | CETRUMINANTIA | ST(0.506969434346495,0.00916624957035315,-1.51323271125153,6.29387336904601)
STEM-MOSCHIDAE|'B(0.195,x)'            | CETRUMINANTIA | ST(0.506969434346495,0.00916624957035315,-1.51323271125153,6.29387336904601)

BOVINI|'B(0.102,x)'                    | BOVIDAE       | ST(0.16170125373343,0.0311120061766834,5.86301508375978,624.148265147333)
TRAGELAPHINI|'B(0.0549,x)'             | BOVIDAE       | ST(0.16170125373343,0.0311120061766834,5.86301508375978,624.148265147333)
REDUNCINI|'B(0.05111,x)'               | BOVIDAE       | ST(0.16170125373343,0.0311120061766834,5.86301508375978,624.148265147333)
HIPPOTRAGINI-ALCELAPHINI|'B(0.0648,x)' | BOVIDAE       | ST(0.16170125373343,0.0311120061766834,5.86301508375978,624.148265147333)
ALCELAPHINI|'B(0.0505,x)'              | BOVIDAE       | ST(0.16170125373343,0.0311120061766834,5.86301508375978,624.148265147333)
CAPRINAE|'B(0.089,x)'                  | BOVIDAE       | ST(0.16170125373343,0.0311120061766834,5.86301508375978,624.148265147333)

```

We can calculate the 2.5% quantile for each ST calibration as it follows: 

```R
# Artiofabula
sn::qst(0.025, 0.553999686032752,0.00768136937401687,-0.964454809591589,7.34267207273677)
#[1] 0.5325987

# Cetruminantia
sn::qst(0.025, 0.506969434346495,0.00916624957035315,-1.51323271125153,6.29387336904601)
#[1] 0.48019

# Bovidae
sn::qst(0.025, 0.16170125373343,0.0311120061766834,5.86301508375978,624.148265147333)
#[1] 0.1588611
```

In the end, the soft-bound calibrations for the nodes listed above are the following:

```
SUINA|'B(0.347,0.53)'                    

WHIPPOMORPHA|'B(0.524,0.48)'         
MYSTICETI|'B(0.1599,0.48)'               
ODONTOCETI|'B(0.2304,0.48)'              
DELPHINIDA|'B(0.1599,0.48)'              
PHOCOENIDAE-MONODONTIDAE|'B(0.076,0.48)' 
HIPPOPOTAMIDAE|'B(0.0774,0.48)'          
PHYSETEROIDEA|'B(0.1382,0.48)'           
BALAENOPTERIDAE|'B(0.073,0.48)'          
 
GIRAFFIDAE|'B(0.14,0.48)'                
STEM-MOSCHIDAE|'B(0.195,0.48)'           

BOVINI|'B(0.102,0.16)'                   
TRAGELAPHINI|'B(0.0549,0.16)'            
REDUNCINI|'B(0.05111,0.16)'              
HIPPOTRAGINI-ALCELAPHINI|'B(0.0648,0.16)'
ALCELAPHINI|'B(0.0505,0.16)'             
CAPRINAE|'B(0.089,0.16)'     
```

### Generating calibrated trees
For this subtree, we found that there were some conflict with one calibrated node: Whippomorpha. 
Therefore, we decided to estimate the timetrees by calibrating this node first with a soft bound 
and then with a low bound: `B(0.478,0.507)` and `L(0.507,0,0.01)`.

The R scripts, tree files, and text files used to generate the calibrated phylogeny can be 
found in this directory. Those files which file name includes the tag `v2` are used to generate 
the calibrated tree with Whippomorpha being calibrated with a low bound (`L(0.507,0,0.01)`), while 
the rest of files are used to obtain the calibrated tree in which Whippomorpha is calibrated with 
a soft bound calibration (`B(0.478,0.507)`).

Note that the R scripts can generates dummy alignments that can be used 
when running `MCMCtree` without the data to reduce disk space (see section 3). 
This "dummy" alignment is saved [here](../../../01_alignments/01_mammal_dummy_alns/laurasiatheria_cetartiodactyla).

After running this script, you will have the following files:

```
00_Filter_trees 
     |- RAxML_tree
     |         |- laurasiatheria_cetartiodactyla.tree             # File not used. Best-scoring ML tree obtained with RAxML
     |         
     |- 431sp_laurasiatheria_cetartiodactyla_MCMCtree_calib.tree  # File output by the R script - Whippomorpha: `B(0.478,0.507)`
     |- 431sp_laurasiatheria_cetartiodactyla_spnameslist.txt      # File output by the R script
     |- 431sp_lcetartiodactyla_MCMCtree_calib_v2.tree             # File output by the R script - Whippomorpha: `L(0.507,0,0.01)`
     |- Calibrations_LaurCetart.txt                               # Input file used by the R script (soft-bound). It matches the tag names
     |                                                            # in input tree with corresponding calibrations to be replaced
     |- Calibrations_LaurCetart_v2.txt                            # Input file used by the R script (low-bound). It matches the tag names
     |                                                            # in input tree with corresponding calibrations to be replaced
     |- Calibrations_Lcetartiodactyla.R                           # R script - Whippomorpha: `B(0.478,0.507)`
     |- Calibrations_Lcetartiodactyla_v2.R                        # R script - Whippomorpha: `L(0.507,0,0.01)`
     |- laurasiatheria_cetartiodactyla_rooted_calibnames.tree     # Input file used by the R script
     |- laurasiatheria_cetartiodactyla_rooted_baseml.tree         # File manually generated after running R script 
                                                                  # to be used by BASEML (calibrations manually removed)
```

Note that we have manually generated the
[`laurasiatheria_cetartiodactyla_rooted_baseml.tree`](00_Filter_trees/laurasiatheria_cetartiodactyla_rooted_baseml.tree),
which does 
not contain the calibrations. This file will be used when running `BASEML` to compute 
the Hessian and the gradient that are needed by `MCMCtree` to run the approximate 
likelihood.

## 2. Check if calibrations are in conflict
The trees described above were used to find 
if there were any conflicts with the calibations used.
You can download the directories 
with the results obtained when running `MCMCtree` without the data
[here](https://www.dropbox.com/s/ngwgzfm0thsqsim/SeqBayesS2_check_conflict_artiodactyla.zip?dl=0) 
for the tree in which a soft-bound calibration is used for Whippomorpha and 
[here](https://www.dropbox.com/s/xr68808m1vqnmnz/SeqBayesS2_check_conflict_artiodactyla_v2.zip?dl=0)
for when a low-bound calibration is used instead for this same node.

Once you download them, you should unzip its content and save the 
two directories inside 
[`01_Check_conflict`](01_Check_conflict) or 
[`01_Check_conflict_v2`](01_Check_conflict_v2), respectively for each subtree with different 
calibrations for node Whippomorpha. The file architecture is the following:

```
01_Check_conflict*                       # Either `01_Check_conflict` or `01_Check_conflict_v2` 
      |- 0X_Prior_*                      # Provided in the zip file, not in this repository due to lack of space
      |- outRdata                        # Directory  with R objects generated by the R script
      |- 00_Check_STanalitycalVSprior.R  # R script 
      |- *[pdf|png]                      # Files output by the R script, provided here
      |- *tsv                            # File output by the R script, provided here
```

Please read all the comments and explanations in
[the R script provided in this directory](01_Check_conflict/00_Check_STanalitycalVSprior.R) 
to understand each step that we followed to avoid having conflicting calibrations in
the tree topology. Sometimes, we might need to adjust the ST calibrations and/or maximum
bounds if the neighbouring calibrations are in conflict (e.g., there are truncation issues). 

In a nutshell:   

   1. First, we run `MCMCtree` without using the data (i.e., 
   without using the alignment, hence the "dummy" alignment used here) and fixing the
   tree topology where only the skew-_t_ (ST) calibrations have been added.   
   2. For each calibrated node, we plot the corresponding analytical ST distribution
   (the one that we have told `MCMCtree` to use) against the corresponding posterior density
   inferred by `MCMCtree` when no data are used (data described in step 1). In addition,
   we add to this plot the posterior density of this node that was inferred by `MCMCtree`
   when using the first data set (72-taxon data set).   
   3. To check for conflict, we do the following for each calibrated node with an 
   ST calibration:   
      * Estimate mean times and quantiles (2.5% and 97.5%) from the posterior density
	  inferred by `MCMCtree` when the data are not used and the fixed tree topology has only
	  ST distributions.   
	  * Estimate mean times and quantiles from the posterior density inferred with
	  data set 1 (72-taxon data set) for the same node.   
	  * Check how much the former deviate from the latter.   
	  * If deviation is lower than ~5%, proceed with step 4.   
   4. If checks in step 3 are ok, we run `MCMCtree` without the data alignment but
   the tree topology now has both the ST calibrations and the calibrations with soft
   bounds (i.e., calibrations that have a minimum and a maximum bound with a 2.5% tail
   probability in each side).   
   5. Then, we generate the same plot as described in step 2.    
   6. Last, we check again for possible conflict as described in step 3. If deviation
   is lower than ~5% for all calibrated nodes, this is the end of the checks. Otherwise, we need 
   to adjust the location and scale parameters of the ST calibrations until no conflict
   is observed by subtracting the corresponding deviation (more details in the R script
   if this adjustment is taking place).   

In this case (laurasiatheria_cetartiodactyla), we did not have to adjust any calibration as there were no
conflicts encountered (see plot below but also other plots within this directory):

## Artiodactyla -- subtree with Whippomorpha calibrated with a soft-bound calibration

**When using only ST calibrations**   
Calibrations used:   
   * Mammalia: ST(1.642,0.425,12.652,1714.565)   
   * Artiodactyla: ST(0.577,0.007,-0.634,7.509)   
   * Artiofabula: ST(0.554,0.008,-0.964,7.343)   
   * Cetruminantia: ST(0.507,0.009,-1.513,6.294)   
   * Bovidae: ST(0.1617,0.0311,5.8630,624.1483)    
   * *Ovis aries*-*Capra hircus*: ST(0.0686,0.0172,3.2247,169.6716 )   

<p align="center">
  <img width="1000" height="600" src="01_Check_conflict/00_Only_ST_L.cetartiodactyla_MCMCruns.png">
</p>

**When using both ST and soft bound calibrations**   
Calibrations used:   
   * Mammalia: ST(1.642,0.425,12.652,1714.565)   
   * Artiodactyla: ST(0.577,0.007,-0.634,7.509)   
   * Artiofabula: ST(0.554,0.008,-0.964,7.343)   
   * Cetruminantia: ST(0.507,0.009,-1.513,6.294)   
   * **Bovidae: ~ST(0.1617,0.0311,5.8630,624.1483)~ -- CONFLICT, adjusted**   
   **to ST(0.1395,0.0246,5.8630,624.1483) for next round**    
   * **_Ovis aries_-_Capra hircus_: ~ST(0.0686,0.0172,3.2247,169.6716)~ -- CONFLICT, adjusted**   
   **to ST(0.0594,0.0157,3.2247,169.6716) for next round**    
   * Monotremata: B(0.2446,1.332)    
   * Tachyglossidae: B(0.0258,1.332)   
   * Suina: B(0.347,0.53)   
   * Whippomorpha: B(0.478,0.507)    
   * Cetacea: B(0.3613,0.56)   
   * Mysticeti: B(0.1599,0.48)   
   * Odontoceti: B(0.2304,0.48)   
   * Delphinida: B(0.1599,0.48)   
   * Phocoenidae-Monodontidae: B(0.076,0.48)   
   * Hippopotamidae: B(0.0774,0.48)    
   * Giraffidae: B(0.14,0.48)   
   * Bovini: B(0.102,0.16)   
   * Tragelaphini: B(0.0549,0.16)   
   * Reduncini: B(0.05111,0.16)   
   * Hippotragini-Alcelaphini: B(0.0648,0.16)   
   * Alcelaphini: B(0.0505,0.16)   
   * Caprinae: B(0.089,0.16)   
   * Cervidae: B(0.17235,0.284)   
   * Stem-Moschidae: B(0.195,0.48)   
   * Neobalaeninae: B(0.2304,0.339)   
   * Balaenopteridae: B(0.073,0.48)   
   * Physeteroidea: B(0.1382,0.48)   

<p align="center">
  <img width="1000" height="600" src="01_Check_conflict/01_SBnST_L.cetartiodactyla_MCMCruns.png">
</p>

**When using both ST and soft bound calibrations - 2nd round**   
Calibrations used:   
   * Mammalia: ST(1.642,0.425,12.652,1714.565)   
   * Artiodactyla: ST(0.577,0.007,-0.634,7.509)   
   * Artiofabula: ST(0.554,0.008,-0.964,7.343)   
   * Cetruminantia: ST(0.507,0.009,-1.513,6.294)   
   * **Bovidae: ~ST(0.1395,0.0246,5.8630,624.1483)~ -- CONFLICT, adjusted**   
   **to ST(0.1489,0.0245,5.8630,624.1483) for next round**    
   * _Ovis aries_-_Capra hircus_: ST(0.0594,0.0157,3.2247,169.6716)   
   * Monotremata: B(0.2446,1.332)    
   * Tachyglossidae: B(0.0258,1.332)   
   * Suina: B(0.347,0.53)   
   * Whippomorpha: B(0.478,0.507)    
   * Cetacea: B(0.3613,0.56)   
   * Mysticeti: B(0.1599,0.48)   
   * Odontoceti: B(0.2304,0.48)   
   * Delphinida: B(0.1599,0.48)   
   * Phocoenidae-Monodontidae: B(0.076,0.48)   
   * Hippopotamidae: B(0.0774,0.48)    
   * Giraffidae: B(0.14,0.48)   
   * Bovini: B(0.102,0.16)   
   * Tragelaphini: B(0.0549,0.16)   
   * Reduncini: B(0.05111,0.16)   
   * Hippotragini-Alcelaphini: B(0.0648,0.16)   
   * Alcelaphini: B(0.0505,0.16)   
   * Caprinae: B(0.089,0.16)   
   * Cervidae: B(0.17235,0.284)   
   * Stem-Moschidae: B(0.195,0.48)   
   * Neobalaeninae: B(0.2304,0.339)   
   * Balaenopteridae: B(0.073,0.48)   
   * Physeteroidea: B(0.1382,0.48)   
   
<p align="center">
  <img width="1000" height="600" src="01_Check_conflict/02_SBnST_tweak1_L.cetartiodactyla_MCMCruns.png">
</p>

**When using both ST and soft bound calibrations - 3rd round**   
Calibrations used:   
   * Mammalia: ST(1.642,0.425,12.652,1714.565)   
   * Artiodactyla: ST(0.577,0.007,-0.634,7.509)   
   * Artiofabula: ST(0.554,0.008,-0.964,7.343)   
   * Cetruminantia: ST(0.507,0.009,-1.513,6.294)   
   * Bovidae: ST(0.1489,0.0245,5.8630,624.1483)
   * _Ovis aries_-_Capra hircus_: ST(0.0594,0.0157,3.2247,169.6716)   
   * Monotremata: B(0.2446,1.332)    
   * Tachyglossidae: B(0.0258,1.332)   
   * Suina: B(0.347,0.53)   
   * Whippomorpha: B(0.478,0.507)    
   * Cetacea: B(0.3613,0.56)   
   * Mysticeti: B(0.1599,0.48)   
   * Odontoceti: B(0.2304,0.48)   
   * Delphinida: B(0.1599,0.48)   
   * Phocoenidae-Monodontidae: B(0.076,0.48)   
   * Hippopotamidae: B(0.0774,0.48)    
   * Giraffidae: B(0.14,0.48)   
   * Bovini: B(0.102,0.16)   
   * Tragelaphini: B(0.0549,0.16)   
   * Reduncini: B(0.05111,0.16)   
   * Hippotragini-Alcelaphini: B(0.0648,0.16)   
   * Alcelaphini: B(0.0505,0.16)   
   * Caprinae: B(0.089,0.16)   
   * Cervidae: B(0.17235,0.284)   
   * Stem-Moschidae: B(0.195,0.48)   
   * Neobalaeninae: B(0.2304,0.339)   
   * Balaenopteridae: B(0.073,0.48)   
   * Physeteroidea: B(0.1382,0.48)   
   
<p align="center">
  <img width="1000" height="600" src="01_Check_conflict/03_SBandSTtweak2_L.cetartiodactyla_MCMCruns.png">
</p>

**Deviations (main 72-taxa VS laurasiatheria_cetartiodactyla data sets)**   
<p align="center">
  <img width="1000" height="600" src="01_Check_conflict/03_SBnST_tweak2_L.cetartiodactyla_meanquant.png">
</p>

The final tree topology can be found in the
[`final_tree_topology`](02_Final_tree_topology)
directory.


## Artiodactyla -- subtree with Whippomorpha calibrated with a low-bound calibration

**When using only ST calibrations**   
Calibrations used:   
   * Mammalia: ST(1.642,0.425,12.652,1714.565)   
   * Artiodactyla: ST(0.577,0.007,-0.634,7.509)   
   * Artiofabula: ST(0.554,0.008,-0.964,7.343)   
   * Cetruminantia: ST(0.507,0.009,-1.513,6.294)   
   * Bovidae: ST(0.1617,0.0311,5.8630,624.1483)    
   * *Ovis aries*-*Capra hircus*: ST(0.0686,0.0172,3.2247,169.6716 )   

<p align="center">
  <img width="1000" height="600" src="01_Check_conflict_v2/00_Only_ST_L.cetartiodactyla_MCMCruns.png">
</p>

**When using both ST and soft bound calibrations**   
Calibrations used:   
   * Mammalia: ST(1.642,0.425,12.652,1714.565)   
   * Artiodactyla: ST(0.577,0.007,-0.634,7.509)   
   * Artiofabula: ST(0.554,0.008,-0.964,7.343)   
   * Cetruminantia: ST(0.507,0.009,-1.513,6.294)~ -- CONFLICT, adjusted**   
   **to ST(0.494,0.014,-1.513,6.294) for next round**       
   * **Bovidae: ~ST(0.1617,0.0311,5.8630,624.1483)~ -- CONFLICT, adjusted**   
   **to ST(0.140,0.0247,5.8630,624.1483) for next round**    
   * **_Ovis aries_-_Capra hircus_: ~ST(0.0686,0.0172,3.2247,169.6716)~ -- CONFLICT, adjusted**   
   **to ST(0.0596,0.0159,3.2247,169.6716) for next round**    
   * Monotremata: B(0.2446,1.332)    
   * Tachyglossidae: B(0.0258,1.332)   
   * Suina: B(0.347,0.53)   
   * Whippomorpha: L(0.507,0,0.01)    
   * Cetacea: B(0.3613,0.56)   
   * Mysticeti: B(0.1599,0.48)   
   * Odontoceti: B(0.2304,0.48)   
   * Delphinida: B(0.1599,0.48)   
   * Phocoenidae-Monodontidae: B(0.076,0.48)   
   * Hippopotamidae: B(0.0774,0.48)    
   * Giraffidae: B(0.14,0.48)   
   * Bovini: B(0.102,0.16)   
   * Tragelaphini: B(0.0549,0.16)   
   * Reduncini: B(0.05111,0.16)   
   * Hippotragini-Alcelaphini: B(0.0648,0.16)   
   * Alcelaphini: B(0.0505,0.16)   
   * Caprinae: B(0.089,0.16)   
   * Cervidae: B(0.17235,0.284)   
   * Stem-Moschidae: B(0.195,0.48)   
   * Neobalaeninae: B(0.2304,0.339)   
   * Balaenopteridae: B(0.073,0.48)   
   * Physeteroidea: B(0.1382,0.48)   

<p align="center">
  <img width="1000" height="600" src="01_Check_conflict_v2/01_SBnST_L.cetartiodactyla_MCMCruns.png">
</p>

**When using both ST and soft bound calibrations - 2nd round**   
Calibrations used:   
   * Mammalia: ST(1.642,0.425,12.652,1714.565)   
   * Artiodactyla: ST(0.577,0.007,-0.634,7.509)   
   * Artiofabula: ST(0.554,0.008,-0.964,7.343)   
   * Cetruminantia: ST(0.494,0.014,-1.513,6.294)~ -- CONFLICT, adjusted**   
   **to ST(0.479,0.016,-1.513,6.294) for next round**       
   * **Bovidae: ~ST(0.140,0.0247,5.8630,624.1483)~ -- CONFLICT, adjusted**   
   **to ST(0.149,0.0247,5.8630,624.1483) for next round**    
   * **_Ovis aries_-_Capra hircus_: ~ST(0.0596,0.0159,3.2247,169.6716)~ -- CONFLICT, adjusted**   
   **to ST(0.0582,0.014,3.2247,169.6716) for next round**    
   * Monotremata: B(0.2446,1.332)    
   * Tachyglossidae: B(0.0258,1.332)   
   * Suina: B(0.347,0.53)   
   * Whippomorpha: B(0.478, 0.507)    
   * Cetacea: B(0.3613,0.56)   
   * Mysticeti: B(0.1599,0.48)   
   * Odontoceti: B(0.2304,0.48)   
   * Delphinida: B(0.1599,0.48)   
   * Phocoenidae-Monodontidae: B(0.076,0.48)   
   * Hippopotamidae: B(0.0774,0.48)    
   * Giraffidae: B(0.14,0.48)   
   * Bovini: B(0.102,0.16)   
   * Tragelaphini: B(0.0549,0.16)   
   * Reduncini: B(0.05111,0.16)   
   * Hippotragini-Alcelaphini: B(0.0648,0.16)   
   * Alcelaphini: B(0.0505,0.16)   
   * Caprinae: B(0.089,0.16)   
   * Cervidae: B(0.17235,0.284)   
   * Stem-Moschidae: B(0.195,0.48)   
   * Neobalaeninae: B(0.2304,0.339)   
   * Balaenopteridae: B(0.073,0.48)   
   * Physeteroidea: B(0.1382,0.48)   
   
<p align="center">
  <img width="1000" height="600" src="01_Check_conflict_v2/02_SBnST_tweak1_L.cetartiodactyla_MCMCruns.png">
</p>

**When using both ST and soft bound calibrations - 3rd round**   
Calibrations used:   
Calibrations used:   
   * Mammalia: ST(1.642,0.425,12.652,1714.565)   
   * Artiodactyla: ST(0.577,0.007,-0.634,7.509)   
   * Artiofabula: ST(0.554,0.008,-0.964,7.343)   
   * Cetruminantia: ST(0.479,0.016,-1.513,6.294)       
   * Bovidae: ST(0.149,0.0247,5.8630,624.1483)    
   * _Ovis aries_-_Capra hircus_: ST(0.0582,0.014,3.2247,169.6716)    
   * Monotremata: B(0.2446,1.332)    
   * Tachyglossidae: B(0.0258,1.332)   
   * Suina: B(0.347,0.53)   
   * Whippomorpha: B(0.478, 0.507)    
   * Cetacea: B(0.3613,0.56)   
   * Mysticeti: B(0.1599,0.48)   
   * Odontoceti: B(0.2304,0.48)   
   * Delphinida: B(0.1599,0.48)   
   * Phocoenidae-Monodontidae: B(0.076,0.48)   
   * Hippopotamidae: B(0.0774,0.48)    
   * Giraffidae: B(0.14,0.48)   
   * Bovini: B(0.102,0.16)   
   * Tragelaphini: B(0.0549,0.16)   
   * Reduncini: B(0.05111,0.16)   
   * Hippotragini-Alcelaphini: B(0.0648,0.16)   
   * Alcelaphini: B(0.0505,0.16)   
   * Caprinae: B(0.089,0.16)   
   * Cervidae: B(0.17235,0.284)   
   * Stem-Moschidae: B(0.195,0.48)   
   * Neobalaeninae: B(0.2304,0.339)   
   * Balaenopteridae: B(0.073,0.48)   
   * Physeteroidea: B(0.1382,0.48)   
   
<p align="center">
  <img width="1000" height="600" src="01_Check_conflict_v2/03_SBandSTtweak2_L.cetartiodactyla_MCMCruns.png">
</p>

**NOTE:** We decided to stop at this stage the manual adjustments as we decided that 
the timetree inferred with the calibrated tree topology with the soft-bound calibration
for node Whippomorpha was going to be stitched onto the backbone tree at the end 
of this sequential Bayesian dating approach. This step might take one to two days to 
finish so we decided to stop at this stage - the bias for mean time estimate for the 
calibrated nodes was lower than 5%, only the bias for the 97.5% quantile did not pass 
this threshold filter (i.e., -8.1%). We believe that the conflict is not large enough so 
we proceed with timetree inference with the prior calibrations detailed above.

**Deviations (main 72-taxa VS laurasiatheria_cetartiodactyla data sets)**   
<p align="center">
  <img width="1000" height="600" src="01_Check_conflict_v2/03_SBnST_tweak2_L.cetartiodactyla_meanquant.png">
</p>

The final tree topology can be found in the
[`final_tree_topology`](02_Final_tree_topology)
directory.

--- 

The next step is to run `MCMCtree` with the final tree topology and the 5-partitions 
alignment! Before that, however, we need to run `BASEML` to calculate the Hessian and 
the gradient, which are needed for the approximate likelihood calculation used by 
`MCMCtree` to speed up the Bayesian inference of divergence times.