# Laurasiatheria the rest - phylogeny

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

MUSTELIDAE-PROCYONIDAE|'B(0.2642,x)'  | ARCTOIDEA | ST(0.397114575053209,0.0313051999152975,-1.41518497370544,427.391036840726)
LOBODONTINI|'B(0.0505,x)'             | ARCTOIDEA | ST(0.397114575053209,0.0313051999152975,-1.41518497370544,427.391036840726)
PHOCIDAE|'B(0.1382,x)'                | ARCTOIDEA | ST(0.397114575053209,0.0313051999152975,-1.41518497370544,427.391036840726)
OTARIOIDEA|'B(0.1597,x)'              | ARCTOIDEA | ST(0.397114575053209,0.0313051999152975,-1.41518497370544,427.391036840726)

FELIFORMIA|'B(0.19535,x)'             | CARNIVORA | ST(0.53833781677108,0.0270936425266216,-2.11521868565242,70.3753474590593)
VIVERRINAE-GENETTINAE|'B(0.2044,x)'   | CARNIVORA | ST(0.53833781677108,0.0270936425266216,-2.11521868565242,70.3753474590593)
```

We can calculate the 2.5% quantile for each ST calibration as it follows: 

```
# Arctoidea
sn::qst(0.025, 0.397114575053209,0.0313051999152975,-1.41518497370544,427.391036840726)
#[1] 0.3267025

# Carnivora
sn::qst(0.025, 0.53833781677108,0.0270936425266216,-2.11521868565242,70.3753474590593)
#[1] 0.4762834
```

In the end, the soft-bound calibrations for the nodes listed above are the following:

```
MUSTELIDAE-PROCYONIDAE|'B(0.2642,0.33)
LOBODONTINI|'B(0.0505,0.33)'          
PHOCIDAE|'B(0.1382,0.33)'             
OTARIOIDEA|'B(0.1597,0.33)'                       

FELIFORMIA|'B(0.19535,0.48)'          
VIVERRINAE-GENETTINAE|'B(0.2044,0.48)'
```

### Generating calibrated trees 
We use the R script [`Calibrations_Ltherest.R`](00_Filter_trees/Calibrations_Ltherest.R)
to generate the phylogeny for this data subset. Note that we use the
[`laurasiatheria_therest_rooted_calibnames.tree`](00_Filter_trees/laurasiatheria_therest_rooted_calibnames.tree)
file, where tag names have been manually added in the 
nodes that are to be calibrated. These tag names are later replaced with the
corresponding calibrations specified in the 
[`laurasiathera_therest_calibrations.txt`](00_Filter_trees/laurasiathera_therest_calibrations.txt)
file. 
In addition, this R script generates dummy alignments that can be used 
when running `MCMCtree` without the data to reduce disk space (see next section 3). 
This "dummy" alignment is saved [here](../../../01_alignments/01_mammal_dummy_alns/laurasiatheria_therest)

After running this script, you will have the following files:

```
00_Filter_trees
        |- RAxML_tree
        |         |- laurasiatheria_therest.tree              # File not used. Best-scoring ML tree obtained with RAxML
        |         
        |- 659sp_laurasiatheria_therest_MCMCtree_calib.tree   # File output by the R script
        |- 659sp_laurasiatheria_therest_spnameslist.txt       # File output by the R script
        |- Calibrations_Ltherest.R                            # R script
        |- laurasiathera_therest_calibrations.txt             # Input file used by the R script. It matches the tag names
        |                                                     # in input tree with corresponding calibrations to be replaced
        |- laurasiatheria_therest_rooted_baseml.tree          # File manually generated after running R script 
        |                                                     # to be used by BASEML (calibrations manually removed)
        |- laurasiatheria_therest_rooted_calibnames.tree      # Input file used by the R script
```

Note that we have manually generated the
[`laurasiatheria_therest_rooted_baseml.tree`](00_Filter_trees/laurasiatheria_therest_rooted_baseml.tree),
which does 
not contain the calibrations. This file was used when running `BASEML` to compute 
the Hessian and the gradient that are needed by `MCMCtree` to run the approximate 
likelihood before we had to add new taxa to the alignment (see below).

## 2. Tree topology change 
After an extra data filtering when we added four extra taxa (i.e., *Pteropus vampyrus*,
*Myotis lucifugus*, *Vicugna pacos*, and *Capra hircus*), the tree topology changed to include the placement of these 
taxa (see the details in
[this `README.md` file](../filter_aln/README.md),
section `EXTRA FILTERING -- ADDING TAXA TO THE ALIGNMENT`, if you did not go through the data filtering before,
which explains why we added these four taxa 
and how this was done).

## 3. Check if calibrations are in conflict
The tree described above (659 taxa with *Pteropus vampyrus*,
*Myotis lucifugus*, *Vicugna pacos*, and *Capra hircus*) was used to find 
if there were any conflicts with the calibations used.
You can download the directories 
with the results obtained when running `MCMCtree` without the data
[here](https://www.dropbox.com/s/frurexfq2ntyxu9/SeqBayesS2_check_conflict_laurasiatheria_therest.zip?dl=0).
Once you download them, you should unzip its content and save them 
inside the 
[`01_Check_conflict`](01_Check_conflict)
directory so the file architecture is the following:

```
01_Check_conflict 
      |- 00_Prior_onlyST                 # Provided in the zip file, not in this repository due to lack of space
      |- 01_Prior_SBandST                # Provided in the zip file, not in this repository due to lack of space
      |- 02_Prior_SBandST_tweak1         # Provided in the zip file, not in this repository due to lack of space
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

In this case (laurasiatheria_therest), we did not have to adjust any calibration as there were no
conflicts encountered (see plot below but also other plots within this directory):

**When using only ST calibrations**   
Calibrations used:   
   * Mammalia: ST(1.642,0.425,12.652,1714.565)    
   * Laurasiatheria: ST(0.694,0.006,0.431,4.953)   
   * Erinaceidae-Soricidae: ST(0.596,0.017,-1.243,19.856)   
   * Scrotifera: ST(0.678,0.006,0.37,4.464)   
   * Chiroptera: ST(0.596,0.016,-1.239,13.572)   
   * Fereungulata: ST(0.671,0.005,0.355,4.329)   
   * Carnivora: ST(0.538,0.027,-2.115,70.375)   
   * Felidae: ST(0.1482,0.0212,0.4170,367.2184)   
   * Pantherinae: ST(0.0699,0.0152,1.1198,331.4972)   
   * Caniformia: ST(0.467,0.031,-1.966,117.012)    
   * Arctoidea: ST(0.397,0.031,-1.415,427.391)   
   * Euungulata: ST(0.660,0.005,0.341,4.207)   
   * Artiodactyla: ST(0.577,0.007,-0.634,7.509)   
   
<p align="center">
  <img width="1000" height="600" src="01_Check_conflict/00_Only_ST_L.therest_MCMCruns.png">
</p>

**When using both ST and soft bound calibrations**   
Calibrations used:   
   * Mammalia: ST(1.642,0.425,12.652,1714.565)    
   * Laurasiatheria: ST(0.694,0.006,0.431,4.953)   
   * Erinaceidae-Soricidae: ST(0.596,0.017,-1.243,19.856)   
   * Scrotifera: ST(0.678,0.006,0.37,4.464)   
   * Chiroptera: ST(0.596,0.016,-1.239,13.572)   
   * Fereungulata: ST(0.671,0.005,0.355,4.329)   
   * Carnivora: ST(0.538,0.027,-2.115,70.375)   
   * **Felidae: ~ST(0.1482,0.0212,0.4170,367.2184)~ -- CONFLICT, adjusted**   
   **to ST(0.14,0.02,0.4170,367.2184) for next round**    
   * **Pantherinae: ~ST(0.0699,0.0152,1.1198,331.4972)~ -- CONFLICT, adjusted**   
   **to ST(0.0666,0.0154,1.1198,331.4972) for next round**     
   * Caniformia: ST(0.467,0.031,-1.966,117.012)    
   * Arctoidea: ST(0.397,0.031,-1.415,427.391)   
   * Euungulata: ST(0.660,0.005,0.341,4.207)   
   * Artiodactyla: ST(0.577,0.007,-0.634,7.509)   
   * Monotremata: B(0.24459,1.332)   
   * Tachyglossidae: B(0.0258,1.332)   
   * Perissodactyla: B(0.555,0.6609)   
   * Ceratomorpha: B(0.48078,0.6609)   
   * Prionodontidae-Felidae: B(0.281,0.6609)   
   * Herpestidae-Eupleridae: B(0.1597,0.339)   
   * Mustelidae-Procyonidae: B(0.2642,0.33)   
   * Feliformia: B(0.19535,0.48)   
   * Viverrinae-Genettinae: B(0.2044,0.48)   
   * Lobodontini: B(0.0505,0.33)   
   * Phocidae: B(0.1382,0.33)   
   * Otarioidea: B(0.1597,0.33)   
   * Pinnipedia: B(0.2045,0.2729)   
   
<p align="center">
  <img width="1000" height="600" src="01_Check_conflict/01_SBnST_L.therest_MCMCruns.png">
</p>

**When using both ST and soft bound calibrations - 1st round**   
Calibrations used:   
   * Mammalia: ST(1.642,0.425,12.652,1714.565)    
   * Laurasiatheria: ST(0.694,0.006,0.431,4.953)   
   * Erinaceidae-Soricidae: ST(0.596,0.017,-1.243,19.856)   
   * Scrotifera: ST(0.678,0.006,0.37,4.464)   
   * Chiroptera: ST(0.596,0.016,-1.239,13.572)   
   * Fereungulata: ST(0.671,0.005,0.355,4.329)   
   * Carnivora: ST(0.538,0.027,-2.115,70.375)   
   * Felidae: ST(0.14,0.02,0.4170,367.2184)
   * Pantherinae: ST(0.0666,0.0154,1.1198,331.4972)
   * Caniformia: ST(0.467,0.031,-1.966,117.012)    
   * Arctoidea: ST(0.397,0.031,-1.415,427.391)   
   * Euungulata: ST(0.660,0.005,0.341,4.207)   
   * Artiodactyla: ST(0.577,0.007,-0.634,7.509)   
   * Monotremata: B(0.24459,1.332)   
   * Tachyglossidae: B(0.0258,1.332)   
   * Perissodactyla: B(0.555,0.6609)   
   * Ceratomorpha: B(0.48078,0.6609)   
   * Prionodontidae-Felidae: B(0.281,0.6609)   
   * Herpestidae-Eupleridae: B(0.1597,0.339)   
   * Mustelidae-Procyonidae: B(0.2642,0.33)   
   * Feliformia: B(0.19535,0.48)   
   * Viverrinae-Genettinae: B(0.2044,0.48)   
   * Lobodontini: B(0.0505,0.33)   
   * Phocidae: B(0.1382,0.33)   
   * Otarioidea: B(0.1597,0.33)   
   * Pinnipedia: B(0.2045,0.2729)   
   
<p align="center">
  <img width="1000" height="600" src="01_Check_conflict/02_SBandSTtweak1_L.therest_MCMCruns.png">
</p>

**Deviations (main 72-taxa VS laurasiatheria_therest data sets)**   
<p align="center">
  <img width="1000" height="600" src="01_Check_conflict/02_SBnST_tweak1_L.therest_meanquant.png">
</p>

The final tree topology can be found in the
[`final_tree_topology`](02_Final_tree_topology)
directory.

--- 

The next step is to run `MCMCtree` with the final tree topology and the 5-partitions 
alignment! Before that, however, we need to run `BASEML` to calculate the Hessian and 
the gradient, which are needed for the approximate likelihood calculation used by 
`MCMCtree` to speed up the Bayesian inference of divergence times.