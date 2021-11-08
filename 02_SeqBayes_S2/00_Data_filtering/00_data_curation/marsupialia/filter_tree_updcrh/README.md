# Marsupialia - phylogeny

## 1. Get tree topology and add calibrations

### Setting maximum age for nodes that only had a minimum assigned following the fossil record
We note that there are several nodes that need to have the maximum age established. For that purpose, 
we found the ancestor of this node, Eometatheria, that was calibrated with an ST distribution and computed the 
2.5% quantile as it follows: Eometatheria. 

```
# List of nodes that require a maximum age to be established 
# to use a soft-bound calibration
DASYUROMORPHIA|'B(0.1597,X)'
VOMBATIFORMES|'B(0.255,X)'
PHALANGERIDAE-BURRAMYIDAE|'B(0.25,X)'
PETAURIDAE-PSEUDOCHEIRIDAE|'B(0.255,X)'
MACROPODOIDEA|'B(0.247,X)'
```

We can calculate the 2.5% quantile for the ST calibration for Eometatheria
(`ST(0.45923291469171,0.0682847927346442,0.0248906633902748,623.433453465944)`)
as it follows: 

```
# Use `sn` package
sn::qst(0.025, 0.45923291469171,0.0682847927346442,0.0248906633902748,623.433453465944)
# 0.3265238
```

In the end, the soft-bound calibrations for the nodes listed above are the following:

```
DASYUROMORPHIA|'B(0.1597,0.33)'
VOMBATIFORMES|'B(0.255,0.33)'
PHALANGERIDAE-BURRAMYIDAE|'B(0.25,0.33)'
PETAURIDAE-PSEUDOCHEIRIDAE|'B(0.255,0.33)'
MACROPODIDAE-POTORIDAE|'B(0.247,0.33)'
```

### Generating calibrated trees 
We use the R script [`Calibrations_Marsupialia.R`](00_Filter_trees/Calibrations_Marsupialia.R)
to generate the phylogeny for this data subset. Note that we use the
[`marsupialia_rooted_calibnames.tree`](00_Filter_trees/marsupialia_rooted_calibnames.tree)
file, where tag names have been manually added in the 
nodes that are to be calibrated. These tag names are later replaced with the
corresponding calibrations specified in the 
[`Calibrations_Marsupialia.txt`](00_Filter_trees/Calibrations_Marsupialia.txt)
file. 
In addition, this R script generates dummy alignments that can be used 
when running `MCMCtree` without the data to reduce disk space (see next section 3). 
This "dummy" alignment is saved [here](../../../01_alignments/01_mammal_dummy_alns/marsupialia).

After running this script, you will have the following files:

```
00_Filter_trees 
     |- RAxML_tree
     |         |- marsupialia.tree             # File not used. Best-scoring ML tree obtained with RAxML
     |         
     |- 307sp_Marsupialia_MCMCtree_calib.tree  # File output by the R script
     |- 307sp_Marsupialia_spnameslist.txt      # File output by the R script
     |- Calibrations_Marsupialia.R             # R script
     |- Calibrations_Marsupialia.txt           # Input file used by the R script. It matches the tag names
     |                                         # in input tree with corresponding calibrations to be replaced
     |- marsupialia_rooted_baseml.tree         # File manually generated after running R script 
     |                                         # to be used by BASEML (calibrations manually removed)
     |- marsupialia_rooted_calibnames.tree     # Input file used by the R script
```

Note that we manually generated the
[`marsupialia_rooted_baseml.tree`](00_Filter_trees/marsupialia_rooted_baseml.tree),
which does not contain the calibrations. 

## 2. Manual change due to conflict with literature
The pyhlogeny in
[`307sp_Marsupialia_MCMCtree_calib.tree`](00_Filter_trees/307sp_Marsupialia_MCMCtree_calib.tree)
has a different tree topology if compared to the inferred best-scoring ML tree,
[`marsupialia.tree`](00_Filter_trees/RAxML_tree/marsupialia.tree). Specifically, we 
updated the Marsupialia tree topology according to the literature, which follows 
the next clustering:

```
# The topology from Mitchell et al. 2014 MB&E shows a pretty robust dataset
# & topology, with an Australidelphia like the tree shown below as 
# it supp fig. s1 [here](https://eur01.safelinks.protection.outlook.com/?url=https%3A%2F%2Facademic.oup.com%2Fmbe%2Farticle%2F31%2F9%2F2322%2F2925703&amp;data=02%7C01%7C%7C03e288fc3d954af8220408d80e43f743%7C569df091b01340e386eebd9cb9e25814%7C0%7C0%7C637275032363930685&amp;sdata=41BkZWuHS8z9jeHyoxXWRKnvxg1QTl8eSkhFFvwJOpc%3D&amp;reserved=0)
(Dromiciops,((notoryctids,(peramelians,dasyuromorphs)),diprotodonts))
```

## 3. Check if calibrations are in conflict
The tree with the topology described above was used to check 
if there were any conflicts with the calibations used.
You can download the directories 
with the results obtained when running `MCMCtree` without the data
[here](https://www.dropbox.com/s/n7zm9mt73jg2jwl/SeqBayesS2_check_conflict_marsupialia.zip?dl=0).
Once you download them, you should unzip its content and save the 
directories inside the 
[`01_Check_conflict`](01_Check_conflict)
directory so the file architecture is the following:

```
01_Check_conflict 
      |- 00_Prior_onlyST                   # Provided in the zip file, not in this repository due to lack of space
      |- 01_Prior_SBandST                  # Provided in the zip file, not in this repository due to lack of space
      |- outRdata                          # Directory  with R objects generated by the R script
      |- 00_Check_STanalitycalVSprior.R    # R script 
      |- *[pdf|png]                        # Files output by the R script, provided here
      |- *tsv                              # File output by the R script, provided here
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

In this case (marsupialia), there is no need for adjustments (see plots below but also
other plots within this directory):

**When using only ST calibrations**   
Calibrations used:   
   * Mammalia: ST(1.642,0.425,12.652,1714.565)  
   * Didelphimorphia-Australidelphia: ST(0.562,0.08,0.035,299.837)   
   * Eometatheria: ST(0.459,0.068,0.025,623.433)   
   
<p align="center">
  <img width="1000" height="600" src="01_Check_conflict/00_Only_ST_Marsupialia_MCMCruns.png">
</p>


**When using both ST and soft bound calibrations**   
Calibrations used:   
   * Mammalia: ST(1.642,0.425,12.652,1714.565)  
   * Didelphimorphia-Australidelphia: ST(0.562,0.08,0.035,299.837)   
   * Eometatheria: ST(0.459,0.068,0.025,623.433)   
   * Monotremata: B(0.24459,1.332)   
   * Tachyglossidae: B(0.0258,1.332)   
   * Marsupialia: B(0.476,1.313)   
   * Paucituberculata: B(0,0.1597)   
   * Didelphidae: B(0.11608,0.281)   
   * Dasyuromorphia: B(0.1597,0.33)   
   * Peramelidae: B(0.0436,0.238)     
   * Vombatiformes: B(0.255,0.33)   
   * Phalangeridae-Burramyidae: B(0.25,0.33)   
   * Petauridae-Pseudocheiridae: B(0.255,0.33)   
   * Macropodidae-Potoridae: B(0.247,0.33)   
   
<p align="center">
  <img width="1000" height="600" src="01_Check_conflict/01_SBnST_Marsupialia_MCMCruns.png">
</p>

**Deviations (main 72-taxa VS Marsupialia data sets)**   
<p align="center">
  <img width="1000" height="600" src="01_Check_conflict/01_SBnST_Marsupialia_meanquant.png">
</p>

The final tree topology can be found in the
[`final_tree_topology`](02_Final_tree_topology)
directory.

--- 

The next step is to run `MCMCtree` with the final tree topology and the 5-partitions 
alignment! Before that, however, we need to run `BASEML` to calculate the Hessian and 
the gradient, which are needed for the approximate likelihood calculation used by 
`MCMCtree` to speed up the Bayesian inference of divergence times.