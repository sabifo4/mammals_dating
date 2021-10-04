# Chiroptera - phylogeny

## 1. Get tree topology and add calibrations
We use the R script [`Calibrations_Lchiroptera.R`](00_Filter_trees/Calibrations_Lchiroptera.R)
to generate the phylogeny for this data subset. Note that we use the
[`laurasiatheria_chiroptera_calibnames.tree`](00_Filter_trees/laurasiatheria_chiroptera_calibnames.tree)
file, where tag names have been manually added in the 
nodes that are to be calibrated. These tag names are later replaced with the
corresponding calibrations specified in the 
[`Calibrations_LaurChiroptera.txt`](00_Filter_trees/Calibrations_LaurChiroptera.txt)
file. 

After running this script, you will have the following files:

```
00_Filter_trees
      |- extra_filtering
      |- RAxML_tree
      |         |- Laurasiatheria_chiroptera.tree              # File not used. Best-scoring ML tree obtained with RAxML
      |         
      |- 885sp_Laurasiatheria_chiroptera_MCMCtree_calib.tree   # File output by the R script
      |- 885sp_Laurasiatheria_chiroptera_spnameslist.txt       # File output by the R script
      |- Calibrations_LaurChiroptera.txt                       # Input file used by the R script. It matches the tag names
      |                                                        # in input tree with corresponding calibrations to be replaced
      |- Calibrations_LChiroptera.R                            # R script
      |- laurasiatheria_chiroptera_calibnames.tree             # Input file used by the R script
      |- laurasiatheria_chiroptera_rooted_baseml.tree          # File manually generated after running R script 
                                                               # to be used by BASEML (calibrations manually removed)
```

Note that we manually generated the
[`laurasiatheria_chiroptera_rooted_baseml.tree`](00_Filter_trees/laurasiatheria_chiroptera_rooted_baseml.tree),
which does 
not contain the calibrations. This file was used when running `BASEML` to compute 
the Hessian and the gradient that are needed by `MCMCtree` to run the approximate 
likelihood before we realised that this data subset needed to be further split 
into two.

## 2. Generating subtree -- splitting the main tree into two
After partitioning the big data subset (885 bat species in "chiroptera" data subset) into two 
data subsets (see the details in
[this `README.md` file](../filter_aln/README.md),
section `# EXTRA FILTERING -- DATA SUBSETTING`, if you did not go through the data filtering before,
which explains why we further partition "chiroptera").

The updated files to be used by `BASEML` and the calibrated trees before the checks 
shown in the next step can be found
[here](00_Filter_trees/extra_filtering).

We manually generated the "dummy alignments" by including the extra taxa added in each subtree. 
Then, we saved them in the corresponding directories
(find [here](../../../01_alignments/01_mammal_dummy_alns/chiroptera_subt1)
the one for the first subtree and 
[here](../../../01_alignments/01_mammal_dummy_alns/chiroptera_subt2)
the one for the second subtree) and used them in the subsequent steps.

## 3. Check if calibrations are in conflict
The trees described above were used to find 
if there were any conflicts with the calibations used in both subtrees.
You can download the directories 
with the results obtained when running `MCMCtree` without the data
[here](https://www.dropbox.com/s/69mbdkvdri0oeii/SeqBayesS2_check_conflict_chirosubt1.zip?dl=0)
for the first subtree and 
[here](https://www.dropbox.com/s/yv5uhhajg71ibr9/SeqBayesS2_check_conflict_chirosubt2.zip?dl=0)
for the second subtree.
Once you download them, you should unzip their content inside the corresponding
[`01_Check_conflict`](01_Check_conflict)
directory for each subtree:

```
01_Check_conflict 
      |- 01_Check_conflict_chiroptera_subtree1  # Save the unzipped content of chiroptera subtree 1
      |- 01_Check_conflict_chiroptera_subtree2  # Save the unzipped content of chiroptera subtree 2
```

Please read all the comments and explanations in
the R scripts provided in the corresponding subdirectories above
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
	  * If deviation is <0.6%, proceed with step 4.   
   4. If checks in step 3 are ok, we run `MCMCtree` without the data alignment but
   the tree topology now has both the ST calibrations and the calibrations with soft
   bounds (i.e., calibrations that have a minimum and a maximum bound with a 2.5% tail
   probability in each side).   
   5. Then, we generate the same plot as described in step 2.    
   6. Last, we check again for possible conflict as described in step 3. If deviation
   is <0.6% for all calibrated nodes, this is the end of the checks. Otherwise, we need 
   to adjust the location and scale parameters of the ST calibrations until no conflict
   is observed by subtracting the corresponding deviation (more details in the R script
   if this adjustment is taking place).   

## Chiroptera subtree 1
**When using only ST calibrations**   
Calibrations used:   
   * Mammalia: ST(1.642,0.425,12.652,1714.565)   
   * Chiroptera: ST(0.596,0.016,-1.239,13.572)   
   
<p align="center">
  <img width="1000" height="600" src="01_Check_conflict/01_Check_conflict_chiroptera_subtree1/00_Only_ST_L.chiroptera_subt1_MCMCruns.png">
</p>

**When using both ST and soft bound calibrations**   
Calibrations used:   
   * Mammalia: ST(1.642,0.425,12.652,1714.565)   
   * Chiroptera: ST(0.596,0.016,-1.239,13.572)   
   * Monotremata: B(0.1289,1.345)   
   * Tachyglossidae: B(0.0258,1.345)   
   * Hipposideridae-Rhinolophidae: B(0.38,0.56)   
   * Craseonycteridae-Megadermatidae: B(0.339,0.478)   
  
<p align="center">
  <img width="1000" height="600" src="01_Check_conflict/01_Check_conflict_chiroptera_subtree1/01_SBnST_L.chiroptera_subt1_MCMCruns.png">
</p>

**Deviations (main 72-taxa VS chiroptera_subt1 data sets)**   
<p align="center">
  <img width="1000" height="600" src="01_Check_conflict/01_Check_conflict_chiroptera_subtree1/01_SBnST_L.chiroptera_subt1_meanquant.png">
</p>

The final tree topology can be found in the
[`final_tree_topology`](/02_Final_tree_topology)
directory.

## Chiroptera subtree 2

**When using only ST calibrations**   
Calibrations used:   
   * Mammalia: ST(1.642,0.425,12.652,1714.565)   
   * Chiroptera: ST(0.596,0.016,-1.239,13.572)   
   
<p align="center">
  <img width="1000" height="600" src="01_Check_conflict/01_Check_conflict_chiroptera_subtree2/00_Only_ST_L.chiroptera_subt2_MCMCruns.png">
</p>

**When using both ST and soft bound calibrations**   
Calibrations used:   
   * Mammalia: ST(1.642,0.425,12.652,1714.565)   
   * Chiroptera: ST(0.596,0.016,-1.239,13.572)   
   * Monotremata: B(0.1289,1.345)   
   * Tachyglossidae: B(0.0258,1.345)   
   * Molossidae-V-M: B(0.38,0.56)   
   * Molossidae-V-M-N: B(0.38,0.56)   
  
<p align="center">
  <img width="1000" height="600" src="01_Check_conflict/01_Check_conflict_chiroptera_subtree2/01_SBnST_L.chiroptera_subt2_MCMCruns.png">
</p>

**Deviations (main 72-taxa VS chiroptera_subt2 data sets)**   
<p align="center">
  <img width="1000" height="600" src="01_Check_conflict/01_Check_conflict_chiroptera_subtree2/01_SBnST_L.chiroptera_subt2_meanquant.png">
</p>

The final tree topology can be found in the
[`final_tree_topology`](02_Final_tree_topology)
directory.

--- 

The next step is to run `MCMCtree` with the final tree topology and the 5-partitions 
alignment! Before that, however, we need to run `BASEML` to calculate the Hessian and 
the gradient, which are needed for the approximate likelihood calculation used by 
`MCMCtree` to speed up the Bayesian inference of divergence times.