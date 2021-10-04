# Rodentia the rest - phylogeny

## 1. Get tree topology and add calibrations
We use the R script [`Calibrations_Rtherest.R`](02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_therest/filter_tree/00_Filter_trees/Calibrations_Rtherest.R)
to generate the phylogeny for this data subset. Note that we use the
[`rodentia_therest_rooted_calibnames.tree`](02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_therest/filter_tree/00_Filter_trees/rodentia_therest_rooted_calibnames.tree)
file, where tag names have been manually added in the 
nodes that are to be calibrated. These tag names are later replaced with the
corresponding calibrations specified in the 
[`Calibrations_rodtherest.txt`](02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_therest/filter_tree/00_Filter_trees/Calibrations_rodtherest.txt)
file. 
In addition, this R script generates dummy alignments that can be used 
when running `MCMCtree` without the data to reduce disk space (see next section 3). 
This "dummy" alignment is saved [here](/02_SeqBayes_S2/00_Data_filtering/01_alignments/01_mammal_dummy_alns/rodentia_therest).

After running this script, you will have the following files:

```
00_Filter_trees
      |- RAxML_tree
      |         |- rodentia_therest.tree               # File not used. Best-scoring ML tree obtained with RAxML
      |- extra_filtering         
      |         
      |- 1314sp_Rodentia_therest_MCMCtree_calib.tree   # File output by the R script
      |- 1314sp_Rodentia_therest_spnameslist.txt       # File output by the R script
      |- Calibrations_Rtherest.R                       # R script
      |- Calibrations_rodtherest.txt                   # Input file used by the R script. It matches the tag names
      |                                                # in input tree with corresponding calibrations to be replaced
      |- rodentia_therest_rooted_baseml.tree           # File manually generated after running R script 
      |                                                # to be used by BASEML (calibrations manually removed)
      |- rodentia_therest_rooted_calibnames.tree       # Input file used by the R script
```

Note that we manually generated the
[`rodentia_therest_rooted_baseml.tree`](02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_therest/filter_tree/00_Filter_trees/rodentia_therest_rooted_baseml.tree),
which does 
not contain the calibrations. This file was used before we further split this data set into 
two data subsets (see below).

## 2. Generating subtree -- splitting the main tree into two
After partitioning the big data subset (1314 rodent species in "rodentia the rest") into two 
data subsets (see the details in
[this `README.md` file](02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_therest/filter_aln/README.md),
section `# EXTRA FILTERING -- DATA SUBSETTING`, if you did not go through the data filtering before,
which explains why we further partition "rodentia the rest").

The updated files to be used by `BASEML` and the calibrated trees before the checks 
shown in the next step can be found
[here](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_therest/filter_tree/00_Filter_trees/extra_filtering).

We manually generated the "dummy alignments" by including the extra taxa added in each subtree. 
Then, we saved them in the corresponding directories
([here](02_SeqBayes_S2/00_Data_filtering/01_alignments/01_mammal_dummy_alns/rodentia_subt1),
for the first subtree and 
[here](02_SeqBayes_S2/00_Data_filtering/01_alignments/01_mammal_dummy_alns/rodentia_subt2)
for the second subtree) and used them in the subsequent steps.

## 3. Check if calibrations are in conflict
The trees described above were used to find 
if there were any conflicts with the calibations used in both subtrees.
You can download the directories 
with the results obtained when running `MCMCtree` without the data
[here](https://www.dropbox.com/s/d2s8qfwhaso34kl/SeqBayesS2_check_conflict_rodsubt1.zip?dl=0)
for the first subtree and 
[here](https://www.dropbox.com/s/zog16ga689jhabv/SeqBayesS2_check_conflict_rodsubt2.zip?dl=0)
for the second subtree.
Once you download them, you should unzip its content and save the 
directories inside the corresponding
[`01_Check_conflict`](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_therest/filter_tree/01_Check_conflict)
directory for each subtree:

```
01_Check_conflict 
      |- 01_Check_conflict_rodsubtree1  # Save the unzipped content of rodentia subtree 1
      |- 01_Check_conflict_rodsubtree2  # Save the unzipped content of rodentia subtree 2
```

Please read all the comments and explanations in
the R scripts provided in the corresponding subdirectories above
to understand each step that we followed to avoid having conflicting calibrations in
the tree topology. Sometimes, we might need to adjust the ST calibrations and/or maximum
bounds if the neighbouring calibrations are in conflict (e.g., there are truncation issues). 

In a nutshell:   

   * 1. First, we run `MCMCtree` without using the data (i.e., 
   without using the alignment, hence the "dummy" alignment used here) and fixing the
   tree topology where only the skew-_t_ (ST) calibrations have been added.   
   * 2. For each calibrated node, we plot the corresponding analytical ST distribution
   (the one that we have told `MCMCtree` to use) against the corresponding posterior density
   inferred by `MCMCtree` when no data are used (data described in step 1). In addition,
   we add to this plot the posterior density of this node that was inferred by `MCMCtree`
   when using the first data set (72-taxon data set).   
   * 3. To check for conflict, we do the following for each calibrated node with an 
   ST calibration:   
      * Estimate mean times and quantiles (2.5% and 97.5%) from the posterior density
	  inferred by `MCMCtree` when the data are not used and the fixed tree topology has only
	  ST distributions.   
	  * Estimate mean times and quantiles from the posterior density inferred with
	  data set 1 (72-taxon data set) for the same node.   
	  * Check how much the former deviate from the latter.   
	  * If deviation is <0.6%, proceed with step 4.   
   * 4. If checks in step 3 are ok, we run `MCMCtree` without the data alignment but
   the tree topology now has both the ST calibrations and the calibrations with soft
   bounds (i.e., calibrations that have a minimum and a maximum bound with a 2.5% tail
   probability in each side).   
   * 5. Then, we generate the same plot as described in step 2.    
   * 4. Last, we check again for possible conflict as described in step 3. If deviation
   is <0.6% for all calibrated nodes, this is the end of the checks. Otherwise, we need 
   to adjust the location and scale parameters of the ST calibrations until no conflict
   is observed by subtracting the corresponding deviation (more details in the R script
   if this adjustment is taking place).   

## Rodentia subtree 1
In this case, we manually adjusted two ST calibrations (nodes labelled "*M. auratus*-*C. griseus*" and *Murinae*) 
that had been conflicting with other ST calibrations 
in previous analyses with the full "rodentia the rest" data set before we splitted it into two because `MCMCtree` 
could not deal with such a large amount of taxa. 
It was only *a posterior* that we realised that, if we followed the procedure you will read in other 
tutorials with other data subsets, we could adjust them but we do not really need to do this as the deviation 
is ~5% from the mean in Murinae and 5.5% in the low quantile for the other ST calibration (i.e., they 
are lower than our ~6% threshold, see 
[tsv file](02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_therest/filter_tree/01_Check_conflict/01_Check_conflict_rodsubtree1/01_SBnST_Rod_subtree1_sumstats.tsv)). 

We plotted both our manually adjusted calibrations based on 
previous analyses with the whole data set against the ST calibrations without adjustments:

<p align="center">
  <img width="700" height="600" src="02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_therest/filter_tree/01_Check_conflict/01_Check_conflict_rodsubtree1/fig_comp_1.png">
</p>   

>>Manually adjusted ST calibration (black) against ST calibration (red) for node   
>>`M. auratus - C. griseus`. In green and purple you have the densities for neighbouring nodes,   
>>`P. maniculatus - M. ochrogaster` and `Muridae`; respectively. All densities are adjusted to 1.   

<p align="center">
  <img width="700" height="600" src="02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_therest/filter_tree/01_Check_conflict/01_Check_conflict_rodsubtree1/fig_comp_2.png">
</p>   

>> Manually adjusted ST calibration (black) against ST calibration (red) for node   
>>`Murinae`. In green and purple you have the densities for neighbouring nodes,   
>> `N. galili_-Muridae` and `Muridae`; respectively. All densities are adjusted to 1.   
 
The differences are not large and should not be affecting downstream analyses.
The final tree topology used in the subsequent analyses can be found in the
[`final_tree_topology`](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_therest/filter_tree/02_Final_tree_topology)
directory.

## Rodentia subtree 2

**When using only ST calibrations**   
Calibrations used:   
   * Mammalia: ST(1.642,0.425,12.652,1714.565)   
   * Muridae: ST(0.174,0.035,3.117,381.97)   
   * Cricetidae: ST(0.137,0.03,3.464,240.317)   
   * Murinae: ST(0.0889,0.0241,3.7938,547.4454)     
   * *Mus pahari*-rest of Mus: ST(0.043,0.013,3.575,94.222)   
   * *Mus caroli*-rest of Mus: ST(0.022,0.007,3.426,86.245)   
   * _Mus spretus_-_Mus musculus_: ST(0.0104,0.0033,3.3802,83.5293)   
   
<p align="center">
  <img width="1000" height="600" src="02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_therest/filter_tree/01_Check_conflict/01_Check_conflict_rodsubtree2/00_Only_ST_Rod_subtree2_MCMCruns.png">
</p>

**When using both ST and soft bound calibrations**   
Calibrations used:   
   * Mammalia: ST(1.642,0.425,12.652,1714.565)   
   * Muridae: ST(0.174,0.035,3.117,381.97)   
   * Cricetidae: ST(0.137,0.03,3.464,240.317)   
   * **Murinae: ~ST(0.0889,0.0241,3.7938,547.4454)~ -- CONFLICT, adjusted**   
   **to ST(0.085,0.021,3.794,547.445) for next round**   
   * *Mus pahari*-rest of Mus: ST(0.043,0.013,3.575,94.222)   
   * *Mus caroli*-rest of Mus: ST(0.022,0.007,3.426,86.245)   
   * **_Mus spretus_-_Mus musculus_: ~ST(0.0104,0.0033,3.3802,83.5293)~ -- CONFLICT, adjusted**   
   **to ST(0.0101,0.003,3.38,83.529) for next round**   
   
   
<p align="center">
  <img width="1000" height="600" src="02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_therest/filter_tree/01_Check_conflict/01_Check_conflict_rodsubtree2/01_SBnST_Rod_subtree2_MCMCruns.png">
</p>

**When using both ST and soft bound calibrations - 2nd round**   
Calibrations used:   
   * Mammalia: ST(1.642,0.425,12.652,1714.565)   
   * Muridae: ST(0.174,0.035,3.117,381.97)   
   * Cricetidae: ST(0.137,0.03,3.464,240.317)   
   * Murinae: ST(0.085,0.021,3.794,547.445)   
   * *Mus pahari*-rest of Mus: ST(0.043,0.013,3.575,94.222)   
   * *Mus caroli*-rest of Mus: ST(0.022,0.007,3.426,86.245)   
   * _Mus spretus_-_Mus musculus_: ST(0.0101,0.003,3.38,83.529)   
   
<p align="center">
  <img width="1000" height="600" src="02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_therest/filter_tree/01_Check_conflict/01_Check_conflict_rodsubtree2/02_SBandSTtweak1_Rod_subtree2_MCMCruns.png">
</p>


**Deviations (main 72-taxa VS rodentia_therest data sets)**   
<p align="center">
  <img width="1000" height="600" src="02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_therest/filter_tree/01_Check_conflict/01_Check_conflict_rodsubtree2/02_SBnSTtweak1_Rod_subtree2_meanquant.png">
</p>

The final tree topology can be found in the
[`final_tree_topology`](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_therest/filter_tree/02_Final_tree_topology)
directory.

--- 

The next step is to run `MCMCtree` with the final tree topology and the 5-partitions 
alignment! Before that, however, we need to run `BASEML` to calculate the Hessian and 
the gradient, which are needed for the approximate likelihood calculation used by 
`MCMCtree` to speed up the Bayesian inference of divergence times.