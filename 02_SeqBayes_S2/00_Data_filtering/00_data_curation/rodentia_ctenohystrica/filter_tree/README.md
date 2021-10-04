# Rodentia ctenohystrica - phylogeny

## 1. Get tree topology and add calibrations
We use the R script [`Calibrations_Rctenohystrica.R`](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_ctenohystrica/filter_tree/00_Filter_trees/Calibrations_Rctenohystrica.R)
to generate the phylogeny for this data subset. Note that we use the
[`rodentia_ctenohystrica_rooted_calibnames.tree`](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_ctenohystrica/filter_tree/00_Filter_trees/rodentia_ctenohystrica_rooted_calibnames.tree)
file, where tag names have been manually added in the 
nodes that are to be calibrated. These tag names are later replaced with the
corresponding calibrations specified in the 
[`Calibrations_ctenohystrica.txt`](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_ctenohystrica/filter_tree/00_Filter_trees/Calibrations_ctenohystrica.txt)
file. 
In addition, this R script generates dummy alignments that can be used 
when running `MCMCtree` without the data to reduce disk space (see next section 3). 
This "dummy" alignment is saved [here](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/01_alignments/01_mammal_dummy_alns/rodentia_ctenohystrica) 
(now saved inside `before_updating_topology` directory; see below).

After running this script, you will have the following files:

```
00_Filter_trees
        |- RAxML_tree
        |         |- rodentia_ctenohystrica.tree              # File not used. Best-scoring ML tree obtained with RAxML
        |-extra_analyses         
        |         
        |- 208sp_rodentia_ctenohystrica_MCMCtree_calib.tree   # File output by the R script
        |- 208sp_rodentia_ctenohystrica_spnameslist.txt       # File output by the R script
        |- Calibrations_ctenohystrica.txt                     # Input file used by the R script. It matches the tag names
        |                                                     # in input tree with corresponding calibrations to be replaced
        |- Calibrations_Rctenohystrica.R                      # R script
        |- rodentia_ctenohystrica_rooted_baseml.tree          # File manually generated after running R script 
        |                                                     # to be used by BASEML (calibrations manually removed)
        |- rodentia_ctenohystrica_rooted_calibnames.tree      # Input file used by the R script
```

Note that we have manually generated the
[`rodentia_ctenohystrica_rooted_baseml.tree`](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_ctenohystrica/filter_tree/00_Filter_trees/rodentia_ctenohystrica_rooted_baseml.tree),
which does 
not contain the calibrations. This file was used before new taxa were added in the tree topology (see below).

## 2. Tree topology change 
After an extra data filtering when we added two extra taxa (i.e., *Ictidomys tridecemlineatus* 
and *Fukomys damarensis*), the tree topology changed to include the placement of these 
two taxa (see the details in
[this `README.md` file](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_ctenohystrica/filter_aln/README.md),
section `EXTRA FILTERING -- ADDING TAXA TO THE ALIGNMENT`, if you did not go through the data filtering before,
which explains why we added these two taxa 
and how this was done).

The updated file to be used by `BASEML` and the calibrated tree before the checks 
shown in the next step can be found
[here](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_ctenohystrica/filter_tree/00_Filter_trees/extra_analyses). 
The "dummy" alignments have also been updated in their corresponding directory 
[here](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/01_alignments/01_mammal_dummy_alns/rodentia_ctenohystrica),
where the previous version has been saved in a directory called `before_updating_topology`.


## 3. Check if calibrations are in conflict
The tree described above (210 taxa with *Ictidomys tridecemlineatus* 
and *Fukomys damarensis*) was used to find 
if there were any conflicts with the calibations used.
You can download the directories 
with the results obtained when running `MCMCtree` without the data
[here](https://www.dropbox.com/s/hqr0dvcvmvxc5tn/SeqBayesS2_check_conflict_ctenohystrica.zip?dl=0).
Once you download them, you should unzip its content and save the 
two directories inside the 
[`01_Check_conflict`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_ctenohystrica/filter_tree/01_Check_conflict)
directory so the file architecture is the following:

```
01_Check_conflict 
      |- 00_Prior_onlyST          # Provided in the zip file, not in this repository due to lack of space
      |- 01_Prior_SBandST         # Provided in the zip file, not in this repository due to lack of space
      |- outRdata                 # Directory  with R objects generated by the R script
      |- 00_Check_STanalitycalVSprior.R  # R script 
      |- *[pdf|png]                      # Files output by the R script, provided here
      |- *tsv                            # File output by the R script, provided here
```

Please read all the comments and explanations in
[the R script provided in this directory](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_ctenohystrica/filter_tree/01_Check_conflict/00_Check_STanalitycalVSprior.R) 
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

In this case (rodentia_ctenohystrica), we did not have to adjust any calibration as there were no
conflicts encountered (see plot below but also other plots within this directory):

**When using only ST calibrations**   
Calibrations used:   
   * Mammalia: ST(1.642,0.425,12.652,1714.565)   
   * Rodentia: ST(0.606,0.005,0.178,5.22)   
   * Caviomorpha-Phiomorpha: ST(0.414,0.014,-0.164,18.225)   
   * Phiomorpha: ST(0.324,0.016,-0.255,17.835)   
   * Caviomorpha: ST(0.356,0.014,0.022,18.214)   
   * *Cavia porcellus*-*Cavia aperea*: ST(0.086,0.01,0.671,36.22)   
   * *Chinchilla lanigera*-*Octodon degus*: ST(0.317,0.014,0.203,19.215)   
   
<p align="center">
  <img width="1000" height="600" src="https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_ctenohystrica/filter_tree/01_Check_conflict/00_Only_ST_RodCtenohystrica_MCMCruns.png">
</p>

**When using both ST and soft bound calibrations**   
Calibrations used:   
   * Mammalia: ST(1.642,0.425,12.652,1714.565)   
   * Rodentia: ST(0.606,0.005,0.178,5.22)   
   * Caviomorpha-Phiomorpha: ST(0.414,0.014,-0.164,18.225)   
   * Phiomorpha: ST(0.324,0.016,-0.255,17.835)   
   * Caviomorpha: ST(0.356,0.014,0.022,18.214)   
   * *Cavia porcellus*-*Cavia aperea*: ST(0.086,0.01,0.671,36.22)   
   * *Chinchilla lanigera*-*Octodon degus*: ST(0.317,0.014,0.203,19.215)   
   * Monotremata: B(0.1289,1.345)   
   * Tachyglossidae: B(0.0258,1.345)   
   * Abrocomidae: B(0.01778,0.09112)   
   
<p align="center">
  <img width="1000" height="600" src="https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_ctenohystrica/filter_tree/01_Check_conflict/01_SBnST_RodCtenohystrica_MCMCruns.png">
</p>

**Deviations (main 72-taxa VS rodentia_ctenohystrica data sets)**   
<p align="center">
  <img width="1000" height="600" src="https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_ctenohystrica/filter_tree/01_Check_conflict/01_SBnST_RodCtenohystrica_meanquant.png">
</p>

The final tree topology can be found in the
[`final_tree_topology`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_ctenohystrica/filter_tree/02_Final_tree_topology)
directory.

--- 

The next step is to run `MCMCtree` with the final tree topology and the 5-partitions 
alignment! Before that, however, we need to run `BASEML` to calculate the Hessian and 
the gradient, which are needed for the approximate likelihood calculation used by 
`MCMCtree` to speed up the Bayesian inference of divergence times.