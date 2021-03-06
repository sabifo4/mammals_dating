# Euarchonta - phylogeny

## 1. Get tree topology and add calibrations
### Setting maximum age for nodes that only had a minimum assigned following the fossil record
We note that the node for Primatomorpha needs to have the maximum age established. For that purpose, 
we found the ancestor of this node that was calibrated with an ST distribution and computed the 
2.5% quantile as it follows:

```r 
# The ancestor is "Euarchontoglires", which ST 
# distributions is 'ST(0.694522263625083,0.00703784739084366,0.319536189703347,7.61865839547998)'
# We use the `sn` package as it follows:
sn::qst(0.025, 0.694522263625083,0.00703784739084366,0.319536189703347,7.61865839547998)
# 0.6809366
```

According to this, the soft-bound calibration for Primatomorpha is `B(0.64645,0.68)`.

### Generating calibrated trees
We use the R script [`Calibrations_Euarchonta.R`](00_Filter_trees/Calibrations_Euarchonta.R)
to generate the phylogeny for this data subset. Note that we use the
[`euarchonta_rooted_calibnames.tree`](00_Filter_trees/euarchonta_rooted_calibnames.tree)
file, where tag names have been manually added in the 
nodes that are to be calibrated. These tag names are later replaced with the
corresponding calibrations specified in the 
[`Calibrations_Euarchonta.txt`](00_Filter_trees/Calibrations_Euarchonta.txt)
file. 
In addition, this R script generates dummy alignments that can be used 
when running `MCMCtree` without the data to reduce disk space (see next section 3). 
This "dummy" alignment is saved [here](../../../01_alignments/01_mammal_dummy_alns/euarchonta).

After running this script, you will have the following files:

```
00_Filter_trees 
     |- RAxML_tree
     |         |- euarchonta.tree             # File not used. Best-scoring ML tree obtained with RAxML
     |         
     |- 486sp_Euarchonta_MCMCtree_calib.tree  # File output by the R script
     |- 486sp_Euarchonta_spnameslist.txt      # File output by the R script
     |- Calibrations_Euarchonta.R             # R script
     |- Calibrations_Euarchonta.txt           # Input file used by the R script. It matches the tag names
     |                                        # in input tree with corresponding calibrations to be replaced
     |- euarchonta_rooted_baseml.tree         # File manually generated after running R script 
     |                                        # to be used by BASEML (calibrations manually removed)
     |- euarchonta_rooted_calibnames.tree     # Input file used by the R script
```

Note that we have manually generated the
[`euarchonta_rooted_baseml.tree`](00_Filter_trees/euarchonta_rooted_baseml.tree),
which does 
not contain the calibrations. This file will be used when running `BASEML` to compute 
the Hessian and the gradient that are needed by `MCMCtree` to run the approximate 
likelihood.


## 2. Check if calibrations are in conflict
Now, we had to check if there were any conflicts with the calibations used.
You can download the directories 
with the results obtained when running `MCMCtree` without the data
[here](https://www.dropbox.com/s/l17ct03llwhrttj/SeqBayesS2_check_conflict_euarchonta.zip?dl=0).
Once you download them, you should unzip its content and save the 
directories inside the 
[`01_Check_conflict`](01_Check_conflict)
directory so the file architecture is the following:

```
01_Check_conflict 
      |- 00_Prior_onlyST         # Provided in the zip file, not in this repository due to lack of space
      |- 01_Prior_SBandST        # Provided in the zip file, not in this repository due to lack of space
      |- 02_Prior_SBandSTtweak1  # Provided in the zip file, not in this repository due to lack of space
      |- 03_Prior_SBandSTtweak2  # Provided in the zip file, not in this repository due to lack of space
      |- outRdata                # Directory  with R objects generated by the R script
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

In this case (euarchonta), we had to undergo several rounds of adjusting calibrations as there were
conflicts encountered (see plots below but also other plots within this directory):


**When using only ST calibrations**   
Calibrations used:   
   * Mammalia: ST(1.642,0.425,12.652,1714.565)  
   * Euarchontoglires: ST(0.695,0.007,0.32,7.619)   
   * Primates: ST(0.655,0.01,-1.355,178.316)   
   * Strepsirrhini: ST(0.548,0.026,-2.506,66.983)    
   * Propithecus-Microcebus: ST(0.37,0.033,-0.876,275.655)    
   * Haplorrhini: ST(0.622,0.011,-1.196,166.803)    
   * Anthropoidea: ST(0.415,0.021,-1.14,156.796)    
   * Aotidae-Callitrichidae: ST(0.1996,0.0266,-1.8564,48.1352)    
   * Cebidae: ST(0.1754,0.0237,-1.7218,40.5053)    
   * Catarrhini: ST(0.314,0.018,-1.23,314.095)    
   * Cercopithecoidea: ST(0.182,0.012,-0.054,88.157)    
   * Cercopithecinae: ST(0.1363,0.0093,-0.0000,10.0000)    
   * Papionini: ST(0.1,0.009,0.641,145.414)    
   * Papio-Mandrillus: ST(0.087,0.007,0.195,50.335)    
   * *Cercocebus atys*-*Mandrillus leucophaeus*: ST(0.073,0.006,0.229,40.768)    
   * Genus macaca: ST(0.053,0.007,1.101,200.596)    
   * *Macaca fascicularis*-*Macaca mulatta*: ST(0.039,0.005,1.207,108.063)    
   * Colobinae: ST(0.127,0.011,0.584,209.998)  
   * _R. roxellana_-_R. bieti: ST(0.0363,0.0061,1.5249,104.3879)   
   * Hominoidea: ST(0.236,0.016,-1.223,135.248)   
   * Hominidae: ST(0.21,0.015,-1.248,109.524)   
   * Homininae: ST(0.122,0.012,-4.859,295.449)   
   * Hominini: ST(0.101,0.01,-7.603,93.226)   
   * *P. paniscus*-*P. troglodites*: ST(0.039,0.003,-0.337,47.276)     

<p align="center">
  <img width="1000" height="600" src="01_Check_conflict/00_Only_ST_Euarchonta_MCMCruns.png">
</p>


**When using both ST and soft bound calibrations**   
Calibrations used:   
   * Mammalia: ST(1.642,0.425,12.652,1714.565)  
   * Euarchontoglires: ST(0.695,0.007,0.32,7.619)   
   * Primates: ST(0.655,0.01,-1.355,178.316)   
   * Strepsirrhini: ST(0.548,0.026,-2.506,66.983)    
   * Propithecus-Microcebus: ST(0.37,0.033,-0.876,275.655)    
   * Haplorrhini: ST(0.622,0.011,-1.196,166.803)    
   * Anthropoidea: ST(0.415,0.021,-1.14,156.796)    
   * **Aotidae-Callitrichidae: ~ST(0.1996,0.0266,-1.8564,48.1352)~ -- CONFLICT, adjusted**
   **to ST(0.1933,0.0289,-1.8564,48.1352) for next round**    
   * **Cebidae: ~ST(0.1754,0.0237,-1.7218,40.5053)~ -- CONFLICT, adjusted**
   **to ST(0.1892,0.0207,-1.7218,40.5053) for next round**   
   * Catarrhini: ST(0.314,0.018,-1.23,314.095)    
   * Cercopithecoidea: ST(0.182,0.012,-0.054,88.157)    
   * Cercopithecinae: ST(0.1363,0.0093,-0.0000,10.0000)    
   * Papionini: ST(0.1,0.009,0.641,145.414)    
   * Papio-Mandrillus: ST(0.087,0.007,0.195,50.335)    
   * *Cercocebus atys*-*Mandrillus leucophaeus*: ST(0.073,0.006,0.229,40.768)    
   * Genus macaca: ST(0.053,0.007,1.101,200.596)    
   * *Macaca fascicularis*-*Macaca mulatta*: ST(0.039,0.005,1.207,108.063)    
   * Colobinae: ST(0.127,0.011,0.584,209.998)  
   * **_R. roxellana_-_R. bieti: ~ST(0.0363,0.0061,1.5249,104.3879)~ -- CONFLICT, adjusted**   
   **to ST(0.0385,0.0077,1.5249,104.3879) for next round**   
   * Hominoidea: ST(0.236,0.016,-1.223,135.248)   
   * Hominidae: ST(0.21,0.015,-1.248,109.524)   
   * Homininae: ST(0.122,0.012,-4.859,295.449)   
   * Hominini: ST(0.101,0.01,-7.603,93.226)   
   * **_P. paniscus_ - _P. troglodites_: ~ST(0.039,0.003,-0.337,47.276)~ -- CONFLICT, adjusted**   
   **to ST(0.038,0.004,-0.337,47.276) for next round**   
   * Monotremata: B(0.2446,1.332)   
   * Tachyglossidae: B(0.0258,1.332)   
   * Platyrrhini: B(0.2045,0.377)   
   * Primatomorpha: B(0.64645,0.68)   
   * Scandentia: B(0.38,0.66)   
   
<p align="center">
  <img width="1000" height="600" src="01_Check_conflict/01_SBnST_Euarchonta_MCMCruns.png">
</p>


**When using both ST and soft bound calibrations - 2nd round**   
Calibrations used:   
   * Mammalia: ST(1.642,0.425,12.652,1714.565)  
   * Euarchontoglires: ST(0.695,0.007,0.32,7.619)   
   * Primates: ST(0.655,0.01,-1.355,178.316)   
   * Strepsirrhini: ST(0.548,0.026,-2.506,66.983)    
   * Propithecus-Microcebus: ST(0.37,0.033,-0.876,275.655)    
   * Haplorrhini: ST(0.622,0.011,-1.196,166.803)    
   * Anthropoidea: ST(0.415,0.021,-1.14,156.796)    
   * **Aotidae-Callitrichidae: ~ST(0.1933,0.0289,-1.8564,48.1352)~ -- CONFLICT, adjusted**   
   **to ST(0.1861,0.0329,-1.8564,48.1352) for next round**    
   * Cebidae: ST(0.1754,0.0237,-1.7218,40.5053) 
   * Catarrhini: ST(0.314,0.018,-1.23,314.095)    
   * Cercopithecoidea: ST(0.182,0.012,-0.054,88.157)    
   * Cercopithecinae: ST(0.1363,0.0093,-0.0000,10.0000)    
   * Papionini: ST(0.1,0.009,0.641,145.414)    
   * Papio-Mandrillus: ST(0.087,0.007,0.195,50.335)    
   * *Cercocebus atys*-*Mandrillus leucophaeus*: ST(0.073,0.006,0.229,40.768)    
   * Genus macaca: ST(0.053,0.007,1.101,200.596)    
   * *Macaca fascicularis*-*Macaca mulatta*: ST(0.039,0.005,1.207,108.063)    
   * Colobinae: ST(0.127,0.011,0.584,209.998)  
   * **_R. roxellana_-_R. bieti: ~ST(0.0385,0.0077,1.5249,104.3879)~ -- CONFLICT, adjusted**   
   **to ST(0.0375,0.0067,1.5249,104.3879) for next round**    
   * Hominoidea: ST(0.236,0.016,-1.223,135.248)   
   * Hominidae: ST(0.21,0.015,-1.248,109.524)   
   * Homininae: ST(0.122,0.012,-4.859,295.449)   
   * Hominini: ST(0.101,0.01,-7.603,93.226)   
   * **_P. paniscus_ - _P. troglodites_: ~ST(0.038,0.004,-0.337,47.276)~ -- CONFLICT, adjusted**   
   **to ST(0.0375,0.032,-0.337,47.276) for next round**   
   * Monotremata: B(0.2446,1.332)   
   * Tachyglossidae: B(0.0258,1.332)   
   * Platyrrhini: B(0.2045,0.377)   
   * Primatomorpha: B(0.64645,0.68)   
   * Scandentia: B(0.38,0.66)   
   
<p align="center">
  <img width="1000" height="600" src="01_Check_conflict/02_SBandSTtweak1_Euarchonta_MCMCruns.png">
</p>


**When using both ST and soft bound calibrations - 3rd round**   
Calibrations used:   
   * Mammalia: ST(1.642,0.425,12.652,1714.565)  
   * Euarchontoglires: ST(0.695,0.007,0.32,7.619)   
   * Primates: ST(0.655,0.01,-1.355,178.316)   
   * Strepsirrhini: ST(0.548,0.026,-2.506,66.983)    
   * Propithecus-Microcebus: ST(0.37,0.033,-0.876,275.655)    
   * Haplorrhini: ST(0.622,0.011,-1.196,166.803)    
   * Anthropoidea: ST(0.415,0.021,-1.14,156.796)    
   * Aotidae-Callitrichidae: ST(0.1861,0.0329,-1.8564,48.1352)   
   * Cebidae: ST(0.1754,0.0237,-1.7218,40.5053) 
   * Catarrhini: ST(0.314,0.018,-1.23,314.095)    
   * Cercopithecoidea: ST(0.182,0.012,-0.054,88.157)    
   * Cercopithecinae: ST(0.1363,0.0093,-0.0000,10.0000)    
   * Papionini: ST(0.1,0.009,0.641,145.414)    
   * Papio-Mandrillus: ST(0.087,0.007,0.195,50.335)    
   * *Cercocebus atys*-*Mandrillus leucophaeus*: ST(0.073,0.006,0.229,40.768)    
   * Genus macaca: ST(0.053,0.007,1.101,200.596)    
   * *Macaca fascicularis*-*Macaca mulatta*: ST(0.039,0.005,1.207,108.063)    
   * Colobinae: ST(0.127,0.011,0.584,209.998)  
   * _R. roxellana_-_R. bieti: ST(0.0375,0.0067,1.5249,104.3879)   
   * Hominoidea: ST(0.236,0.016,-1.223,135.248)   
   * Hominidae: ST(0.21,0.015,-1.248,109.524)   
   * Homininae: ST(0.122,0.012,-4.859,295.449)   
   * Hominini: ST(0.101,0.01,-7.603,93.226)   
   * _P. paniscus_ - _P. troglodites_: ST(0.0375,0.032,-0.337,47.276)   
   * Monotremata: B(0.2446,1.332)   
   * Tachyglossidae: B(0.0258,1.332)   
   * Platyrrhini: B(0.2045,0.377)   
   * Primatomorpha: B(0.64645,0.68)   
   * Scandentia: B(0.38,0.66)   
   
<p align="center">
  <img width="1000" height="600" src="01_Check_conflict/03_SBandSTtweak2_Euarchonta_MCMCruns.png">
</p>

**Deviations (main 72-taxa VS Euarchonta data sets)**   
<p align="center">
  <img width="1000" height="600" src="01_Check_conflict/03_SBnSTtweak2_Euarchonta_meanquant.png">
</p>

The final tree topology can be found in the
[`final_tree_topology`](02_Final_tree_topology)
directory.

--- 

The next step is to run `MCMCtree` with the final tree topology and the 5-partitions 
alignment! Before that, however, we need to run `BASEML` to calculate the Hessian and 
the gradient, which are needed for the approximate likelihood calculation used by 
`MCMCtree` to speed up the Bayesian inference of divergence times.