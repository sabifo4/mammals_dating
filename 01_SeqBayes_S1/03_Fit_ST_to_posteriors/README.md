# 1. Fit ST-distributions 
Once we have summarised the posterior estimates, we can fit skew-_t_ (ST) distributions to each of the
mean posterior time estimates for the 71 internal nodes of the mammal phylogeny. To do this, we have
written the R script [`00_Fit_skewT.R`](/01_SeqBayes_S1/03_Fit_ST_to_posteriors/00_Fit_skewT.R),
which uses the R function `sn::st.mple` ([Azzalini, 2020](http://azzalini.stat.unipd.it/SN))
for that purpose. In order to reproduce the results, a seed number (i.e., `12345`) has been set.

Note that the method that optimises the finding of the best fitted ST distributions to the posterior
estimates for each node is the `BFGS` ([Nash, 1990](https://www.taylorfrancis.com/books/9781315139784)).
We use a `for` loop limited to 50 iterations where, during each iteration, the 
`BFGS` optimization approach tries to find the best ST distribution fitted to the corresponding node. 
If this approach was to fail after 50 iterations, then the method by default in that R function (`nlminb`)
is used. 

> **NOTE**: We first used the `nlminb` method ([Azzalini, 2020](http://azzalini.stat.unipd.it/SN)) as a
>first choice of optimization search. Nevertheless, it was difficult to converge for some of the deeper
>nodes in the mammals phylogeny, while the `BFGS` approach would find it either in the first iteration 
>or after 2 or 3 searches. Therefore, we decided to use first the `BFGS` optimization algorithm 
>and, if it failed, then `nlminb` would be used. 

The file architecture is organised as follows:   

```
03_Fit_ST_to_posteriors/ 
           |- 00_fitST/
           |       |- logs/
           |       |- plots/ 
           |       |- Rdata/   <-- Only 2 out of the 4 objects are provided in the 
           |       |               repository (lack of space). Objects `prior.divtimes.RData` 
           |       |               and `post.divtimes.RData` can be generated if the script 
           |       |               is re-run or downloading this file: https://www.dropbox.com/s/fiyunh6puaro3kr/SeqBayesS1_FitST_old_Rdata_priorandpost.zip?dl=0	   
           |       |- Rout/     			   
           |
           |- 01_trees_with_ST_calibs/
           |       |- Rinp/
           |       |- Rout/ 
           |
           |- 02_MCMCtree_prior_ST/
           |       |- MCMCtree/
           |       |- plots/ 
           |
           |- 00_Fit_skewT.R
           |- 01_Add_STcalibs_to_tree.R
           |- 02_Eval_skewT.R
           |- README.md
```

When running the code to fit ST distributions to each posterior node age, the following takes place:   

   * [`logs`](/01_SeqBayes_S1/03_Fit_ST_to_posteriors/00_fitST/logs):
   the [`log_file_convergence_BFGS.txt`](/01_SeqBayes_S1/03_Fit_ST_to_posteriors/00_fitST/logs/log_file_convergence_BFGS.txt)
   file is generated. It contains a record of the trials needed to fit a ST distribution to a given node.
   If only one iteration has been needed, as it has happened in this analysis, only the header
   `Working with node t_nXXX...` is printed out.   
   * [`plots`](/01_SeqBayes_S1/03_Fit_ST_to_posteriors/00_fitST/plots):
   we generate one plot per node in which the fitted ST distribution is drawn on top of the 
   posterior distribution. The resulting images are used to evaluate how good the ST fits the posterior
   distribution.   
   * [`Rdata`](/01_SeqBayes_S1/03_Fit_ST_to_posteriors/00_fitST/Rdata):
   we generate two R objects, [`ST.fitted.dists.RData`](/01_SeqBayes_S1/03_Fit_ST_to_posteriors/00_fitST/Rdata/ST.fitted.dists.RData)
   and [`ST.fitted.objects.RData`](/01_SeqBayes_S1/03_Fit_ST_to_posteriors/00_fitST/Rdata/ST.fitted.objects.RData)
   with the resulting ST distributions. This will allow us to later load
   them in other R scripts without the need to re-run the search. Objects
   `prior.divtimes.RData` and `post.divtimes.RData` can be generated if the script
   is re-run or downloading [this file](https://www.dropbox.com/s/fiyunh6puaro3kr/SeqBayesS1_FitST_old_Rdata_priorandpost.zip?dl=0).   
   * [`Rout`](/01_SeqBayes_S1/03_Fit_ST_to_posteriors/00_fitST/Rout):
   this directory contains a tab separated file with the values of the parameters of each ST 
   distribution, [`ST.fitted.dists.G2.40.tsv`](/01_SeqBayes_S1/03_Fit_ST_to_posteriors/00_fitST/Rout/ST.fitted.dists.G2.40.tsv).
   The output file `ST.fitted.dists.G2.40.tsv` remains in that main directory as it is subsequently used by
   the script [`01_Add_STcalibs_to_tree.R`](/01_SeqBayes_S1/03_Fit_ST_to_posteriors/01_Add_STcalibs_to_tree.R)
   to add the corresponding ST calibrations to the newick tree.

# 2. Add ST calibrations in `MCMCtree` format 
We then use the R script [`01_Add_STcalibs_to_tree.R`](/01_SeqBayes_S1/03_Fit_ST_to_posteriors/01_Add_STcalibs_to_tree.R)
to add the ST calibrations in `MCMCtree` format to the tree file. This script will be using the following input files:   

   * [`FigTree_72sp_nodelabels.tree`](/01_SeqBayes_S1/03_Fit_ST_to_posteriors/01_trees_with_ST_calibs/Rinp/FigTree_72sp_nodelabels.tree): the tree file with numerical labels that indiciate the node positions. This newick tree 
   was obtained from the `out.txt` file output by `MCMCtree` during the estimation of posterior times.   
   * [`Node_calibrations_mammalia_72sp.txt`](/01_SeqBayes_S1/03_Fit_ST_to_posteriors/01_trees_with_ST_calibs/Rinp/Node_calibrations_mammalia_72sp.txt): mapping file with three columns: (i) the name of the clades for which 
   calibrations had been used to in the first Bayesian analysis, (ii) the soft bound calibrations used in this 
   analysis, and (iii) the numbers corresponding to the node position that match the newick tree described above.   
   * [`ST.fitted.dists.G2.40.tsv`](/01_SeqBayes_S1/03_Fit_ST_to_posteriors/00_fitST/Rout/ST.fitted.dists.G2.40.tsv): fitted ST distributions using the `sn::st.mple` R function 
   as described in step 1 above. It contains several columns: (i) the number corresponding to 
   the node position in `MCMCtree` format, i.e., `t_nXX`; (ii) four columns with the four parameters of the 
   ST distribution, and (iii) two columns with the ST distribution written in `MCMCtree` format to be appended in the tree with 
   all the decimals or with just three decimals. 

Last, two output files are created:   
   * [`72sp_MAMMALS_atlantogenata_tarver2016_71STcalib.tree`](/01_SeqBayes_S1/03_Fit_ST_to_posteriors/01_trees_with_ST_calibs/Rout/72sp_MAMMALS_atlantogenata_tarver2016_71STcalib.tree): this file has all the nodes  
   that were calibrated with soft bound distributions in the first Bayesian analysis with the 
   corresponding fitted ST distributions.   
   * [`72sp_MAMMALS_atlantogenata_tarver2016_71STcalib_rounded.tree`](/01_SeqBayes_S1/03_Fit_ST_to_posteriors/01_trees_with_ST_calibs/Rout/72sp_MAMMALS_atlantogenata_tarver2016_71STcalib_rounded.tree): same content that the file 
   described above has, but the values for the parameters have been rounded to three decimals.   

# 3. Run prior 
The tree with the ST distributions generated in the previous step (71 ST distributions) is used 
again by `MCMCtree` to collect samples from the prior. This is done to check that the ST distributions fitted
in R are sensible and correspond to our expectations (i.e., there is not a mismatch between the "user-specified" priors and the "effective" priors that `MCMCtree` will be using as node calibrations).
All analyses regarding this step can be found in the directory [`02_MCMCtree_prior_ST`](/01_SeqBayes_S1/03_Fit_ST_to_posteriors/02_MCMCtree_prior_ST). Its structure is the following:   

   * [**`MCMCtree`**](/01_SeqBayes_S1/03_Fit_ST_to_posteriors/02_MCMCtree_prior_ST/MCMCtree)   
   This directory contains the input files used to run `MCMCtree` [here](/01_SeqBayes_S1/03_Fit_ST_to_posteriors/02_MCMCtree_prior_ST/MCMCtree/inp_files)
   as well as (i) a file with the inferred dated tree in newick format, (ii) the seed number for each MCMC run, and (iii) another file with the corresponding samples collected for the parameters of estimated.
   The resulting phylogeny and the samples collected for the two chains ran can be found here too. Note that there are other files output by `MCMCtree` but,
   due to the limited size a file can have in a repository, they have not been included here.   
   * [**`plots`**](/01_SeqBayes_S1/03_Fit_ST_to_posteriors/02_MCMCtree_prior_ST/plots)   
   The R script [`02_Eval_skewT.R`](https://github.com/sabifo4/mammals_dating/blob/main/01_SeqBayes_S1/03_Fit_ST_to_posteriors/02_Eval_skewT.R) is then used to generate several plots in which we add the fitted 
   ST distributions previously calculated (i.e., "analytical ST", which were saved in `RData` files [here](/01_SeqBayes_S1/03_Fit_ST_to_posteriors/00_fitST/Rdata) and are now loaded
   in this step) against the estimated distributions that result from averaging across all the collected samples gathered by `MCMCtree` when sampling from the prior (i.e., the estimates on prior times
   when calibrating the tree with the analytical ST distributions) and from the posterior (i.e., the posterior estimates inferred during the Bayesian dating analysis described [here](/01_SeqBayes_S1/02_MCMCtree)).   

## [[IMPORTANT NOTES]] 
**NOTE 1: Unfortunately, we did not keep track of the R version with which we ran the scripts here, so you** 
**might not get exactly the same results that you see here if you re-run the scripts with the R**
**function `sn::st.mple`.**   

**NOTE 2: If you re-run the scripts with the results obtained with `MCMCtree` (see [here](/01_SeqBayes_S1/02_MCMCtree)**
**to find the link to download the data) with the same R version we did (we are unsure about which one this was,**
**but it was one of the released R versions between 2018 and 2020), you will get the results we provide in**
**this directory when you use the data in [`00_main_tree_T2/01_MCMCtree_posterior_old`](https://www.dropbox.com/s/qxsgfe0gbwxro9p/SeqBayesS1_MCMCtree_mainT2_posterior_old.zip?dl=0).**
**Note that removing 11 genes out of the final filtered genes in the alignment did not have an impact on the**
**estimated posterior divergence times as we show** 
**[here](https://github.com/sabifo4/mammals_dating/blob/main/01_SeqBayes_S1/02_MCMCtree/00_MCMCtree_analyses/00_main_tree_T2/plot_oldtimesVSnewtimes/00_Check_oldpostVSnewpost-I.pdf)** 
**and [here](https://github.com/sabifo4/mammals_dating/blob/main/01_SeqBayes_S1/02_MCMCtree/00_MCMCtree_analyses/00_main_tree_T2/plot_oldtimesVSnewtimes/00_Check_oldpostVSnewpost-II.pdf)**
**(see also further comments [here](/01_SeqBayes_S1/02_MCMCtree)).**
**Consequently, we did not repeat all the steps described above to fit the ST distributions again to the alignment without**
**the 11 genes because we would not get any significant differences. You can also find in directory**
**[`00_fitST_new`](/01_SeqBayes_S1/03_Fit_ST_to_posteriors/00_fitST_new)**
**the results with the generated fitted ST distributions to the posterior densities obtained when using the alignment**
**without these 11 genes -- you can see that there are no significant differences [here](https://github.com/sabifo4/mammals_dating/blob/main/01_SeqBayes_S1/03_Fit_ST_to_posteriors/02_MCMCtree_prior_ST_new/plots/Compare_OLDvsNEW_ST.png)**
**and in individual plots [here](/01_SeqBayes_S1/03_Fit_ST_to_posteriors/02_MCMCtree_prior_ST_new/plots/ST_comparison)**
**(found within [`02_MCMCtree_prior_ST_new/plots`](/01_SeqBayes_S1/03_Fit_ST_to_posteriors/02_MCMCtree_prior_ST_new/plots) directory)**
**between the ST distributions.**

---
The plots show that the ST distributions are sensible (see image below), and so they can be used to continue with the sequential
Bayesian dating analysis. The ST distributions estimated with the `sn::st.mple` R function are plotted in red, those estimated when
sampling from the prior with `MCMCtree` in light green, and when sampling from the posterior with `MCMCtree` in black:

<p align="center">
  <img width="500" height="500" src="https://github.com/sabifo4/mammals_dating/blob/main/figs/FigS6.png">
</p>

>>**Fig S6. Assessment of the skew-_t_ (ST) distributions fitted to the internal nodes of the 72-species mammal tree.**