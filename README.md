# A Species-Level Timeline of Mammal Evolution Integrating Phylogenomic Data

## What will you find here?
In this repository, you can find a description of the several steps we have carried out to infer the divergence times of a phylogeny with 
4,705 mammal taxa using our new Bayesian sequential-subtree dating approach. In a nutshell, this approach consists of using the inferred posterior divergence times on
a first data set (in this study, this is an alignment with 15,268 filtered genes from ENSEMBL ordered from slow- to fast-evolving and partitioned into
four block 4, for 72 mammal taxa) as calibration priors for a subsequent dating analyses with a much larger data set (in this study, this is a
concatenated alignment of 168 nuclear, 12 mitochondrial protein-coding, and 2 mitochondrial non-coding genes for 4,705 mammal taxa).
The advantage of this method is how much it can reduce the amount of computational time needed to date such a large phylogeny.

## How is the content structured?
All the details about our study is structured as it follows:   

   1. [**00_Data_collection**](00_Data_collection)   
   Every project starts with the data collection. You can follow the link above to read all the details about this step.   
   2. [**01_SeqBayes_S1**](01_SeqBayes_S1)      
   Once the data are gathered, they need to be processed and filtered before proceeding with the Bayesian dating approach. Therefore,
   this directory contains several subdirectories in which all these steps are detailed:   
      * [00_Gene_filtering](01_SeqBayes_S1/00_Gene_filtering)   
      Here, you will find all the analyses we carried out during the gene filtering step and how we obtained the partitioned alignment
	  with 72 taxa.   
      * [01_BASEML](01_SeqBayes_S1/01_BASEML)   
      This link will redirect you to the description of the steps taken to calculate the Hessian and the gradient of the partitioned alignment.   
      * [02_MCMCtree](01_SeqBayes_S1/02_MCMCtree)   
      Here, you can find the results for the Bayesian dating approach using the approximate likelihood as implemented in `MCMCtree`. Besides, you will also find the analyses to evaluate chain convergence.   
      * [03_Fit_ST_to_posteriors](01_SeqBayes_S1/03_Fit_ST_to_posteriors)   
      After averaging the posterior time estimates across the samples collected throughout the different MCMCs, we can fit skew-_t_ (ST) distributions to each of the nodes of the mammal phylogeny. You can find here the details of this approach.   
   3. [**02_SeqBayes_S2**](02_SeqBayes_S2)   
   Once we have the ST distributions fitted to each node, we can proceed to use them as prior calibrations when dating the larger mammal data set
   with ~5,000 taxa (data set 2, number of taxa prior to data curation). The tasks are also divided in different steps as it was done in the first step:   
      * [00_Data_filtering](02_SeqBayes_S2/00_Data_filtering)   
      This directory contains a summary of the procedure followed to (i) curate the alignments for each data subset when working
	  with the second data set and (ii) generate the corresponding phylogenies. For each data subset, a directory has been generated with a 
	  `README.md` file in which a step-by-step guideline details all the filtering steps followed. 
      * [01_BASEML](02_SeqBayes_S2/01_BASEML)   
      Once the alignments and trees are ready, we can then compute the Hessian and the gradient for the subsequent Bayesian inference of divergence times.
	  Click the link above to be redirected to the steps and results for this step.    
      * [02_MCMCtree](02_SeqBayes_S2/02_MCMCtree)   
      Here you will find the final results regarding the inferred posterior times for the mammal phylogeny with 4,705 taxa.   
      * [03_Generate_final_mammal_tree](02_SeqBayes_S2/03_Generate_final_mammal_tree)   
      Here you will find the final mammal Tree of Life that we generated at the end of our Bayesian sequential-subtree dating approach.   
   4. [**03_Extra_analyses**](03_Extra_analyses)   
   This directory contains a description of extra analyses that were carried out after the main analyses were finished. The
   main purpose was to compare the divergence times inferred with the 72-sp alignment (data 1) to those estimated for these same
   72 taxa when using (i) the 182 loci for these taxa, (ii) only the nuclear (1st+2nd CPs) for these taxa, and (iii) only mitochondrial
   (1st+ 2nd CPs) for these taxa (everything extracted from data 2).   
   5. [**calibrations**](calibrations)   
   This directory contains an excel file with a summary of the calibration priors used in both Bayesian dating analyses and
   a pdf file with the justification for the calibrations used.   
   6. [**figs**](figs)   
   This directory contains figures that are embedded in some of the `README.md` files in this GitHub repository.   
   7. [**src**](src)   
   This directory contains the software needed to carry out the tasks described above.   

## Do you have a compiled version of your data?
Yes, we do! We have uploaded a zip file with the data that you need to reproduce our analyses in 
[`FigShare`](https://figshare.com/), which you can access
[here](https://figshare.com/s/4718bffeae304f754350).

In this zip file you will find several directories with the following content:   

   * `aln`: here you will find the alignments used in both steps of our Bayesian
   sequential-subtree approach.   
   * `clocktest`: here you will find the results of the Bayesian model selection analysis of rate model.   
   * `inBV`: here you will find the `in.BV` files generated by `BASEML` for each data alignment (both from 
   the first and the second step). These files are needed to run `MCMCtree` with the approximate likelihood.   
   * `paleodb`: here you will find a `csv` file with the details about the mining of mammal genera
   from [PaleoDB](https://paleobiodb.org/).   
   * `timetrees`: here you will find the timetrees estimated by `MCMCtree` when using the alignments for
   each data set (both from the first and the second step of our Bayesian approach). 
   * `trees`: here you will find the calibrated tree topologies (input tree files needed by `MCMCtree`) for 
   each of the data set (both from the first and the second step of our Bayesian approach).   

In addition, please note that, in this GitHub repository, we have `README.md` files where we explain in detail 
each analysis we carried out during this study in a step-by-step tutorial manner.
Apart from providing you with code snippets and links to the 
scripts and source code we used, we also include in these `README.md` files the links to zip files
that contain the data that you will need to reproduce each analysis as well as the results we obtained. 

---

We hope you can find this repository useful. Happy Bayesian inference! :)
