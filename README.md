# Dense Phylogenomic Dating Provides a High-Resolution Timeline of Mammal Evolution

## What will you find here?
In this repository, you can find a description of the several steps we have carried out to infer the divergence times of a phylogeny with 
4,705 mammal taxa using a sequential Bayesian approach. In a nutshell, this approach consists of using the inferred posterior divergence times on
a first data set (in this study, this is an alignment with 15,268 filtered genes from ENSEMBL ordered from slow- to fast-evolving and partitioned into
four block 4, for 72 mammal taxa) as calibration priors for a subsequent dating analyses with a much larger data set (in this study, this is a
concatenated alignment of 168 nuclear, 12 mitochondrial protein-coding, and 2 mitochondrial non-coding genes for 4,705 mammal taxa).
The advantage of this method is how much it can reduce the amount of computational time needed to date such a large phylogeny.

## How is the content structured?
All the details about our study is structured as it follows:   

   1. [**00_Data_collection**](https://github.com/sabifo4/mammals_dating/tree/main/00_Data_collection)   
   Every project starts with the data collection. You can follow the link above to read all the details about this step.   
   2. [**01_SeqBayes_S1**](https://github.com/sabifo4/mammals_dating/tree/main/01_SeqBayes_S1)      
   Once the data are gathered, they need to be processed and filtered before proceeding with the Bayesian dating approach. Therefore,
   this directory contains several subdirectories in which all these steps are detailed:   
      * [00_Gene_filtering](https://github.com/sabifo4/mammals_dating/tree/main/01_SeqBayes_S1/00_Gene_filtering)   
      Here, you will find all the analyses we carried out during the gene filtering step and how we obtained the partitioned alignment
	  with 72 taxa.   
      * [01_BASEML](https://github.com/sabifo4/mammals_dating/tree/main/01_SeqBayes_S1/01_BASEML)   
      This link will redirect you to the description of the steps taken to calculate the Hessian and the gradient of the partitioned alignment.   
      * [02_MCMCtree](https://github.com/sabifo4/mammals_dating/tree/main/01_SeqBayes_S1/02_MCMCtree)   
      Here, you can find the results for the Bayesian dating approach using the approximate likelihood as implemented in `MCMCtree`. Besides, you will also find the analyses to evaluate chain convergence.   
      * [03_Fit_ST_to_posteriors](https://github.com/sabifo4/mammals_dating/tree/main/01_SeqBayes_S1/03_Fit_ST_to_posteriors)   
      After averaging the posterior time estimates across the samples collected throughout the different MCMCs, we can fit skew-_t_ (ST) distributions to each of the nodes of the mammal phylogeny. You can find here the details of this approach.   
   3. [**02_SeqBayes_S2**](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2)   
   Once we have the ST distributions fitted to each node, we can proceed to use them as prior calibrations when dating the larger mammal data set
   with ~5,000 taxa (data set 2, number of taxa prior to data curation). The tasks are also divided in different steps as it was done in the first step:   
      * [00_Data_filtering](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering)   
      This directory contains a summary of the procedure followed to (i) curate the alignments for each data subset when working
	  with the second data set and (ii) generate the corresponding phylogenies. For each data subset, a directory has been generated with a 
	  `README.md` file in which a step-by-step guideline details all the filtering steps followed. 
      * [01_BASEML](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/01_BASEML)   
      Once the alignments and trees are ready, we can then compute the Hessian and the gradient for the subsequent Bayesian inference of divergence times.
	  Click the link above to be redirected to the steps and results for this step.    
      * [02_MCMCtree](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/02_MCMCtree)   
      Here you will find the final results regarding the inferred posterior times for the mammal phylogeny with 4,705 taxa.   
   4. [**calibrations**](https://github.com/sabifo4/mammals_dating/tree/main/calibrations)   
   This directory contains an excel file with a summary of the calibration priors used in both Bayesian dating analyses.   
   5. [**figs**](https://github.com/sabifo4/mammals_dating/tree/main/figs)   
   This directory contains figures that are embedded in some of the `README.md` files in this GitHub repository.   
   6. [**src**](https://github.com/sabifo4/mammals_dating/tree/main/src)   
   This directory contains the software needed to carry out the tasks described above.   


---

We hope you can find this useful! :)