# A Species-Level Timeline of Mammal Evolution Integrating Phylogenomic Data

## What will you find here?
In this repository, you can find a description of the several steps we have carried out to infer the divergence times of a phylogeny with 
4,705 mammal taxa using our new Bayesian sequential-subtree dating approach. In a nutshell, this approach consists of using the inferred posterior divergence times on
a first dataset (in this study, this is a 72-taxon alignment generated with 15,268 filtered mammal genomes from ENSEMBL ordered from slow- to fast-evolving and partitioned into
four blocks) as calibration priors for a subsequent dating analysis with a much larger dataset (in this study, this is a
concatenated alignment of 168 nuclear, 12 mitochondrial protein-coding, and 2 mitochondrial non-coding genes for 4,705 mammal taxa).
The advantage of this method is how much it can reduce the amount of computational time needed to date such a large phylogeny.

## How is the content structured?
We have divided this repository into the following sections:   

   1. [**00_Data_collection**](00_Data_collection)   
   Every project starts with the data collection. You can follow the link above to read all the details about this step.   
   2. [**01_SeqBayes_S1**](01_SeqBayes_S1)      
   Once the data are gathered, they need to be processed and filtered before proceeding with Bayesian inference. Therefore,
   this directory contains several subdirectories in which all these steps are detailed:   
      * [00_Gene_filtering](01_SeqBayes_S1/00_Gene_filtering)   
      Here, you will find all the analyses we carried out to filter the data and the steps we followed to then obtain the partitioned alignment
	  with 72 taxa.   
      * [01_BASEML](01_SeqBayes_S1/01_BASEML)   
      Here, you will find the step-by-step tutorial we followed to calculate the Hessian and the gradient of the partitioned alignment.   
      * [02_MCMCtree](01_SeqBayes_S1/02_MCMCtree)   
      Here, you can find the results of the Bayesian dating approach using the approximate likelihood as implemented in `MCMCtree`. Besides, you will also find the analyses to evaluate chain convergence.   
      * [03_Fit_ST_to_posteriors](01_SeqBayes_S1/03_Fit_ST_to_posteriors)   
      After averaging over the posterior time estimates sampled during the independent MCMCs, we can fit skew-_t_ (ST) distributions to each of the nodes of the 72-taxon phylogeny. Here, you can find all the details about this approach.   
   3. [**02_SeqBayes_S2**](02_SeqBayes_S2)   
   Once we have the ST distributions fitted to each node, we can proceed to use them as prior calibrations when dating the larger mammal dataset
   with ~5,000 taxa (dataset 2, number of taxa prior to data curation). The tasks are also divided in different steps as it was done in the first step:   
      * [00_Data_filtering](02_SeqBayes_S2/00_Data_filtering)   
      This directory contains a summary of the procedure followed to (i) curate the alignments for each data subset when working
	  with the second dataset and (ii) generate the corresponding phylogenies. For each data subset, a directory has been generated with a 
	  `README.md` file in which a step-by-step guideline details all the filtering steps followed. 
      * [01_BASEML](02_SeqBayes_S2/01_BASEML)   
      Once the alignments and trees are ready, we can compute the Hessian and the gradient for the subsequent Bayesian inference of divergence times.
	  Click the link above to be redirected to the steps and results of this step.    
      * [02_MCMCtree_updcrh](02_SeqBayes_S2/02_MCMCtree_updcrh)   
      Here, you will find the final results regarding the inferred posterior times for the mammal phylogeny with 4,705 taxa.   
      * [03_Generate_final_mammal_tree](02_SeqBayes_S2/03_Generate_final_mammal_tree)   
      Here, you will find the final mammal Tree of Life that we generated at the end of our Bayesian sequential-subtree dating approach.   
   4. [**03_Extra_analyses**](03_Extra_analyses)   
   This directory contains a description of extra analyses that were carried out after the main analyses were finished. The
   main purpose was to compare the divergence times inferred with the 72-taxon alignment (dataset 1) to those estimated 
   when using (i) 182 loci, (ii) only nuclear (1st+2nd CPs) loci, and (iii) only mitochondrial
   (1st+ 2nd CPs) loci for the same 72 mammal taxa (genes extracted from dataset 2 for these 72 mammal taxa).   
   5. [**calibrations**](calibrations)   
   This directory contains a pdf file with the justification for the calibrations used.   
   6. [**figs**](figs)   
   This directory contains figures that are embedded in some of the `README.md` files in this GitHub repository.   
   7. [**src**](src)   
   This directory contains the software needed to carry out the tasks described above.   

## Do you have a compiled version of your data?
Yes, we do! We have uploaded a zip file with the data that you need to reproduce our analyses in 
[`FigShare`](https://figshare.com/), which you can access at DOI 10.6084/m9.figshare.14885691.

In this zip file you will find several directories with the following content:   

   * `aln`: here you will find the alignments used in both steps of our Bayesian
   sequential-subtree approach. Directory `aln/00_step_01/` contains the alignments for the
   72 genomes in phylip format. Directory `aln/01_step_02/` contains the uncompressed, raw subtree
   alignments for the 4,705 taxa, with missing species represented as sequences of gaps.
   Directory `aln/01_step_02_patterns/` contains the same alignments after processing, that is,
   with missing species removed and with the alignments compressed into site patterns
   (see `MCMCtree` and `PAML` documentation for alignment formats). The processed alignments
   are the ones used by MCMCtree to calculate the likelihood during estimation of gradient and
   Hessian. Note each alignment file contains several alignment blocks, with each block corresponding
   to an alignment partition. If you load an alignment file into an alignment editor, make sure
   your editor allows you to see all partition blocks and not just the first one.   
   * `clocktest`: here, you will find the full results for the Bayesian selection of relaxed-clock
   model.   
   * `inBV`: here you will find the `in.BV` files generated by `BASEML`, which contain the estimated gradient
   and Hessian for each alignment
   partition, which are required to estimate the divergence times under the approximate likelihood method in
   `MCMCtree`. This directory is subdivided into `step01` and `step02` subdirectories corresponding to
   the data for the 72 genomes and the 4,705-taxon subtrees, respectively.
   * `paleodb`: here you will find a `csv` file with the details about the mining of mammal genera
   from [PaleoDB](https://paleobiodb.org/).   
   * `timetrees`: here you will find the timetrees estimated with `MCMCtree` in Nexus format.
   These are suitable for plotting with `FigTree`. For the 4,705 taxa, both the separate subtrees and the fully
   stitched 4,705-taxon tree are provided.   
   * `trees`: here you will find the calibrated tree topologies (input tree files needed by `MCMCtree`) 
   in Newick format for each of the dataset (both from the first and the second step of our Bayesian approach).
   Subdirectories `step01` and `step02 `contain the corresponding trees for the 72 genomes and 4,705 taxa, respectively.
   These trees are required, together with the alignment and `in.BV` files to estimate the divergence times with `MCMCtree.   

In addition, please note that, in this GitHub repository, we have `README.md` files where we explain in detail 
each analysis we carried out during this study in a step-by-step tutorial manner.
Apart from providing you with code snippets and links to the 
scripts and source code we used, we also include in these `README.md` files the links to zip files
that you can use to reproduce each analysis (as well as to check if you have gotten the same results we did!). 

---

We hope you can find this repository useful. Happy Bayesian inference! :)
