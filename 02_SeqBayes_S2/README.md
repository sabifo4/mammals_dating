# 1. Quick summary of tasks accomplished in step 1

## 1.1 Generating alignment with 72 species and 15K loci (genomic data)
[`Here`](/01_SeqBayes_S1/00_Gene_filtering),
you will find the steps followed to obtain 
the final partitioned genomic alignment for the 72 mammal taxa used in the first step of the sequential 
Bayesian dating analysis:   

   1) Filtering 15,904 genes shared acrossed the 72 mammal taxa of the phylogeny to 15,268 genes.   
   2) Ordering the filtered genes from slow- to fast-evolving.   
   3) Generating the 4-partitions alignment with the ordered genes.   
   
## 1.2. Running `BASEML` and `MCMCtree`
Before starting with the Bayesian clock-dating analysis, we computed the Hessian and the gradient for each 
of the four partitions with `BASEML` and under each of the seven tree hypotheses
so we could later use the approximate likelihood calculation implemented in `MCMCtree` 
(see [dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613)). 
We remind you of the proceudre we followed to run `BASEML`:   

   1) Prepare one folder for each of the partitions in which you have a (1) file with the tree topology, (2) a file 
      with the alignment of the corresponding partition, and (3) the control file to run `MCMCtree`.
	  As we want to run the Bayesian clock-dating analysis under 7 different tree hypotheses, you should have 
	  a file architecture like the following in your working directory (one for each tree hypothesis):
	  
	  ```
	  TreeX
	     |
	     |-p01 
	     |  |- alignment_part1.aln
	     |  |- control_file.ctl 
	     |  |- tree.tree
	     |
	     |-p02 
	     |  |- alignment_part2.aln
	     |  |- control_file.ctl 
	     |  |- tree.tree
	     |
	     |-p03
	     |  |- alignment_part3.aln
	     |  |- control_file.ctl 
	     |  |- tree.tree
	     |
	     |-p04 
	        |- alignment_part4.aln
	        |- control_file.ctl 
	        |- tree.tree
	  ```
   2) Now, make sure that you set the option `usedata = 3` in the control files.   
   3) You can execute `MCMCtree` within each directory or use bash script as in a pipeline. Note, however, that we will not let 
      `MCMCtree` run until it finishes its job. We can kill each run once the files `tmp0001*` have been 
	  generated. We will be interested in keeping `tmp0001.ctl`, `tmp0001.trees`, and `tmp0001.txt` files; which are the control file 
	  the tree file, and the alignment file, respectively. These files are needed to run `BASEML`. The 
	  following file architecture should be the one you should then have for the analyses carried out under each tree hypothesis:   
	  
	  ```
	  TreeX
	     |
	     |-p01 
	     |  |- tmp0001.trees
	     |  |- tmp0001.ctl 
	     |  |- tmp0001.txt
	     |
	     |-p02 
	     |  |- tmp0001.trees
	     |  |- tmp0001.ctl 
	     |  |- tmp0001.txt
	     |
	     |-p03
	     |  |- tmp0001.trees
	     |  |- tmp0001.ctl 
	     |  |- tmp0001.txt
	     |
	     |-p04 
	        |- tmp0001.trees
	        |- tmp0001.ctl 
	        |- tmp0001.txt
	  ```
	  
	  Now, before running `BASEML` in each subdirectory, you need to modify the `tmp0001.ctl` file and 
	  set `method = 1`. Make sure all your `tmp0001.ctl` files have the following content:   
	  ```
		seqfile = tmp0001.txt
		treefile = tmp0001.trees
		outfile = tmp0001.out
		noisy = 3
		model = 4
		fix_alpha = 0
		alpha = 0.5
		ncatG = 5
		Small_Diff = 0.1e-6
		getSE = 2
		method = 1   
	  ```
	  
	  After that, do not change any file name and just execute `BASEML` (we recommend using a pipeline to speed up the analysis).   
	  
   4) Once `BASEML` has finished for each partition under each tree hypothesis, you want to save the Hessian and 
      the gradient computed, which can be found in one of the output files: `rst2`. You will need 
      to generate an `in.BV` file for each tree hypothesis with the content of each of the `rst2` files output for each partition.
      **IMPORTANT**:  The alignment that will be used in `MCMCtree` will have the partitions
      ordered as `p01>p02>p03>p04` (i.e., genes ordered from slow- to fast-evolving in these partitions).  
      Therefore, make sure that you concatenate the `rst2` files in the correct order, that is, `p01>p02>p03>p04`. We want that the Hessian and the gradients are properly matched 
      to the corresponding partition when `MCMCtree` uses the order found in both the `in.BV` and the alignment files.
      For instance, you could use a bash code like the following:
	  
	  ```sh
	  for i in `seq 1 4`
	  do
	  cat p0$i/rst2 >> in.BV 
	  printf "\n\n" >> in.BV 
	  done
	  ```
	  
	  Then, you should place the `in.BV` file in every directory where `MCMCtree` will run under the approximate 
	  likelihood and the autocorrelated-rates model (GBM). The directory for each tree hypothesis should look like the following:   
	  
	  ```
	  TreeX
	     |
	     |- GBM_model
	     	 |- mcmctree.ctl 
	     	 |- alignment_4partitions.aln 
	     	 |- in.BV 
	     	 |- 72sp_calibrated.tree  
	     
	  ```
	  
	  The `alignment_4partitions.aln` file should have the alignments with the 4 partitions in the correct 
	  order (`p01>p02>p03>p04`) and the header for each partition should be in phylip format.
	  The `72sp_calibrated.tree` should contain the calibrated tree topology (Newick format) evaluated at the moment (which is fixed 
	  in `MCMCtree`). 
	  The `mcmctree.ctl` should look like the following:   
	  ```
	     seed = -1
          seqfile = alignment_4partitions.aln
         treefile = 72sp_calibrated.tree
         mcmcfile = mcmc.txt
          outfile = out.txt
            ndata = 4
          seqtype = 0    * 0: nucleotides; 1:codons; 2:AAs
          usedata = 2    * 0: no data (prior); 1:exact likelihood;
                         * 2:approximate likelihood; 3:out.BV (in.BV)
            clock = 3    * 1: global clock; 2: independent rates; 3: correlated rates
            model = 4    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
            alpha = 0.5  * alpha for gamma rates at sites
             ncatG = 5    * No. categories in discrete gamma
         cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?
           BDparas = 1 1 0.1    * birth, death, sampling
       rgene_gamma = 2 40   * gammaDir prior for rate for genes
      sigma2_gamma = 1 10    *s gammaDir prior for sigma^2     (for clock=2 or 3)
             print = 1   * 0: no mcmc sample; 1: everything except branch rates 2: everything
            burnin = 100000
            sampfreq = 1000 
           nsample = 20000
	   ```   
	   
	   Note that, in order to be able to gather enough samples and ensuring convergence, `MCMCtree` was run twice.
	   The `seqfile` and `treefile` used here match the example given above. Make sure that you use the names
	   you have for your alignment and tree files. Note that the `burnin`, `sampfreq`, and `nsample` options 
	   vary when running `MCMCtree` when sampling from the prior and when sampling from the posterior. You can 
	   find the analyses under the prior 
	   [here](https://www.dropbox.com/s/09u8l81dgw166do/SeqBayesS1_MCMCtree_6treehyp_MCMCtree.zip?dl=0),
	   under the posterior with the main tree hypothesis 
	   [here](https://www.dropbox.com/s/kh73rbu9cnxts3r/SeqBayesS1_MCMCtree_mainT2_posterior.zip?dl=0),
	   and under the posterior with the rest of tree hypotheses 
	   [here](https://www.dropbox.com/s/09u8l81dgw166do/SeqBayesS1_MCMCtree_6treehyp_MCMCtree.zip?dl=0).
	  
## 1.3. Fitting skew-_t_ distributions to resulting posterior densities
Once `MCMCtree` finished, we used the results obtained under the 
main tree hypothesis 
and fitted a skew-_t_ (ST) distribution to each of the estimated posterior
distributions for each node. The whole procedure is detailed
[here](/01_SeqBayes_S1/03_Fit_ST_to_posteriors)).

# 2. Second part of the sequential Bayesian dating approach 

## 2.1. Outline 
The mammal phylogeny that contains over 5K mammal taxa and 182 genes is too big to be analysed as in 
a unique alignment. Therefore, we had to work with separate data subsets with less taxa to reduce
computing time so we could run the clock-dating analysis in `MCMCtree`. 

For that purpose, taxa were separated according to the following taxonomical
groups: Afrotheria, Xenarthra, Marsupialia, Euarchonta, Lagomorpha, Laurasiatheria, Rodentia, and Monotremata.
The number of taxa within Laurasiatheria and Rodentia,
however, was still too large to run a Bayesian clock-dating analysis in `MCMCtree` within a reasonable
amount of time despite using the approximate likelihood. Therefore, these data sets were further divided
into other data subsets. For Laurasiatheria, we generated data subsets with Artiodactyla (which you might find 
as "laurasiatheria cetartiodactyla" in other tutorials and data sets provided in this GitHub repository), Chiroptera
(which you might find 
as "laurasiatheria chiroptera" in other tutorials and data sets provided in this GitHub repository),
and a third data subset with the rest of the taxa belonging to Laurasiatheria that were not included in the 
previous data subsets. Note that, however, after running the analyses with 
the Chiroptera data subset, we realised that it had to be further subset: we generated "chiroptera subtree 1" and 
"chiroptera subtree 2" data subsets. In total, we had four data subsets for Laurasiatheria.
When it comes to the Rodentia data set, we first generated 3 data subsets: Sciuridae and related (which you might find 
as "rodentia squirrel" or "squirrel" in other tutorials and data sets provided in this GitHub repository),
Ctenohystrica, and a third data subset with the rest of rodents. As it happened with the Chiroptera data
subset, the data subset with the rest of rodentia had too many taxa for `MCMCtree` to run. Therefore, we 
further subset it into two: "rodentia subtree 1" and "rodentia subtree 2" data subsets. In total, we had 
4 data subsets for Rodentia. The directories where all the steps followed to curate each of these data subsets
are explained (see 
[here](/02_SeqBayes_S2/00_Data_filtering/00_data_curation))
are the following:   

   * [Afrotheria](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/afrotheria)   
   * [Chiroptera](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/chiroptera)   
   * [Euarchonta](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/euarchonta)   
   * [Lagomorpha](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/lagomorpha)   
   * [Laurasiatheria_cetartiodactyla](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/laurasiatheria_cetartiodactyla)   
   * [Laurasiatheria_therest](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/laurasiatheria_therest)   
   * [Marsupialia](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/marsupialia)   
   * [Rodentia_ctenohystrica](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_ctenohystrica)   
   * [Rodentia_squirrel](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_squirrel)   
   * [Rodentia_therest](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_therest)   
   * [Xenarthra](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/xenarthra)   

Note that all these data subsets (the steps followed to further subset Laurasiatheria and Rodentia data subsets 
as mentioned above is included in the corresponding directories listed above) had the corresponding inferred 
phylogeny rooted using Monotremata as an outgroup, which includes *Tachyglossus aculeatus*, *Zaglossus bruijni*, and *Ornithorhynchus anatinus*.

The 182 genes available for the taxa in these data subsets are either mitochondrial or nuclear data.
We generated the following alignments for each data subset:   

   1) `mt_3cp`: mitochondrial coding genes, only the concatenated third codon positions of the 
                genes found in the taxa included in the corresponding data subset.   
   2) `mt_12cp`: mitochondrial coding genes, only the concatenated first and second codon positions of the 
                genes found in the taxa included in the corresponding data subset.   
   3) `nt_3cp`: nuclear coding data, only the concatenated third codon positions of the 
                genes found in the taxa included in the corresponding data subset.   
   4) `nt_12cp`: nuclear coding data, only the concatenated first and second codon positions of the 
                genes found in the taxa included in the corresponding subtree.   
   5) `mt_rna`: mitochondrial RNA data, concatenated genes found in the taxa included in
                the corresponding data subset.   
   6) `concatenated alignment`: all the previous 5 partitions concatenated, one after the other, for the taxa included in
                the corresponding data subset. This data subset is not used in any subsequent analysis.
   7) **`5-partitions alignment`**: all the first 5 partitions separated in individual blocks, one block per partition.
                This is the alignment that will be used in `MCMCtree` to infer 
                the 4.7K mammal timetree.

The tutorials available in the directories for each data subset listed before explain how to generate them and provide you 
with the links to download the alignments for each data subset. 

In addition, for each data subset, we had the best-scoring maximum-likelihood (ML) tree estimated by
[`RAxML`](https://github.com/stamatak/standard-RAxML).
This ML tree was then pruned according to the filtering steps mentioned in the tutorials which link is 
provided above for each data subset. The pruned ML tree was then calibrated by first manually
placing tags in the nodes to be calibrated and then using and R script to replace these tags with 
the corresponding calibrations (see the `filter_tree/00_Filter_trees` content inside the directories for 
each data subset for more information).

[Here](/02_SeqBayes_S2/00_Data_filtering/01_alignments) 
you will find the instructions to run a perl script that will count the missing data for each 
individual partition of each data subset. 

## 2.2. Run `BASEML` to estimate the Hessian and the gradient
Once the data subsets are filtered and the corresponding calibrated phylogenies are obtained,
we can run `BASEML` to estimate the Hessian and the gradient needed to run `MCMCtree` with the 
approximate likelihood. You may check the tutorial used in the first step
[here](/01_SeqBayes_S1/01_BASEML/02_Hessian)
to remember how we carry out this analysis. 

[`Here`](/02_SeqBayes_S2/01_BASEML),
you will find a description of how we ran `BASEML` with the data subsets for this second step 
of the sequential Bayesian dating approach. We also provide you with the links to download the data 
you will need to reproduce our results as well as the results we obtained. 

# 2.3. Run `MCMCtree` to infer species divergence times for each data subset
Once the Hessian and the gradient have been estimated for each data subset, we can run `MCMCtree`
with the approximate likelihood 
calculation ([dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613)). 

Before running `MCMCtree`, however, we ran safety checks to make sure that the calibrations that we 
were using were not in conflict. For each of the directories with the information to filter each 
data subset that you can find 
[here](/02_SeqBayes_S2/00_Data_filtering/00_data_curation),
you will find a description of these checks inside the `filter_tree` directory. The R script and the 
output plots will be inside `filter_tree/01_Check_conflict` dirctories for each data subset.

Once all the checks were finished, we ran `MCMCtree` with each data subset. You can download the results 
obtained [here](https://www.dropbox.com/s/1vjkggr4ujrnfha/SeqBayesS2_MCMCtree.zip?dl=0),
while you can find 
[here](/02_SeqBayes_S2/02_MCMCtree) 
the output files that we generated for each data subset to assess chain convergence and summarise the 
estimated divergence times.
