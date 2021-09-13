# 1. Quick summary overview of tasks accomplished so far 

## 1.1 Generating alignment with 72 sp 
In the [`01_SeqBayes_S1/00_Gene_filtering`](https://github.com/sabifo4/mammals/tree/master/01_SeqBayes_S1/00_Gene_filtering) directory,
you will find the `README.md` where the whole process followed to obtain 
the final partitioned alignment for the 72 mammal taxa used in the first step of the Bayesian 
sequential dating analysis:   

   1) Filtering 15,904 genes shared acrossed the 72 mammal taxa of the phylogeny down to 15,431 genes.   
   2) Ordering the filtered genes from slow- to fast-evolving.   
   3) Concatenating the ordered genes and dividing the resulting alignment into 4 partitions.   
   
## 1.2. Running `BASEML` and `MCMCtree` with 72sp
Before starting with the Bayesian dating analysis, we computed the Hessian and the gradient for each 
of the four partitions with `BASEML` so we could use the approximate likelihood calculation implemented in `MCMCtree` 
([dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613)). 
Nevertheless, the format in which the alignment needs to be to be read by `BASEML` had to be first generated.
We remind you of the trick we did to easily obtain this:   

   1) Prepare one folder for each of the partitions in which you have a (1) file with the tree topology, (2) a file 
      with the alignment of the corresponding partition, and (3) the control file to run `MCMCtree`.
	  As we want to run the Bayesian dating analysis under 7 different tree topology hypotheses, you should have 
	  something like this in your working directory but for each of the tree topologies (i.e., the file "tree.tree" 
	  will change as it will contain the topology evaluated in each case):
	  
	  ```
	  Tree1
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
   3) You can run `MCMCtree` within each directory or work with a pipeline. Note, however, that we will not let 
      `MCMCtree` run until it finishes the job. We can kill each run once the files `tmp0001*` have been 
	  generated. We will be interested in keeping `tmp0001.ctl`, `tmp0001.trees`, and `tmp0001.txt`; which are the control file 
	  the tree file, and the alignment file, respectively; needed to run `BASEML`. This means that we will want to have now the 
	  following file architecture for each tree topology:   
	  
	  ```
	  Tree1
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
	  
	  Nevertheless, before running `BASEML` in each subdirectory, you need to modify the `tmp0001.ctl` file and 
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
	  
	  Do not change any file name and then just run `BASEML` (we recommend using a pipeline to speed up things).   
	  
   4) Once `BASEML` has finished for each partition under each tree topology, you want to save the Hessian and 
      the gradient computed, which can be found in one of the output files: `rst2`. You will need 
      to generate an `in.BV` file with the content of each of the `rst2` files corresponding to each partition.
      **IMPORTANT**: Make sure that you concatenate the `rst2` files in the correct order of the partitions, that 
      is, `p01>p02>p03>p04`. This is because the alignment that will be used in `MCMCtree` will have the partitions 
      written in this order in the alignment file. Therefore, we want that the Hessian and the gradients are matched 
      to the corresponding partition when `MCMCtree` uses the order found in both the `in.BV` and the alignment files.
      For instance, you could use a bash code like the following:
	  
	  ```
	  for i in `seq 1 4`
	  do
	  cat p0$i/rst2 >> in.BV 
	  printf "\n\n" >> in.BV 
	  done
	  ```
	  
	  Then, you should place the `in.BV` file in every directory where `MCMCtree` will be run under the approximate 
	  likelihood. The directory should look like this for each tree topology tested:   
	  
	  ```
	  Tree1
	     |
	     |- GBM_model
	     |	 |- mcmctree.ctl 
	     |	 |- alignment_4partitions.aln 
	     |	 |- in.BV 
	     |	 |- 72sp_calibrated.tree  
	     |
	     |- ILN_model
	      	 |- mcmctree.ctl 
	      	 |- alignment_4partitions.aln 
	      	 |- in.BV 
	      	 |- 72sp_calibrated.tree  
	  ```
	  
	  The `alignment_4partitions.aln` should have the alignments with the 4 partitions in the correct 
	  order (`p01>p02>p03>p04`) and the header for each partition should be in phylip format.
	  The `72sp_calibrated.tree` should have the tree topology evaluated at the moment (which is fixed 
	  in `MCMCtree`) with all the node calibrations used included. Note that we ran it under both 
	  the autocorrelated-rates relaxed-clock model (`GBM`) and the log-normal independent-rates relaxed-clock
	  model (`ILN`).
	  The `mcmctree.ctl` should look like the following inside the `GBM_model` directory, just changes
	  `clock = 2` when running it in the `ILN_model` directory:   
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
            mpfreq = 1000 
           nsample = 20000
	   ```   
	   Note that, in order to be able to gather enough samples and ensuring convergence, `MCMCtree` was run twice using these
	   parameters. The `seqfile` and `treefile` used here match the example given above. Make sure that you use the names
	   you have for your alignment and tree files.
	  
## 1.3. Fitting skew-_t_ distributions to resulting posterior time estimates
Once `MCMCtree` finished, we used the estimated posterior parameters obtained under the 
second tree hypothesis (i.e., [Atlantogenata_tarver2016](https://github.com/sabifo4/mammals/tree/master/01_SeqBayes_S1/01_BASEML/01_trees/00_rooted_trees_72sp/atlantogenata_scandentia_primates_tarver2016))
and fitted a skew-_t_ (ST) distribution to each of the estimated posterior
distributions for each node. The whole procedure is detailed in the [`01_SeqBayes_S1/03_Fit_ST_to_posteriors`](https://github.com/sabifo4/mammals_proj/tree/master/01_SeqBayes_S1/03_Fit_ST_to_posteriors)
directory. 

# 2. Second part of the sequential Bayesian dating approach 

## 2.1. Outline 
The mammal phylogeny that contains over 5,030 mammal taxa and 182 genes is too big to be analysed in a unique tree file
or a unique alignment. Therefore, we had to work with separate data subsets with less taxa to reduce computing time so
we could run the dating analysis in `MCMCtree`. As it has been previously mentioned in a previous tutorial file, taxa were separated according to the following taxonomical
groups: Afrotheria, Xenarthra, Marsupialia, Euarchonta, Lagomorpha, Laurasiatheria, Rodentia, or Monotremata. The number of taxa within Laurasiatheria and Rodentia,
however, was still too high to run a Bayesian dating analysis within a reasonable amount of time. Therefore, these two data sets were further divided into three subsets
(three subtrees for taxa within Rodentia and three subtrees for taxa within Laurasiatheria). The tags for the subtrees are the following:   

   * Afrotheria   
   * Marsupialia   
   * Xenarthra   
   * Euarchonta   
   * Lagomorpha   
   * Laurasiatheria_cetartiodactyla   
   * Laurasiatheria_chiroptera   
   * Laurasiatheria_therest   
   * Rodentia_squirrels   
   * Rodentia_ctenohystrica   
   * Rodentia_therest   

Note that all these 11 subtrees were rooted using Monotremata, which includes *Tachyglossus aculeatus*, *Zaglossus bruijni*, and *Ornithorhynchus anatinus*,
as the outgroup.

The filtered genes for these taxa are either mitochondrial or nuclear data. We generated the following 
alignments:   

   1) **`mt_3cp`**: mitochondrial data, only the concatenated third codon positions of the 
                mitochondrial genes found in the taxa included in the corresponding subtree.   
   2) **`mt_12cp`**: mitochondrial data, only the concatenated first and second codon positions of the 
                mitochondrial genes found in the taxa included in the corresponding subtree.   
   3) **`nt_3cp`**: nuclear data, only the concatenated third codon positions of the 
                nuclear genes found in the taxa included in the corresponding subtree.   
   4) **`nt_12cp`**: nuclear data, only the concatenated first and second codon positions of the 
                nuclear genes found in the taxa included in the corresponding subtree.   
   5) **`mt_rna`**: mitochondrial RNA data, concatenated genes found in the taxa included in
                the corresponding subtree.   
   6) **`all`**: all the previous 5 partitions concatenated, one after the other, for the taxa included in
                the corresponding subtree.   
   7) **`all_partitioned`**: all the first 5 partitions separated in individual blocks, one per partition.   

**NOTE: If you want to follow the details given below to reproduce our analyses, you may want to download**
**[this zip] (ADD LINK TO ZIP)  file. Specifically, it contains all the directories and subdirectories that we go through below,** 
**the results we obtained, and the scripts we used.**

If you have downloaded the zip file mentioned above, you might want to add the extracted `00_generate_trees`, 
`01_generate_alignments` directories to this specific location. If you have done this sucessfully, you should have the following directories: 

```
02_SeqBayes_S2 
    |- 00_generate_trees       <-- DIRECTORY EXTRACTED FROM THE ZIP FILE
    |- 01_generate_alignments  <-- DIRECTORY EXTRACTED FROM THE ZIP FILE
    |- 02_BASEML 
    |- 02_MCMCtree	
    |- README.md
```

## 2.2. Generating partitioned and concatenated alignments
Before proceeding to generate the final alignment and tree files, we can divide this step into two:   

   * A. Alignments already checked   
   * B. Alignments pending to check   
   
### A. Alignments already checked
Some of the alignments for the 11 subtrees to be evaluated had already been checked for possible taxonomic issues:   

   * Taxonomic mismatch   
   * Updated species name   
   * Removal of specific taxa (e.g., few genes, misslabeled species, etc.)   
   * Others   
   
For these subtrees (i.e., afrotheria, euarchonta, marsupialia, and xenarthra), an existing alignment 
with the required filtering steps to tackle the identified taxonomic issues had been created
(i.e., `alignment.phylip` file that can be found for each subtree directory inside 
`02_SeqBayes_S2/00_generate_trees/<name_subtree>/RAxML_n_filt`). Nevertheless, this alignment had all partitions concatenated
except for the `nt_3cp` one. This means that the final alignment file, `alignment.phylip`, 
was a concatenated alignment with 4 partitions: `mt_3cp`, `mt_12cp`, `nt_12cp`, and `mt_rna`.
Therefore, the first thing that we had to do was to apply the same filtering steps used to obtain this concatenated 
alignment to the unfiltered alignment with the `nt_3cp` partition (separate alignment phylip that can be 
found inside `02_SeqBayes_S2/00_generate_trees/<name_subtree>/RAxML_n_filt/checked_aln/original_aln_tree`).
The details for the steps carried out for this filtering step are detailed in the `README.md`
file under each of the subdirectories for the corresponding subtree. The architecture file follows 
the following scheme:  

```
02_SeqBayes_S2/
     |- 00_generate_trees/
             |- <name_subtree>/                      # i.e., "afrotheria", "euarchonta", etc.
                         |- RAxML_n_filt/ 
                              |- checked_aln/        # dir where processing of `3nt_cp` aln partition takes place
                              |     |- alignment_nt3cp.phylip        # filtered `3nt_cp` partition 
                              |     |- RAxML_bestTree.BS_ML_GTRCAT   # filtered tree    
                              |     |- *txt 
                              |     |- _other dirs and files_   # number of dirs vary from subtree to subtree, depends on filtering
                              |
                              |- alignment.phylip 
                              |- *.txt 
                              |- README.md      # file detailing step by step all the filtering process
                              |- summary.html   # html file that flags taxa under possible wrong placement in the tree
                              |- *.csv			 
                              |- *.RData
                              |- parse_lineage.R			 							  
```

>> _**There are other files under the `<name_subtree>` directory which are not added here as they are discussed subsequently for the next step**_

After the filtering steps have been applied to the `nt_3cp` partition (see `README.md` in each subdirectory for specific details 
followed to filter each data subset), the final filtered `nt_3cp` partition can be found inside `checked_aln` subdirectory named as 
`alignment_nt3cp.phylip`. The same applies to the filtered tree, `checked_aln/RAxML_bestTree.BS_ML_GTRCAT`.

> **NOTE**: Subtrees `afrotheria` and `euarchonta` needed to undergo a further filtering as there were some 
> subspecies that had not been removed. Therefore, the `README.md` inside the directories for these two  
> subtrees will contain an extra filtering step in which the concatenated filtered alignment (the one that  
> contains the filtered `nt_3cp` partition and the rest of the filtered alignment, details in the `README.md`) 
> was further processed to remove these taxa. For these subtrees, the directory `checked_aln` will contain 
> the concatenated alignment `alignment.phylip` too.

After that, the next step is to generate one alignment for each partition as well as an 
alignment with all the partitions concatenated. A detailed explanation about how to generate them 
can also be found in each `README.md` file for each subtree.
As a summary, just note that we ran the `02_SeqBayes_S2/01_generate_alignments/Concatenate_seqs_for_MCMCtree.R` R script 
from RStudio (updating general settings as discussed in each `README.md`). For more information about what to 
expect, take a look at the `README.md` file for each subtree. The final alignments can be found in a separate 
directory (`02_SeqBayes_S2/01_generate_trees`) which structure looks like this: 

```
02_SeqBayes_S2/ 
   |- 01_generate_alignments/
                  |- 00_mammals_alns 
                          |- <name_substree>/ 
                                  |- <name_subtree>.aln 
                                  |- <name_subtree>_mt3cp.aln 
                                  |- <name_subtree>_mt12cp.aln 
                                  |- <name_subtree>_mtrnacp.aln 
                                  |- <name_subtree>_nt3cp.aln 
                                  |- <name_subtree>_nt12cp.aln 
                                  |- partitions.txt 
								  
```

> **NOTE**: If you want to run on your own the scripts provided in the zip file to reproduce our analyses,
> you can move `02_SeqBayes_S2/01_generate_alignments/00_mammal_alns` and `02_SeqBayes_S2/01_generate_alignments/Rout` to another directory 
> and then run the the R script `Concatenate_seqs_for_MCMCtree.R` from RStudio as detailed in the `README.md` for 
> each subtree. Also, for `afrotheria` and `euarchonta` directories, you may find subdirectories labelled as 
> `<name_subtree>_old`, where alignments without the last filtering (remove subspecies) were saved.
> You will also see other files and directories that are not shown here because they are later discussed for subsequent steps.

### B. Alignments pending to check

The same filters applied to the alignments classified as group "A" were applied for this other subgroup "B". The difference is that we first 
concatenate the 1st+2nd CPs that had been already generated to the 3rd CPs. Note that no filtering was applied, which means 
that both alignment files are in the same state (i.e., identified as "cleaned"). All the details can be found in each 
`README.md` file inside the corresponding directories for each subtree.
In a nutshell, we did the following:   

   1. Concatenate `cleaned` alignment with 1st+2nd CP and `cleaned` alignment with 3rd CP   
   2. Apply filters: remove and rename taxa in both alignment and tree files.   
   3. Generate final filtered and concatenated alignments and pruned tree files.   

The directories for these subtrees look as it follows:

```
02_SeqBayes_S2/ 
    |- 00_generate_trees/
             |- <name_subtree_GROUP.A>/ # i.e., "afrotheria", "euarchonta", "marsupialia", "xenarthra"
             |- newly_checked/	 
                   |- <name_subtree_GROUP.B>/ # i.e., "lagomorpha", "laurasiatheria_X", "rodentia_X"
                         |- RAxML_n_filt/ 
                              |- checked_aln/        # dir where alignment is filtered
                              |     |- alignment.phylip              # filtered alignment with 1st+2nd CP and 3rd CP concatenated
                              |     |- RAxML_bestTree.BS_ML_GTRCAT   # filtered tree   
                              |     |- *txt 
                              |     |- _other dirs_  # number of dirs vary from subtree to subtree, depends on filtering
                              |    
                              |- checks_SAC/       # dir with notes about taxonomic filtering - might or not be present			  
                              |- alignment.phylip  # unfiltered alignment with 1st+2nd CP
                              |- *.txt 
                              |- RAxML_bestTree.BS_ML_GTRCAT # tree generated with unfiltered 1st+2nd CP alignment found at this same level
                              |- README.md         # file detailing step by step all the filtering process
                              |- summary.html
                              |- *.csv	| *.xlsx   # Excel or CSV file containing info about taxonomic filtering		 
                              |- *.RData
                              |- parse_lineage.R		

```
**NOTE**: Rodentia subtrees (labelled as `ctenohystrica`, `squirrel`, `rod_therest`) have an extra directory level:

```
rodentia
 |- RAxML_n_filt/ 
         |- <name_subtree> # i.e., "ctenohystrica", "rod_therest", "squirrel"
         |     |- SAC_check/     # dir with notes about taxonomic filtering - might be deleted eventually       
         |     |- checked_aln/   # dir where alignment is filtered
         |     |     |- alignment.phylip              # filtered alignment with 1st+2nd CP and 3rd CP concatenated
         |     |     |- RAxML_bestTree.BS_ML_GTRCAT   # filtered tree   
         |     |     |- *txt 
         |     |     |- _other dirs_        # number of dirs vary from subtree to subtree, depends on filtering
         |     |
         |     |- *.txt 
         |     |- RAxML_bestTree.BS_ML_GTRCAT # tree generated with the unfiltered 1st+2nd CP alignment found at this same level
         |     |- alignment.phylip            # unfiltered alignment with 1st+2nd CP
         |     |- *.csv	| *.xlsx              # XLSX or CSV file containing info about taxonomic filtering			 
         |     |- summary.html					  
         | 
         |- *.txt 
         |- README.md            # file detailing step by step all the filtering process
         |- *.RData
         |- parse_lineage.R		
```
>> _**There are other files under the `<name_subtree>` directory which are not added here as they are discussed subsequently for the next step**_


## 2.3. Generate newick files for each subtree   

For each of the subtrees, we had the best-scoring maximum-likelihood (ML) tree estimated by [`RAxML`](https://github.com/stamatak/standard-RAxML).
This ML tree was then pruned according to the filtering steps mentioned above. Then, the pruned ML tree was used to manually
add the name of the node calibrations, which were going to be later replaced by the corresponding calibrations using an R script. 

The directory `02_SeqBayes_S2/00_generate_trees` contains one directory per subtree  
following the format detailed above. Inside each of them, you will find the resulting filtered tree file as well as the input and output files the R script uses to add the calibrations 
to the corresponding calibrated nodes. In general, although you might observe some differences depending on additional filtering steps 
these data subsets might have needed, these directories have the following architecture:    

```
02_SeqBayes_S2/ 
    |- 00_generate_trees/
            |- <name_subtree>/  # i.e., "afrotheria", "euarchonta", etc.
                |- RAxML_n_filt/ # dir where alignment and tree files were processed
                |
                |- <name_subtree>.tree # final tree output after filtering 				
                |
                |- <name_subtree>_rooted_baseml.tree # tree to be used by BASEML to compute gradient and Hessian
                |
                |- Xsp_<name_subtree>_MCMCtree_calib.tree  # X corresponds to the number of taxa this 
                |                                          # subtree has. File output by the R script
                |	
                |- <name_subtree>_rooted_calibnames.tree  # Input file used by the R script 
                |                                         # It has the names of calibrations manually added
                |
                |- Xsp_<name_subtree>_spnameslist.tree    # X corresponds to the number of taxa this 
                |                                         # subtree has. Summary output file by the R script 
                |    
                |- Calibrations_<name_subtree>.R   # R script that parses the data   
                |    	
                |- Calibrations_<name_subtree>.txt # Input file used by the R script with the matching nmes 
                                                   # of the calibrations and the actual calibrations in MCMCtree format
```

Each of the R scripts can be run within the corresponding subdirectories. Note that they also generate a dummy alignment 
that can be used in the subsequent steps as an input for `MCMCtree`. This is because the actual alignment is not used when 
`MCMCtree` uses the approximate likelihood to infer divergence times -- it uses the `in.BV` file generated by `BASEML` with  
the gradient and the Hessian (see [this tutorial](https://github.com/sabifo4/mammals/tree/master/01_SeqBayes_S1/01_BASEML/02_Hessian) as a reminder).
Therefore, in order to save space, we can use a dummy alignment that contains the header in 
PHYLIP format and the taxa names. Specifically, each dummy alignment has been saved in `02_SeqBayes_S2/01_generate_alignments/01_mammal_dummy_alns`.

>> Some of the subtrees needed some of the original calibrations to be modified when running `MCMCtree` without the data 
>> to evaluate if the user-specified prior matches the one the dating software uses. For these subtrees
>> (i.e., laurasiatheria_therest, rodentia_therest, euarchonta, laurasiatheria_cetartiodactyla), the directory will have an additional 
>> R script, an additional file with calibrations, and an additional tree file. These files names will end with `STtweak` before 
>> the extension. Note that marsupialia had some of the soft bound calibrations modified, but this was manually changed in the 
>> calibrations file, thus this directory will not have extra files as the change was manually included.

## 2.4. Run `BASEML` to estimate Hessian and gradient
At this stage, we can run `BASEML` (you may check the tutorial used in the first step [here](https://github.com/sabifo4/mammals/tree/master/01_SeqBayes_S1/01_BASEML/02_Hessian))
with each of the subtrees to estimate the Hessian and the gradient. 

You can go to the [`02_BASEML`](https://github.com/sabifo4/mammals/tree/master/02_SeqBayes_S2/02_BASEML) directory,
in which a description of the subsequent analyses with the files output by `BASEML` is given. Note that the files concerning 
this analysis are only available in the zip file that you might have downloaded at the beginning of this tutorial.

# 2.5. Bayesian inference of divergence times
Once the Hessian and gradient have been estimated for each data subset, we can run `MCMCtree` with the approximate likelihood 
calculation ([dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613)). 

You can find all the details about the whole procedure, including the safety checks prior to proceed with the Bayesian dating analysis, 
in the [`02_MCMCtree` directory](https://github.com/sabifo4/mammals/tree/master/02_SeqBayes_S2/02_MCMCtree). Note that the files concerning 
these analyses are only available in the zip file that you might have downloaded at the beginning of this tutorial.
