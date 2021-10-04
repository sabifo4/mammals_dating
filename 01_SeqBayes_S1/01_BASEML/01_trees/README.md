# 1. Generate rooted trees for 7 tree hypothesis
We used the best-scoring maximum-likelihood tree with 72 species 
that was obtained with `RAxML` v8.2.12 ([Stamatakis 2014](https://github.com/stamatak/standard-RAxML)),
see [here](ML_tree_72sp/RAxML_bestTree.concatenated.rooted.tree),
as a guide tree to generate the 7 tree hypotheses that have been most debated with regards to the placement of  
taxa in the mammal tree of life. Note that we used *Ornithorhynchus anatinus* as an outgroup 
to root the trees.

To easily visualize the new taxa that we added to each clade in comparison to the study carried out by 
[dos Reis et al. 2012](https://royalsocietypublishing.org/doi/full/10.1098/rspb.2012.0683?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%3dpubmed)
with 43 mammal taxa, we generated the nexus 
file [`Coloured_diffs_43spVS72sp_nexus`](ML_tree_72sp/Coloured_diffs_43spVS72sp_nexus.tree).
If loaded in [`FigTree`](http://tree.bio.ed.ac.uk/software/figtree/), 
you can see in red the names of the new taxa. In addition, we coloured the branch lengths leading to their 
sister taxa so it is easier to find the common node in the 7 tree hypotheses.

The final trees with 72 species can be found in the following 
directories according to each tree hypothesis. The trees with 43 taxa for each hypothesis 
have been also included in the corresponding directories as they were used as backbones to add the new
taxa:   

   * [Tree 1 - T1](00_rooted_trees_72sp/01_T1): tree hypothesis 1,
   Atlantogenata rooting and *Tupaia belnageri* sister clade to Primates.   
   * [Tree 2 - T2](00_rooted_trees_72sp/02_T2): tree hypothesis 2,
   Atlantogenata rooting and *Tupaia belnageri* sister clade to Glires.
   This is the tree hypothesis used in the subsequent downstream analyses and later referred to as "main tree".   
   * [Tree 3 - T3](00_rooted_trees_72sp/03_T3): tree hypothesis 3,
   same as T1 but Chiroptera and *Equus caballus* are placed as sister clades to Carnivora.   
   * [Tree 4 - T4](00_rooted_trees_72sp/04_T4): tree hypothesis 4,
   same as T3 but *Vicugna pacos* and *Sus scrofa* exchange their placement.   
   * [Tree 5 - T5](00_rooted_trees_72sp/05_T5): tree hypothesis 5,
   topology arranged as in the [`Ensembl` species tree](http://www.ensembl.org/info/about/speciestree.html), that is, same as T4 but *Tupaia belangeri* is sister clade to
   Primates.   
   * [Tree 6 - T6](00_rooted_trees_72sp/06_T6): tree hypothesis 6,
   Epitheria rooting and same placement for the other taxa as T2.   
   * [Tree 7 - T7](00_rooted_trees_72sp/07_T7): tree hypothesis 7,
   Exafroplacentalia rooting and same placement for the other taxa as T2.   

Note that only one tree hypothesis is going to be used during the second step of the 
sequential Bayesian dating: **T2**. To easily track which taxa are differently placed in each of the 
seven tree hypotheses, you can check the PDF and the SVG files generated for each of tree hypothesis 
[here](02_tree_hypotheses_figures/). 
These taxa/clades have been identified with different colours: (1, blue) Xenarthra, (2, orange) Afrotheria,
(3, green) Chiroptera, (4, blue) *Equus caballus*, (5, yellow) *Vicugna pacos*, (6, purple) *Sus scrofa*,
(7, red) *Tupaia belangeri*. 

>>**NOTE**: for a thorough discussion on the different mammal tree hypotheses,
>> you may read [Tarver et al., 2016 ](https://pubmed.ncbi.nlm.nih.gov/26733575/) 
>> and [dos Reis et al. 2012](https://royalsocietypublishing.org/doi/full/10.1098/rspb.2012.0683?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%3dpubmed).

# 2. Calibrate rooted trees for 7 tree hypotheses
Once the 7 rooted trees with 72 mammal species are generated, we need to calibrate the nodes for which
we have fossil information. We calibrated the root of the tree and 31 internal nodes; 32 calibrations 
in total. These are included as "soft bounds", that is, uniform distributions defined by a 
maximum and a minimum age (based on the fossil record) between which the true age is expected to lie.
Besides, the `MCMCtree` user defines left and right tail probabilities to the maximum and minimum bounds
so they are not hard bounds. 
In this study, we set the tail probabilities to $p=0.25$, which means that there is a probability of 2.5%
that the true age is beyond the bound. 

| Calibration name      | Calibration (`MCMCtree`) |
|-----------------------|------------------------|
| MAMMALIA              | 'B(1.649,2.51254)'     |
| THERIA                | 'B(1.2156,1.696)'      |
| PLACENTALIA           | 'B(0.616,1.646)'       |
| EUARCHONTOGLIRES      | 'B(0.616,1.646)'       |
| PRIMATES              | 'B(0.56,0.6611)'       |
| ANTHROPOIDEA          | 'B(0.339,0.6611)'      |
| CATARRHINI            | 'B(0.2444,0.339)'      |
| HOMINIDAE             | 'B(0.1163,0.339)'      |
| HOMININAE             | 'B(0.0533,0.339)'      |
| HOMININI              | 'B(0.065,0.1)'         |
| CERCOPITHECINAE       | 'B(0.0533,0.34)'       |
| PAPIONINI             | 'B(0.053,0.339)'       |
| STREPSIRRHINI         | 'B(0.339,0.6611)'      |
| GLIRES                | 'B(0.56,1.646)'        |
| RODENTIA              | 'B(0.56,0.6611)'       |
| NONSQUIRREL RODENTS   | 'B(0.476,0.592)'       |
| DIPODIDAE-MUROIDEA    | 'B(0.407,0.592)'       |
| MURINAE               | 'B(0.072,0.16)'        |
| LAGOMORPHA            | 'B(0.476,0.6611)'      |
| EUUNGULATA            | 'B(0.524,0.6611)'      |
| ARTIODACTYLA          | 'B(0.505,0.6611)'      |
| CETRUMINANTIA         | 'B(0.505,0.6611)'      |
| BOVIDAE               | 'B(0.16,0.281)'        |
| CARNIVORA             | 'B(0.373,0.6611)'      |
| CANIFORMIA            | 'B(0.373,0.6611)'      |
| CHIROPTERA            | 'B(0.4760,0.6611)'     |
| LIPOTYPHLA            | 'B(0.6160,1.6460)'     |
| XENARTHRA             | 'B(0.476,1.646)'       |
| AFROTHERIA            | 'B(0.56,1.646)'        |
| PAENUNGULATA          | 'B(0.56,1.646)'        |
| MARSUPIALIA           | 'B(0.4760,1.313)'      |
| EOMETATHERIA          | 'B(0.2303,0.56)'       |

These calibrations can also be found in the [`Calibrations_mammalia.txt`](01_calibtrees_72sp/Calibrations_mammalia.txt) text file, which 
has two columns delimited by a `|`. The column on the right has the names of the nodes that are to be calibrated 
and the one on the left the corresponding soft calibrations in `MCMCtree` format (time unit = 100 Mya).
>> **NOTE**: You will see that there are four nodes in the text file that have different tag names if compared to the ones included
>> in the table above: "HOMININAE" is "GORILLA-HUMAN", "STREPSIRRHINI" is "STREPSIRHINI (we had a typo found at the end of the study. 
>> We have not changed it in the file as it might affect downstream analyses), "NONSQUIRREL RODENTS" is "ROD-NOSQUIRREL",
>> "DIPODIDAE-MUROIDEA" is "DIPOD-RATT", "MURINAE" is "MURIDAE" (another typo we found later), "CETRUMINANTIA" is "WHIP-RUM", and "BOVIDAE" is "BOV-ANTIL".

First, we manually included the names of the nodes to be calibrated in the newick trees for each of 
the 7 hypotheses (files that end with `*_calibrated.tree` and are saved in each of the seven directories 
inside the [`01_calibtrees_72sp`](01_calibtrees_72sp)
directory). Then, we used the [`Include_calibrations.R`](01_calibtrees_72sp/Include_calibrations.R)
R script to replace them with the soft calibrations in `MCMCtree` format detailed in the 
[`Calibrations_mammalia.txt`](01_calibtrees_72sp/Calibrations_mammalia.txt)
text file.

To ease the visualization of the final calibrated trees and manually check the calibrations were placed 
in the correct nodes, we generated two files:    

   * `FIGTREE_72sp_<name_tree_hypothesis>_calibrated.pdf`: PDF file with the final tree hypothesis (cladogram format) with the names 
   of the calibrations added to the corresponding nodes. Each directory within the [`01_calibtrees_72sp`](01_calibtrees_72sp)
   contains this specific file.   
   * `FIGTREE_72sp_<name_tree_hypothesis>_calibrated`: File to be opened with [`FigTree`](http://tree.bio.ed.ac.uk/software/figtree/).
   The user can then zoom in and out and change the display of the tree.   
   
