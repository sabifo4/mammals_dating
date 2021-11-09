# The Mammal Tree of Life 
Here, you can find the final mammal Tree of Life after we grafted the different 
data subsets onto the backbone tree (i.e., combining the results from the first step 
of the sequential Bayesian dating approach with the those obtained during the 
second step). The scripts that we have used for that purpose can be found
[here](scripts) and the input trees [here](input_trees). The tasks each script performs are the following:   

   * [4705sp_node_ages.py](scripts/4705sp_node_ages.py): this script builds the 4,705-sp tree by grafting the subtrees onto
   the backbone tree.   
   * [72sp_subtree_graft.py](scripts/72sp_subtree_graft.py): this script creates the additional teres with node ages from CI.   

You can generate the output files we detail above with the mammal timetree 
by running the following commands :

```sh 
# Run from the `input_trees` directory
# to generate the output files you see
# in this directory 
python3 ../scripts/72sp_subtree_graft.py
```

When you run this, you will see the following messages
printed out on the terminal screen:

```
---------- Subtree Chiroptera_subt1-time, 253 taxa ----------
- Applying rule 1
MRCA age: 0 -> 0.560476
Adjusted branch length: 0.594646 -> 0.034170000000000034
Replacing 1 taxa with 252


---------- Subtree Chiroptera_subt2-time, 631 taxa ----------
- Applying rule 1
MRCA age: 0 -> 0.5415110000000001
Adjusted branch length: 0.594646 -> 0.05313499999999993
Replacing 1 taxa with 630


---------- Subtree Artiodactyla-time, 428 taxa ----------
- Applying rule 1
MRCA age: 0.570244 -> 0.5790270000000001
Adjusted branch length: 0.082329 -> 0.0735459999999999
Replacing 2 taxa with 428


---------- Subtree Afrotheria-time, 57 taxa ----------
- Applying rule 1
MRCA age: 0.658755 -> 0.7076720000000001
Adjusted branch length: 0.127826 -> 0.0789089999999999
Replacing 3 taxa with 57


---------- Subtree Euarchonta-time, 483 taxa ----------
- Applying rule 1
Backbone tree has {'carlito_syrichta', 'saimiri_boliviensis_boliviensis'} but at not in subtree
MRCA age: 0.6489510000000001 -> 0.6653889999999999
Adjusted branch length: 0.048929 -> 0.03249100000000016
Replacing 24 taxa with 461

- Applying rule 2
MRCA age: 0 -> 0.5806549999999999
Adjusted branch length: 0.689588 -> 0.10893300000000006
Replacing 1 taxa with 22


---------- Subtree Lagomorpha-time, 85 taxa ----------
- Applying rule 1
MRCA age: 0.476208 -> 0.48311099999999996
Adjusted branch length: 0.18579 -> 0.17888700000000007
Replacing 2 taxa with 85


---------- Subtree Marsupialia-time, 304 taxa ----------
- Applying rule 1
Backbone tree has {'notamacropus_eugenii'} but at not in subtree
MRCA age: 0.5815349999999999 -> 0.5367989999999999
Adjusted branch length: 0.840268 -> 0.885004
Replacing 3 taxa with 304

- Applying rule 2
MRCA age: 0 -> 0.66537
Adjusted branch length: 1.997583 -> 1.332213
Replacing 1 taxa with 3


---------- Subtree Xenarthra-time, 30 taxa ----------
- Applying rule 1
MRCA age: 0.622112 -> 0.635598
Adjusted branch length: 0.164468 -> 0.150982
Replacing 2 taxa with 30


---------- Subtree Sciuridae_and_related-time, 264 taxa ----------
- Applying rule 1
MRCA age: 0 -> 0.586579
Adjusted branch length: 0.60782 -> 0.021241000000000065
Replacing 1 taxa with 264


---------- Subtree Ctenohystrica-time, 207 taxa ----------
- Applying rule 1
Backbone tree has {'heterocephalus_glaber_female'} but at not in subtree
MRCA age: 0.41351499999999997 -> 0.538855
Adjusted branch length: 0.179647 -> 0.054306999999999994
Replacing 6 taxa with 206


---------- Subtree Rodentia_therest_subt1-time, 627 taxa ----------
- Applying rule 1
Backbone tree has {'mus_spretus', 'mus_pahari', 'mus_caroli', 'cricetulus_griseus_chok1gshd'} but at not in subtree
MRCA age: 0.548631 -> 0.5366660000000001
Adjusted branch length: 0.044531 -> 0.056495999999999894
Replacing 12 taxa with 627


---------- Subtree Rodentia_therest_subt2-time, 688 taxa ----------
- Applying rule 1
MRCA age: 0.168637 -> 0.21713299999999994
Adjusted branch length: 0.102742 -> 0.05424600000000007
Replacing 2 taxa with 686


---------- Subtree All_Laurasiatheria, 1962 taxa ----------
- Applying rule 1
MRCA age: 0.697734 -> 0.7378299999999999
Adjusted branch length: 0.053029 -> 0.01293300000000009
Replacing 17 taxa with 1962
```

To clean the directory, you should run the following commands:

```sh
# Run from `input_trees` 
rm Working*tree All_Laurasiatheria.tree
```

Now, you will see the `4705sp.tree` file inside this directory, which 
is the mammal timetree that you also have in this main directory. Now,
to generate the rest of the output files, you can run the following 
commands: 

```sh 
# Run from `input_trees` 
python3 ../scripts/4705sp_node_ages.py
```

Once this script finishes, you will see the *nexus and *nwk files 
that you have in this directory.
An explanation for each of the output files obtained in each step
is given below:   

  * [4705sp.tree](4705sp.tree): this is the final mammal timetree in NEXUS format. It includes the confidence intervals (CIs) for the
  estimated mean divergence times.   
  * [4705sp_colours.tree](4705sp_colours.tree): this is the final mammal timetree in NEXUS format. It has the instructions to plot the 
  stitched tree with `FigTree` in radial format, with 95% CIs of node ages, and with the main clades coloured; all as in the
  main figure in the paper.   
  * [4705sp_mean.nwk](4705sp_mean.nwk): this is the final mammal timetree in newick format. The CIs are not included here. You can
  use these trees in macroevolutionary analysis (e.g., trait regression).   
  * [4705sp_ci025.nexus](4705sp_ci025.nexus): this is a mammal timetree in NEXUS format which mean divergence times correspond to the minimum
  ages that were estimated in the 2.5% CIs included in the final mammal timetree. The CIs
  estimated with the final mammal timetree have also been included.   
  * [4705sp_ci025.nwk](4705sp_ci025.nwk): this is a mammal timetree in newick format which mean divergence times correspond to the minimum
  ages that were estimated in the 2.5% CIs included in the final mammal timetree. CIs are not included. You can
  use these trees in macroevolutionary analysis (e.g., trait regression).   
  * [4705sp_ci975.nexus](4705sp_ci975.nexus): this is a mammal timetree in NEXUS format which mean divergence times correspond to the maximum
  ages that were estimated in the 97.5% CIs included in the final mammal timetree. The CIs
  estimated with the final mammal timetree have also been included.   
  * [4705sp_ci975.nwk](4705sp_ci975.nwk): this is a mammal timetree in newick format which mean divergence times correspond to the minimum
  ages that were estimated in the 97.5% CIs included in the final mammal timetree. The CIs
  estimated with the final mammal timetree have also been included. You can
  use these trees in macroevolutionary analysis (e.g., trait regression).   
  * [4705sp_names.txt](4705sp_names.txt): this file contains a list of the 4,705 mammal species included in the final mammal timetree.   

