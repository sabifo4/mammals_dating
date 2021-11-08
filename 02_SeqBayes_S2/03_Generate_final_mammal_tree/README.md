# The Mammal Tree of Life 
Here, you can find the final mammal Tree of Life after we grafted the different 
data subsets onto the backbone tree (i.e., combining the results from the first step 
of the sequential Bayesian dating approach with the those obtained during the 
second step). The scripts that we have used for that purpose can be found
[here](scripts). The tasks they perform are the following:   

   * [4705sp_node_ages.py](scripts/4705sp_node_ages.py): this script builds the 4,705-sp tree by grafting the subtrees onto
   the backbone tree.   
   * [72sp_subtree_graft.py](scripts/72sp_subtree_graft.py): this script creates the additional teres with node ages from CI.   

An explanation for each output file can be found below:   

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

