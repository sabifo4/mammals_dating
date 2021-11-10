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
python3.8 ../scripts/72sp_subtree_graft.py
```

To clean the directory and remove tree files that you will not need,
you should run the following commands:

```sh
# Run from `input_trees` 
rm Working*tree All_Laurasiatheria.tree
```

Now, you can run the next python script so you can generate 
the *nexus and *nwk tree files that you may use in macroevolutionary analysis
(e.g., trait regression):

```sh 
# Run from `input_trees` 
python3.8 ../scripts/4705sp_node_ages.py
```

>>**NOTES**: We used Python v3.8.12 and BioPython v1.79 on the WLS 
>>and Python3.8.5 and BioPython v1.77 on Mac OSX. 
>>Also, if you have the WLS and use Python v3.8.12 and BioPython v1.77
>>there seems to be an issue with the 95% CIs as they are not printed out.
>>If you upgrade BioPython to v1.79, this is sorted out when using 
>>Python v3.8.12 on the WLS.
>>Last, you might see that some `from` tags in the output nexus files do not 
>>have the same names as we renamed the file names for the subtrees and 
>>the backbone tree so it matches the data in `FigShare`.

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

