# Directory with calibrated trees 
This directory is structured as it follows:   

   * One directory for each tree topology ($7$ directories)   
   * `README.md` (this file)   

Within each directory, you will find the following:   

   * Tree file used to include calibrations with the R script [`Include_calibrations.R`](01_SeqBayes_S1/01_BASEML/01_trees/01_calibtrees_72sp/Include_calibrations.R),
   e.g., `72sp_atlantogenata_tarver2016_calibrated.tree`.   
   * Calibrated tree file in newick format, e.g., `72sp_MCMCtree_atlantogenata_tarver2016_calibrated.tree`.   
   * `FigTree` file, e.g., `FIGTREE_72sp_atlantogenata_tarver2016_calibrated`; in which the calibrated tree has already been loaded with the calibration labels
     and is presented in a cladogram view. Open with [`FigTree`](http://tree.bio.ed.ac.uk/software/figtree/). 
   * Tree with the calibrated labels in pdf format. This is the output you get if
     you load the `FigTree` file in [`FigTree`](http://tree.bio.ed.ac.uk/software/figtree/).   

The main tree topology is T2 ([here](01_SeqBayes_S1/01_BASEML/01_trees/01_calibtrees_72sp/02_T2/72sp_atlantogenata_tarver2016_MCMCtree_calib.tree)). 

# Notes 
Some of the calibrations could not be included in the following
tree hypotheses because of different clustering:   

   * Tree hypothesis 5:   
      * Misses `EUUNGULATA` because _Equus caballus_ is not clustered
	    with Artiodactyla but as an outgroup of carnivores.   
      * Misses `ROD-NOSQUIRREL` because *Ictidomys tridecemlineatus* is
        placed as an outgroup of the group which calibration is `DIPOD-RATT`.   

   * Tree hypotheses 3 and 4:   
      * Miss `EUUNGULATA` because _Equus caballus_ is not clustered
	    with Artiodactyla but as an outgroup of carnivores.   
   
	   