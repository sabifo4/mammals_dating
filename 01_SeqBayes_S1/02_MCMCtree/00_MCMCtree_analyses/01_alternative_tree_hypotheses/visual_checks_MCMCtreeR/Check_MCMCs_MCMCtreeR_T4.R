#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#--------------#
# LOAD PACKAGE #
#--------------#
library( MCMCtreeR )

#-----------------------#
# SET WORKING DIRECTORY #
#-----------------------#
library( rstudioapi ) 
# Get the path to current open R script and find main dir
path_to_file <- getActiveDocumentContext()$path
wd <- paste( dirname( path_to_file ), "/", sep = "" )
setwd( wd )
wd_trees <- gsub( pattern = "01_visual_checks_MCMCtreeR/", replacement = "", x = wd )
source( "Functions_plots_MCMCtreeR.R" )

#-------------#
# START TASKS #
#-------------#

# 1. Load trees and mcmc files for master tree: 04_atlantogenata_tarver2016_laurasiatheriaensembl 
## GBM - posterior
phy.post.GBM <- mcmc.post.GBM <- vector( mode = "list", 10 )
for( i in seq(1:length( phy.post.GBM ) ) ){
        
  cat( "Parsing files for run ", i, "...\n" )
  phy.post.GBM[[ i ]] <-  MCMCtreeR::readMCMCtree( paste( wd_trees, "00_MCMCtree/mcmc",
                                                            i, "/mcmctree_GBM/04_atlantogenata_tarver2016_laurasiatheriaensembl/FigTree.tre", sep = "" ),
                                                     from.file = TRUE )
  mcmc.post.GBM[[ i ]] <- read.table( paste(  wd_trees, "00_MCMCtree/mcmc",
                                                i, "/mcmctree_GBM/04_atlantogenata_tarver2016_laurasiatheriaensembl/mcmc.txt", sep = "" ),
                                        sep = "\t", header = T, stringsAsFactors = F )
}

# 2. Create objects in APE format 
## GBM - posterior
phy.post.GBM.obj <- phy.post.GBM.edge <- mcmc.chain.GBM <- vector( mode = "list", 10 )
for( i in seq(1, length( phy.post.GBM.obj ) ) ){
        cat( "Parsing files from run ", i, "...\n" )
        phy.post.GBM.obj[[ i ]]  <- phy.post.GBM[[ i ]]$apePhy
        phy.post.GBM.edge[[ i ]] <- phy.post.GBM.obj[[ i ]]$edge
        mcmc.chain.GBM[[ i ]]    <- mcmc.post.GBM[[ i ]]
        
}
# 3. Create PHY object with GBM and ILN objects 
phy.all <- vector( mode = "list", length = 10 )
names( phy.all ) <- c( paste( "GBM.post.", 1:10, sep = "" ) )
for( i in seq(1:length( phy.all ) ) ){
        
  print( i )
  cat( "Parsing GBM.post files from run ", i, "...\n" )
  phy.all[[ i ]] <- phy.post.GBM.obj[[ i ]]
  
}


# 4. Extract ages with node age posteriors from column 2
## GBM - posterior
mcmc.node.ages.post.GBM <- all.nodes.post.GBM <- vector( mode = "list", 10 )
for( i in 1:length( mcmc.node.ages.post.GBM ) ){
        mcmc.node.ages.post.GBM[[ i ]] <- mcmc.chain.GBM[[ i ]][, 2:Ntip( phy.all[[ i ]] )] 
        all.nodes.post.GBM[[ i ]]      <- as.numeric( gsub( "t_n", "",
                                                            colnames( mcmc.node.ages.post.GBM[[ i ]] ) ) )
}

# 5. Create a vector of names for each list element as  
#    internal nodes from APE tree, using phy$edge object 
## GBM - posterior
node.ages.names.post.GBM <- vector( mode = "list", 10 )
for( i in 1:length( node.ages.names.post.GBM ) ){
        node.ages.names.post.GBM[[ i ]] <- c( Ntip( phy.all[[ i ]] ) + 1,
                                              phy.post.GBM.edge[[ i ]][ which( phy.post.GBM.edge[[ i ]][,2] > Ntip( phy.all[[ i ]] ) ), 2 ] )
        
        
}


# 6. Find where each posterior node age appears in APE 
#    edge object.
## GBM - posterior
match.nodes.post.GBM <- vector( mode = "list", 10 )
for( i in 1:length( match.nodes.post.GBM ) ){
        match.nodes.post.GBM[[ i ]] <- match( all.nodes.post.GBM[[ i ]],
                                              as.numeric( node.ages.names.post.GBM[[ i ]] ) )
        
}

# 7. Create a list extracting the info from the mcmc
#    chain in APE node order 
## GBM - posterior
node.ages.post.GBM <- vector( mode = "list", 10 )
for( i in 1:length( node.ages.post.GBM ) ){
        node.ages.post.GBM[[ i ]] <- lapply( match.nodes.post.GBM[[ i ]],
                                             function( uu ) mcmc.node.ages.post.GBM[[ i ]][, uu ] )
}

# 8. Name each element in list 
## GBM - posterior
for( i in 1:length( node.ages.names.post.GBM ) ){
        names( node.ages.post.GBM[[ i ]] ) <- node.ages.names.post.GBM[[ i ]]
        
}

# 9. Create node.ages.all object to get all GBM and ILN objects
node.ages.all <- vector( mode = "list", length = 10 )
for( i in 1:length( node.ages.all ) ){
        
  node.ages.all[[ i ]] <- node.ages.post.GBM[[ i ]]
  
}

# 10. Get means for edges
## GBM - posterior
# Initialize matrix to get the edges of the 10 runs 
phy.matrix.post.GBM.edges <- matrix( 0, nrow = 10, ncol = length( phy.all[[1]]$edge.length ) )
for( i in 1:10 ){
        phy.matrix.post.GBM.edges[i,] <- phy.all[[ i ]]$edge.length
}
# Get mean of each edge. Also generate one mean without chains 4, 7, 8 -- according to visualization
phy.mean.edge.post.GBM      <- apply( X = phy.matrix.post.GBM.edges, 2, mean )
phy.mean.edge.post.GBM.filt <- apply( X = phy.matrix.post.GBM.edges[c(1:3,5:6,9:10),], 2, mean ) 
# Create copy of ape object and then replace the edge.length values accordingly
phy.mean.post.GBM                   <- phy.mean.post.GBM.filt <- phy.all[[1]]
phy.mean.post.GBM$edge.length       <- phy.mean.edge.post.GBM
phy.mean.post.GBM.filt$edge.length  <- phy.mean.edge.post.GBM.filt

# 11. Get means for ages
#------------------------------------#
# FUNCTION TO GET THE MEAN DISTANCES # 
#------------------------------------#
# Arguments:
#    len              numeric, number of internal nodes for which samples have been
#                     collected. Default: 71.
#    chains           numeric, indexes of the chains through which iterate to get the 
#                     collected samples for each age. Default: 1:10.
#    node.ages.list   list, list previously created with "chain" amount of entries 
#                     where the node ages each chain have been saved. Default: 10 chains,
#                     although it can be less, e.g., "prior.GBM".
#                     These objects "node.ages.list" could be "node.ages.post.GBM",
#                     or "node.ages.post.ILN", for instance.
#    tmp              numeric, temporary empty vector of length amount of samples 
#                     collected during the MCMCMs. Default: 20,001.
#    out              list, empty list previously created with "len" amount.
#                     E.g., "node.mean.dist.ages.post.GBM"
get_mean_dist_ages <- function( len = 71, chains = 1:10, node.ages.list, tmp ){
        
        # Create a list with 71 node entries, which we will output. We will then iterate
        # over each of the entries for each node sampled during each run. We will then
        # save the corresponding sampled vaules for this node 
        # to a the corresponding entry in the new created list.
        out          <- vector( mode = "list", length( node.ages.list[[ 1 ]] ) )
        names( out ) <- names( node.ages.list[[ 1 ]] ) 
        
        # Iterate over each node. We assume there are len = 71.
        # This will be modified if needed through the argument "len"
        for( j in 1:len ){
                
                # Reset local vars to access tmp vector
                a  <- 1
                b  <- 20001
                cte <- 20000
                
                # Iterate through every run and ten collect values for node_j
                for( i in chains ){
                        tmp[a:b] <- unlist( node.ages.list[[ i ]][j] ) # These are 20,001 values for node_j for run_i  
                        a <- b+1
                        b <- a+cte
                }
                
                # Save all the 20,001*10 values that have been sampled across the 10 runs 
                # and saved in vector "tmp"
                # Note that "out" is the list object with 71 entries, e.g., "node.mean.dist.ages.post.GBM"
                out[[ j ]] <- tmp 
                
        }
        
        # Return out
        return( out )
}

# a) Create tmp vector to allocate memory with 10*20,001 values
tmp <- rep( 0, 10*length( unlist( node.ages.post.GBM[[ 1 ]][1] ) ) )
# b) Loop from node 1 to 71 and then go through the corresponding values sampled for 
#    this node for each run. Use tmp vector to save values
node.mean.dist.ages.post.GBM <- get_mean_dist_ages( len = length( node.ages.post.GBM[[ 1 ]] ),
                                                    chains = 1:10,
                                                    node.ages.list = node.ages.post.GBM,
                                                    tmp = tmp )

# TRY PLOT [[ IT WORKS! It returns the last plot ]]
# mcmc.tree.plot.RETPLOT( phy = phy.all$GBM, node.ages = node.ages.GBM, 
#                  analysis.type = "user", cex.tips = 0.2,
#                  time.correction = 100, 
#                  scale.res = c( "Eon", "Period" ), plot.type = "distributions",
#                  cex.age = 0.4, cex.labels = 0.5, relative.height = 0.08, 
#                  col.tree = "grey40", no.margin = TRUE,
#                  density.col = "blue", add.time.scale = TRUE)

#-------#
# PLOTS #
#-------#

#=============================#
# 0. PLOT: Evalute GBM chains #
#=============================#
## THE NEXT LINES CAN BE RUN AT ONCE
## THE OUTPUT FILE IS SAVED AS 
##    "T4_Explore_GBM_runs.pdf"
# Set transparent colour
transp.col <- adjustcolor( col = "blue", alpha.f = 0.3 )
last.plot  <- mcmc.tree.plot.RETPLOT( phy = phy.all$GBM.post.1, xlim.scale = c(-20,300),
                                      node.ages = node.ages.post.GBM[[ 1 ]], 
                                      show.tip.label = TRUE,
                                      analysis.type = "user", cex.tips = 0.8,
                                      time.correction = 100, 
                                      scale.res = c( "Eon", "Period" ), plot.type = "distributions",
                                      cex.age = 0.6, cex.labels = 0.8, relative.height = 0.08, 
                                      col.tree = "grey40", no.margin = TRUE,
                                      density.col = transp.col, add.time.scale = TRUE, grey.bars = FALSE #,density.border.col = "red"
)

# Add now the rest of GBM distributions with this helper function
col.mod = "grey"
transparency = 0.3 
col.plot <- c( "white", "black", "darkgreen", "red", "darkblue", "orange", "pink",
              "lightblue", "cyan", "brown" )
for( i in 2:10 ){
        add.extra.dists( phy = phy.all[[ i ]], num.models = 1, last.plot = last.plot,
                         node.ages =  node.ages.post.GBM[[ i ]], plot.type = "distributions",
                         time.correction = 100,
                         density.col = "white", density.border.col = col.plot[ i ],
                         distribution.height = 0.8, transparency = transparency
        )
}

# Add legend
legend( locator(1), legend = c( "GBM.post.r1", "GBM.post.r2", "GBM.post.r3", "GBM.post.r4", "GBM.post.r5",
                                "GBM.post.r6", "GBM.post.r7", "GBM.post.r8", "GBM.post.r9", "GBM.post.r10"),
        col = col.plot,
        lty = 1,
        bty = "n"
)
abline( v = 181.8, col = "black", lty = 2 )

# All seems OK!