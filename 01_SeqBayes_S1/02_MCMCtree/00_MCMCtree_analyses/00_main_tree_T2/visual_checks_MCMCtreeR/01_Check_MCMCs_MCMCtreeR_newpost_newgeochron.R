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
wd_trees <- gsub( pattern = "visual_checks_MCMCtreeR/", replacement = "", x = wd )
source( "Functions_plots_MCMCtreeR.R" )

#-------------#
# START TASKS #
#-------------#

# 1. Load trees and mcmc files for master tree with the updated prior calibrations 
#    as of September 2021
## GBM - posterior NEW
phy.post.GBM <- mcmc.post.GBM <- vector( mode = "list", 4 )
for( i in seq(1:length( phy.post.GBM ) ) ){
        
        cat( "Parsing files for run ", i, "...\n" )
        phy.post.GBM[[ i ]] <-  MCMCtreeR::readMCMCtree( paste( wd_trees, "02_MCMCtree_posterior_newchrono/run",
                                                                i, "/mcmctree_GBM/FigTree.tre", sep = "" ),
                                                         from.file = TRUE )
        mcmc.post.GBM[[ i ]] <- read.table( paste( wd_trees, "02_MCMCtree_posterior_newchrono/run",
                                                   i, "/mcmctree_GBM/mcmc.txt", sep = "" ),
                                            sep = "\t", header = T, stringsAsFactors = F )
}

## GBM - posterior OLD (this is the analysis inside `01_MCMCtree_posterior`, which also uses 
##       alignment without 11 genes common in second data set but the set of priors used is 
##       not updated as of Sep 2021)
phy.oldpost.GBM <- mcmc.oldpost.GBM <- vector( mode = "list",  length( c(3,4,7,8,9,11,13,16) ) )
count <- 0 
for( i in c(3,4,7,8,9,11,13,16) ){
        
        cat( "Parsing files for run ", i, "...\n" )
        count <- count + 1
        phy.oldpost.GBM[[ count ]] <-  MCMCtreeR::readMCMCtree( paste( wd_trees, "01_MCMCtree_posterior/02_atlantogenata_tarver2016/run",
                                                                       i, "/mcmctree_GBM/FigTree.tre", sep = "" ),
                                                            from.file = TRUE )
        mcmc.oldpost.GBM[[ count ]] <- read.table( paste( wd_trees, "01_MCMCtree_posterior/02_atlantogenata_tarver2016/run",
                                                          i, "/mcmctree_GBM/mcmc.txt", sep = "" ),
                                                   sep = "\t", header = T, stringsAsFactors = F )
}

# 2. Create objects in APE format 
## GBM - posterior NEW
phy.post.GBM.obj <- phy.post.GBM.edge <- mcmc.chain.GBM <- vector( mode = "list", 4 )
for( i in seq(1, length( phy.post.GBM.obj ) ) ){
        cat( "Parsing files from run ", i, "...\n" )
        phy.post.GBM.obj[[ i ]]  <- phy.post.GBM[[ i ]]$apePhy
        phy.post.GBM.edge[[ i ]] <- phy.post.GBM.obj[[ i ]]$edge
        mcmc.chain.GBM[[ i ]]    <- mcmc.post.GBM[[ i ]]
}

## GBM - posterior OLD 
phy.oldpost.GBM.obj <- phy.oldpost.GBM.edge <- mcmc.chain.oldpost.GBM <- vector( mode = "list", 8 )
for( i in seq( 1, length( phy.oldpost.GBM.obj ) ) ){
        cat( "Parsing files from run ", i, "...\n" )
        phy.oldpost.GBM.obj[[ i ]]     <- phy.oldpost.GBM[[ i ]]$apePhy
        phy.oldpost.GBM.edge[[ i ]]    <- phy.oldpost.GBM.obj[[ i ]]$edge
        mcmc.chain.oldpost.GBM[[ i ]]  <- mcmc.oldpost.GBM[[ i ]]
}

# 3. Create PHY object with GBM and ILN objects 
phy.all <- vector( mode = "list", length = 12 )
names( phy.all ) <- c( paste( "GBM.post.", 1:4, sep = "" ),
                       paste( "GBM.oldpost.", 1:8, sep = "" ) )
k <- 0
for( i in seq( 1:length( phy.all ) ) ){
        
        if( i < 5 ){
                print( i )
                cat( "Parsing GBM.post files from run ", i, "...\n" )
                phy.all[[ i ]] <- phy.post.GBM.obj[[ i ]]
        }
        else{
                k <- k + 1
                print( i )
                cat( "Parsing GBM.oldpost files from run ", k, "...\n" )
                phy.all[[ i ]] <- phy.oldpost.GBM.obj[[ k ]] 
        }
}


# 4. Extract ages with node age posteriors from column 2
## GBM - posterior
mcmc.node.ages.post.GBM <- all.nodes.post.GBM <- vector( mode = "list", 4 )
for( i in 1:length( mcmc.node.ages.post.GBM ) ){
        mcmc.node.ages.post.GBM[[ i ]] <- mcmc.chain.GBM[[ i ]][, 2:Ntip( phy.all[[ i ]] )] 
        all.nodes.post.GBM[[ i ]]      <- as.numeric( gsub( "t_n", "",
                                                            colnames( mcmc.node.ages.post.GBM[[ i ]] ) ) )
}

## GBM - oldpost
mcmc.node.ages.oldpost.GBM <- all.nodes.oldpost.GBM <- vector( mode = "list", 8 )
j <- 4
for( i in 1:8 ){
        j <- j + 1
        cat( "Parsing entry ", j, "in phy.all\n" )
        mcmc.node.ages.oldpost.GBM[[ i ]] <- mcmc.chain.oldpost.GBM[[ i ]][, 2:Ntip( phy.all[[ j ]] )] 
        all.nodes.oldpost.GBM[[ i ]]      <- as.numeric( gsub( "t_n", "",
                                                             colnames( mcmc.node.ages.oldpost.GBM[[ i ]] ) ) )
}

# 5. Create a vector of names for each list element as  
#    internal nodes from APE tree, using phy$edge object 
## GBM - posterior
node.ages.names.post.GBM <- vector( mode = "list", 4 )
for( i in 1:length( node.ages.names.post.GBM ) ){
        node.ages.names.post.GBM[[ i ]] <- c( Ntip( phy.all[[ i ]] ) + 1,
                                              phy.post.GBM.edge[[ i ]][ which( phy.post.GBM.edge[[ i ]][,2] > Ntip( phy.all[[ i ]] ) ), 2 ] )
}

## GBM - oldpost
node.ages.names.oldpost.GBM <- vector( mode = "list", 8 )
j <- 4
for( i in 1:length( node.ages.names.oldpost.GBM ) ){
        j <- j + 1
        node.ages.names.oldpost.GBM[[ i ]] <- c( Ntip( phy.all[[ j ]] ) + 1,
                                              phy.oldpost.GBM.edge[[ i ]][ which( phy.oldpost.GBM.edge[[ i ]][,2] > Ntip( phy.all[[ j ]] ) ), 2 ] )
}

# 6. Find where each posterior node age appears in APE 
#    edge object.
## GBM - posterior
match.nodes.post.GBM <- vector( mode = "list", 4 )
for( i in 1:length( match.nodes.post.GBM ) ){
        match.nodes.post.GBM[[ i ]] <- match( all.nodes.post.GBM[[ i ]],
                                              as.numeric( node.ages.names.post.GBM[[ i ]] ) )
}

## GBM - oldpost
match.nodes.oldpost.GBM <- vector( mode = "list", 8 )
for( i in 1:length( match.nodes.oldpost.GBM ) ){
        match.nodes.oldpost.GBM[[ i ]] <- match( all.nodes.oldpost.GBM[[ i ]],
                                              as.numeric( node.ages.names.oldpost.GBM[[ i ]] ) )
}


# 7. Create a list extracting the info from the mcmc
#    chain in APE node order 
## GBM - posterior
node.ages.post.GBM <- vector( mode = "list", 4 )
for( i in 1:length( node.ages.post.GBM ) ){
        node.ages.post.GBM[[ i ]] <- lapply( match.nodes.post.GBM[[ i ]],
                                             function( uu ) mcmc.node.ages.post.GBM[[ i ]][, uu ] )
}

## GBM - oldpost 
node.ages.oldpost.GBM <- vector( mode = "list", 8 )
for( i in 1:length( node.ages.oldpost.GBM ) ){
        node.ages.oldpost.GBM[[ i ]] <- lapply( match.nodes.oldpost.GBM[[ i ]],
                                             function( uu ) mcmc.node.ages.oldpost.GBM[[ i ]][, uu ] )
}


# 8. Name each element in list 
## GBM - posterior
for( i in 1:length( node.ages.names.post.GBM ) ){
        names( node.ages.post.GBM[[ i ]] ) <- node.ages.names.post.GBM[[ i ]]
        
}

## GBM - oldpost 
for( i in 1:length( node.ages.names.oldpost.GBM ) ){
        names( node.ages.oldpost.GBM[[ i ]] ) <- node.ages.names.oldpost.GBM[[ i ]]
        
}


# 9. Create node.ages.all object to get all GBM and ILN objects
node.ages.all <- vector( mode = "list", length = 12 )
k <- 0
for( i in 1:length( node.ages.all ) ){
        
        if( i < 5 ){
                node.ages.all[[ i ]] <- node.ages.post.GBM[[ i ]]
        }
        else{
                k <- k + 1
                node.ages.all[[ i ]] <- node.ages.oldpost.GBM[[ k ]]
        }
}

# 10. Get means for edges
## GBM - posterior
# Initialize matrix to get the edges of the 16 runs 
phy.matrix.post.GBM.edges <- matrix( 0, nrow = 4, ncol = length( phy.all[[1]]$edge.length ) )
for( i in 1:4 ){
        phy.matrix.post.GBM.edges[i,] <- phy.all[[ i ]]$edge.length
}
# Get mean of each edge. Also generate one mean without chains 1,2,5,6,10,12,14,15 -- according to visualization
# The latter object has been generated and added here after analysing the plots that appear later
# in the script
phy.mean.edge.post.GBM      <- apply( X = phy.matrix.post.GBM.edges, 2, mean )
phy.mean.edge.post.GBM.filt <- apply( X = phy.matrix.post.GBM.edges[c(1,3:4),], 2, mean ) 
# Create copy of ape object and then replace the edge.length values accordingly
phy.mean.post.GBM                   <- phy.mean.post.GBM.filt <- phy.all[[1]]
phy.mean.post.GBM$edge.length       <- phy.mean.edge.post.GBM
phy.mean.post.GBM.filt$edge.length  <- phy.mean.edge.post.GBM.filt

## GBM - oldpost
# Intialize matrix to get the edges of the 7 runs 
phy.matrix.oldpost.GBM.edges <- matrix( 0, nrow = 8, ncol = length( phy.all[[5]]$edge.length ) )
j <- 4
for( i in 1:8 ){
        j <- j + 1
        phy.matrix.oldpost.GBM.edges[i,] <- phy.all[[ j ]]$edge.length
}
# Get mean of each edge.
phy.mean.edge.oldpost.GBM      <- apply( X = phy.matrix.oldpost.GBM.edges, 2, mean )
# Create copy of ape object and then replace the edge.length values accordingly
phy.mean.oldpost.GBM                   <- phy.all[[5]]
phy.mean.oldpost.GBM$edge.length       <- phy.mean.edge.oldpost.GBM

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
#                     although it can be less, e.g., "oldpost.GBM".
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
tmp <- rep( 0, 4*length( unlist( node.ages.post.GBM[[ 1 ]][1] ) ) )
# b) Loop from node 1 to 71 and then go through the corresponding values sampled for 
#    this node for each run. Use tmp vector to save values
node.mean.dist.ages.post.GBM <- get_mean_dist_ages( len = length( node.ages.post.GBM[[ 1 ]] ),
                                                    chains = 1:4,
                                                    node.ages.list = node.ages.post.GBM,
                                                    tmp = tmp )

# c) Now we will create a list with 71 node entries and we will concatenate the sampled values 
#    and I will just pick those chains that seem not to have problems -- visually, they seemed 
#    all OK, but a check with quantiles seems to point that chain 2 had issues, i.e., 1,3,4.
tmp <- rep( 0, 3*length( unlist( node.ages.post.GBM[[ 1 ]][1] ) ) )
node.mean.dist.ages.post.GBM.filt <- get_mean_dist_ages( len = length( node.ages.post.GBM[[ 1 ]] ),
                                                         chains = c(1,3,4),
                                                         node.ages.list = node.ages.post.GBM,
                                                         tmp = tmp )

# d) Now the same for GBM.oldpost 
tmp <- rep( 0, 8*length( unlist( node.ages.oldpost.GBM[[ 1 ]][1] ) ) )
node.mean.dist.ages.oldpost.GBM <- get_mean_dist_ages( len = length( node.ages.oldpost.GBM[[ 1 ]] ),
                                                     chains = 1:8,
                                                     node.ages.list = node.ages.oldpost.GBM,
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

#====================================#
# 0. PLOT: Evaluate posterior chains #
#====================================#
# Saved as 
#   plots/03_NEW-UPDGC-Explore_GBM_posterior_runs.pdf
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
col.plot <- c( "purple", "black", "darkgreen", "red" )
count <- 1 # Purple is for chain plotted first
for( i in c(2:4) ){
        count <- count + 1
        add.extra.dists( phy = phy.all[[ i ]], num.models = 1, last.plot = last.plot,
                         node.ages =  node.ages.post.GBM[[ i ]], plot.type = "distributions",
                         time.correction = 100,
                         density.col = "white", density.border.col = col.plot[ count ],
                         distribution.height = 0.8, transparency = transparency
        )
}

# Add legend
legend( locator(1), legend = c( "GBM.post.r1", "GBM.post.r2", "GBM.post.r3", "GBM.post.r4" ),
        col = col.plot[1:4],
        lty = 1,
        bty = "n"
)
abline( v = 179.5, col = "black", lty = 2 )

#=============================#
# 1. PLOT: Filtered posterior #
#=============================#
## THE NEXT LINES CAN BE RUN AT ONCE
## THE OUTPUT FILE IS SAVED AS 
## This filtering is based on evaluating the 97.5%
## and 2.5% quantiles, not visual.
##    "plots/02_NEW-Filtered_chains-1.3.4.pdf"
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
col.plot <- c( "purple", "black", "darkgreen" )
count <- 1 # Purple is for chain plotted first
for( i in c(3,4) ){
        count <- count + 1
        add.extra.dists( phy = phy.all[[ i ]], num.models = 1, last.plot = last.plot,
                         node.ages =  node.ages.post.GBM[[ i ]], plot.type = "distributions",
                         time.correction = 100,
                         density.col = "white", density.border.col = col.plot[ count ],
                         distribution.height = 0.8, transparency = transparency
        )
}

# Add legend
legend( locator(1), legend = c( "GBM.post.r1", "GBM.post.r3", "GBM.post.r4" ),
        col = col.plot[1:3],
        lty = 1,
        bty = "n"
)
abline( v = 179.5, col = "black", lty = 2 )


#========================================#
# 2. PLOT: Plot post-old VS post-updated #
#========================================#
## THE NEXT LINES CAN BE RUN AT ONCE
## THE OUTPUT FILE IS SAVED AS 
##    "03_NEW-UPDGC-mean.postUPDGCVSmean.post.pdf"
# Set transparent colour
transp.col <- adjustcolor( col = "blue", alpha.f = 0.3 )
last.plot  <- mcmc.tree.plot.RETPLOT( phy = phy.mean.post.GBM.filt, xlim.scale = c(-20,300),
                                      node.ages = node.mean.dist.ages.post.GBM.filt, 
                                      show.tip.label = TRUE,
                                      analysis.type = "user", cex.tips = 0.8,
                                      time.correction = 100, 
                                      scale.res = c( "Eon", "Period" ), plot.type = "distributions",
                                      cex.age = 0.6, cex.labels = 0.8, relative.height = 0.08, 
                                      col.tree = "grey40", no.margin = TRUE,
                                      density.col = transp.col, add.time.scale = TRUE, grey.bars = FALSE
)

# Add now the other distribution
col.mod = "grey"
transparency = 0.3 
col.plot <- c( "darkgreen" )
add.extra.dists( phy = phy.mean.oldpost.GBM, num.models = 1, last.plot = last.plot,
                 node.ages =  node.mean.dist.ages.oldpost.GBM, plot.type = "distributions",
                 time.correction = 100,
                 density.col = "white", density.border.col = col.plot,
                 distribution.height = 0.8, transparency = transparency
)

# Add legend
legend( locator(1), legend = c( "mean.postUPDGC", "mean.post" ),
        col = c( "blue", "darkgreen" ),
        lty = 1,
        bty = "n"
)
abline( v = 183, col = "black", lty = 2 )

#--------------------#
# GENERATE R OBJECTS #
#--------------------#
# 1. Generate list object that can be used to plot
phy.all.UPDGC            <- phy.all
node.ages.post.GBM.UPDGC <- node.ages.post.GBM
save( phy.all.UPDGC, file = "R_objects/phy.all.UPDGC.Rdata" )
save( node.ages.post.GBM.UPDGC, file = "R_objects/node.ages.post.GBM.UPDGC.Rdata" )

# 2. Generate list object with the mean ages of the posterior after 
#    it has been filtered (it keeps only chains without conflict)
phy.mean.post.GBM.UPDGC.filt <- phy.mean.post.GBM.filt
save( phy.mean.post.GBM.UPDGC.filt, file = "R_objects/phy.mean.post.GBM.UPDGC.filt.RData" )
