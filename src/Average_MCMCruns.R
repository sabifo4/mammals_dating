#==========#
# FUNCTION #
#==========#--------------------------------------------------------#
# Average branch lenghts and node ages in object "phy"
# Object "phy" is previously generated and it contains one 
# entry for each MCMC that ran in MCMCtree
# Argument: 
#     phy: List, length equals the number of MCMC runs, 
#          so there is one "phy" object in each entry of the list. 
average.runs <- function ( phy ){
  
  # Generate global vars
  num.runs      <- length( phy )
  num.rows.tree <- length( phy[[ 1 ]]$apePhy$edge.length )
  tmp.mat.bl    <- matrix( 0, ncol = num.runs, nrow = num.rows.tree )
  tmp.mat.mean.ag <- tmp.mat.upb.ag <- tmp.mat.lowb.ag <- 
    matrix( 0, ncol = num.runs, nrow = dim( phy[[ 1 ]]$nodeAges )[1] )
  for( i in 1:num.runs){
    tmp.mat.bl[,i]      <- phy[[ i ]]$apePhy$edge.length
    tmp.mat.mean.ag[,i] <- phy[[ i ]]$nodeAges[,1]
    tmp.mat.lowb.ag[,i] <- phy[[ i ]]$nodeAges[,2]
    tmp.mat.upb.ag[,i]  <- phy[[ i ]]$nodeAges[,3]
  }
  
  # Get mean branch lengths 
  mean.bl <- apply( tmp.mat.bl, 1, mean )
  
  # Get mean node ages, upper and low bounds
  mean.nodes <- apply( tmp.mat.mean.ag, 1, mean )
  lowbound   <- apply( tmp.mat.lowb.ag, 1, mean )
  upbound    <- apply( tmp.mat.upb.ag, 1, mean )
  nodeAges   <- cbind( mean.nodes, lowbound, upbound )
  colnames( nodeAges ) <- colnames( phy[[ 1 ]]$nodeAges )
  rownames( nodeAges ) <- rownames( phy[[ 1 ]]$nodeAges )
  
  # Generate phy with average nodeAges and branch lengths 
  outPhy                    <- phy[[ 1 ]]
  outPhy$apePhy$edge.length <- mean.bl
  outPhy$nodeAges           <- nodeAges
  
  # Return outPhy
  return( outPhy )
  
}
