#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------------#
# SET WORKING DIRECTORY #
#-----------------------#
library( rstudioapi ) 
# Get the path to current open R script and find main dir
path_to_file <- getActiveDocumentContext()$path
wd <- paste( dirname( path_to_file ), "/", sep = "" )
setwd( wd )

#--------------#
# LOAD LIBRARY #
#--------------#
library( sn )

#-----------#
# LOAD DATA #
#-----------#
# Load objects previously created with the fitted ST dists and 
# the posterior divtimes
load( "00_fitST_new/Rdata/ST.fitted.objects.RData" )
load( "00_fitST_new/Rdata/ST.fitted.dists.RData" )
load( "00_fitST_new/Rdata/post.divtimes.RData" )

# Load prior from MCMC when using ST dists as calibs
prior.ST.run1 <- read.table( "02_MCMCtree_prior_ST/MCMCtree/mcmc1/mcmc.txt",
                             header = T, sep = "\t" )
prior.ST.run2 <- read.table( "02_MCMCtree_prior_ST/MCMCtree/mcmc2/mcmc.txt",
                             header = T, sep = "\t" )
# Get parameters and divtimes
divtimes.prior.ST.run1 <- prior.ST.run1[,-c(1,(dim(prior.ST.run1)[2]) )]
divtimes.prior.ST.run2 <- prior.ST.run2[,-c(1,(dim(prior.ST.run2)[2]) )]
divtimes.prior.ST  <- rbind( divtimes.prior.ST.run1, divtimes.prior.ST.run2 )

mean_est_divt_prior.ST  <- apply( X = divtimes.prior.ST, MARGIN = 2, FUN = mean )
quant_est_divt_prior.ST <- apply( X = divtimes.prior.ST, MARGIN = 2, FUN = quantile, probs = c(0.025,0.975) )
quant_est_divt_prior.ST <- t( quant_est_divt_prior.ST )
all.equal( quant_est_divt_prior.ST[1,], quantile( divtimes.prior.ST[,1], probs = c( 0.025, 0.975 ) ) )

#-------#
# PLOTS #
#-------#
# Plot analytical ST VS prior with analytical ST VS Posterior
# png( filename = paste( "02_MCMCtree_prior_ST_new/plots/00_Check_analitST_priorST_posterior.png", sep = "" ),
png( filename = paste( "02_MCMCtree_prior_ST_new/plots/00_Check_analitST_priorST_posterior_SMALL.png", sep = "" ),
     width = 1024, height = 768 )
# pdf( file = paste( "02_MCMCtree_prior_ST_new/plots/00_Check_analitST_priorST_posterior_SMALL.pdf", sep = "" ) )
par( mfrow = c(8,9), mai=c(0.1,0.1,0.1,0.1))
for( i in 1:length( colnames( divtimes.prior.ST ) ) ){
  
  # 1. Get a tmp_node with ST distributions
  tmp_node <- ST.fitted.objects[[ i ]]
  tracemem( tmp_node )
  # 2. Plot fitted ST for each node previously 
  #    computed using the values sampled during the MCMC
  
  # 2.1. Find limit axis
  max_x_st            <- round( max( density( rst( n = 1000, dp = tmp_node$dp.complete ) )$x ) + 0.5 )
  max_x_chain         <- round( max( density( divtimes[,i] )$x ) + 0.5 )
  max_x_chain_priorST <- round( max( density( divtimes.prior.ST[,i] )$x ) + 0.5 )
  x_lim       <- max( max_x_chain, max_x_st, max_x_chain_priorST )
  max_y_st            <- round( max( density( rst( n = 1000, dp = tmp_node$dp.complete ) )$y ) + 0.5 )
  max_y_chain         <- round( max( density( divtimes[,i] )$y ) + 0.5 )
  max_y_chain_priorST <- round( max( density( divtimes.prior.ST[,i] )$y ) + 0.5 )
  y_lim       <- max( max_y_chain, max_y_st, max_y_chain_priorST )
  # 2.2. Plot
  plot( density( divtimes[,i], adj = 1 ),
        xlim = c( 0, x_lim ), ylim = c( 0, y_lim ), #main = colnames( divtimes )[i], cex.main = 0.6,
        main=NULL, xlab = '', ylab = '', cex.axis = 0.7, mgp = c(2.5,0,0),
        #xaxt = "n", yaxt = "n", 
        tck = 0.01, col.main = "white" )
  title( colnames( divtimes )[i], line = -2, cex.main = 0.8, adj = 0.9 )
  lines( density( divtimes.prior.ST[,i], adj = 1 ),
         col = "green" )
  curve( dst( x, xi = tmp_node$dp[1], omega = tmp_node$dp[2],
              alpha = tmp_node$dp[3], nu = tmp_node$dp[4] ),
         from = 0, to = x_lim,
         n = 1e4, add = TRUE, col = "red" )

}

plot( x = 1, y = 1, col = "white", xaxt = "n", yaxt = "n", xlab = '', ylab = '' )
info.legend <- c( "Post. ages with SB calibs",
                  "Prior ages with 71 ST calibs",
                  "Analytical ST calibrations" )
col.legend  <- c( "black", "green", "red" )
#coords.plot <- locator()
legend( x = 0.611875, y = 1.099012, legend = info.legend, col = col.legend,
        lty = 1, bty = 'n', cex = 0.7 )

dev.off()


# Plot analytical ST VS prior with analytical ST
png( filename = paste( "02_MCMCtree_prior_ST_new/plots/00_Check_analitST_priorST.png", sep = "" ),
     width = 1024, height = 768 )
par( mfrow = c(8,9), mai=c(0,0,0,0))
for( i in 1:length( colnames( divtimes.prior.ST ) ) ){
  
  # 1. Get a tmp_node with ST distributions
  tmp_node <- ST.fitted.objects[[ i ]]
  
  # 2. Plot fitted ST for each node previously 
  #    computed using the values sampled during the MCMC
  
  # 2.1. Find limit axis
  max_x_st            <- round( max( density( rst( n = 1000, dp = tmp_node$dp.complete ) )$x ) + 0.5 )
  max_x_chain_priorST <- round( max( density( divtimes.prior.ST[,i] )$x ) + 0.5 )
  x_lim       <- max( max_x_st, max_x_chain_priorST )
  max_y_st            <- round( max( density( rst( n = 1000, dp = tmp_node$dp.complete ) )$y ) + 0.5 )
  max_y_chain_priorST <- round( max( density( divtimes.prior.ST[,i] )$y ) + 0.5 )
  y_lim       <- max( max_y_st, max_y_chain_priorST )
  # 2.2. Plot
  plot( density( divtimes.prior.ST[,i], adj = 1 ),
        xlim = c( 0, x_lim ), ylim = c( 0, y_lim ), main = colnames( divtimes.prior.ST )[i],
        xlab = '', ylab = '', col = "green" )
  curve( dst( x, xi = tmp_node$dp[1], omega = tmp_node$dp[2],
              alpha = tmp_node$dp[3], nu = tmp_node$dp[4] ),
         from = 0, to = x_lim,
         n = 1e4, add = TRUE, col = "red" )
  
}

plot( x = 1, y = 1, col = "white" )
info.legend <- c( "Prior ages using 71 ST calibs",
                  "Analytical ST calibrations" )
col.legend  <- c( "green", "red" )
#coords.plot <- locator()
legend( x = 0.611875, y = 1.099012, legend = info.legend, col = col.legend,
        lty = 1, bty = 'n', cex = 0.8 )

dev.off()

# Repeat the same but extract individual plots too
for( i in 1:length( colnames( divtimes.prior.ST ) ) ){
  
  # 1. Get a tmp_node with ST distributions
  tmp_node <- ST.fitted.objects[[ i ]]
  
  # 2. Plot fitted ST for each node previously 
  #    computed using the values sampled during the MCMC
  png( filename = paste( "02_MCMCtree_prior_ST_new/plots/individual_plots/Compare_priorST.vs.analytST_", colnames( divtimes )[i], ".png", sep = "" ),
       width = 1024, height = 768 )
  # 2.1. Find limit axis
  max_x_st            <- round( max( density( rst( n = 1000, dp = tmp_node$dp.complete ) )$x ) + 0.5 )
  max_x_chain_priorST <- round( max( density( divtimes.prior.ST[,i] )$x ) + 0.5 )
  x_lim       <- max( max_x_st, max_x_chain_priorST )
  max_y_st            <- round( max( density( rst( n = 1000, dp = tmp_node$dp.complete ) )$y ) + 0.5 )
  max_y_chain_priorST <- round( max( density( divtimes.prior.ST[,i] )$y ) + 0.5 )
  y_lim       <- max( max_y_st, max_y_chain_priorST )
  # 2.2. Plot
  plot( density( divtimes.prior.ST[,i], adj = 1 ),
        xlim = c( 0, x_lim ), ylim = c( 0, y_lim ), #main = colnames( divtimes.prior.ST )[i],
        xlab = '', ylab = '', col = "green",
        main = paste( colnames( divtimes )[i], " = ",
                      "ST(", paste0( round(tmp_node$dp, 2),
                                     collapse = "," ),
                      ")", sep = "" ) )
  curve( dst( x, xi = tmp_node$dp[1], omega = tmp_node$dp[2],
              alpha = tmp_node$dp[3], nu = tmp_node$dp[4] ),
         from = 0, to = x_lim,
         n = 1e4, add = TRUE, col = "red" )
  # 2.3. Legend 
  info.legend <- c( "Est. prior times using 71 ST calibs",
                    "Analytical ST calibs" )
  col.legend  <- c( "green", "red" )
  legend( "topright", legend = info.legend, col = col.legend,
          lty = 1, bty = 'n', cex = 0.8 )
  # 2.3. Close graph
  dev.off()
  
}


# Plot OLD vs NEW ST
# Repeat the same but extract individual plots too
load( "00_fitST/Rdata/ST.fitted.objects.RData" )
ST.fitted.objects.old <- ST.fitted.objects
load( "00_fitST/Rdata/ST.fitted.dists.RData" )
ST.fitted.dists.old <- ST.fitted.dists
load( "00_fitST_new/Rdata/ST.fitted.objects.RData" )
load( "00_fitST_new/Rdata/ST.fitted.dists.RData" )
for( i in 1:length( colnames( divtimes.prior.ST ) ) ){
  
  # 1. Get a tmp_node with ST distributions
  tmp_node     <- ST.fitted.objects[[ i ]]
  tmp_node_OLD <- ST.fitted.objects.old[[ i ]]
  
  # 2. Plot fitted ST for each node previously 
  #    computed using the values sampled during the MCMC
  png( filename = paste( "02_MCMCtree_prior_ST_new/plots/ST_comparison/Compare_OLDvsNEW_ST_", colnames( divtimes )[i], ".png", sep = "" ),
       width = 1024, height = 768 )
  # 2.1. Find limit axis
  max_x_st  <- round( max( density( rst( n = 1000, dp = tmp_node$dp.complete ) )$x ) + 0.5 )
  max_x_OLD <- round( max( density( rst( n = 1000, dp = tmp_node_OLD$dp.complete ) )$x ) + 0.5 )
  x_lim       <- max( max_x_st, max_x_OLD )
  max_y_st    <- round( max( density( rst( n = 1000, dp = tmp_node$dp.complete ) )$y ) + 0.5 )
  max_y_OLD   <- round( max( density( rst( n = 1000, dp = tmp_node_OLD$dp.complete ) )$y ) + 0.5 )
  y_lim       <- max( max_y_st, max_y_OLD )
  # 2.2. Plot
  curve( dst( x, xi = tmp_node$dp[1], omega = tmp_node$dp[2],
              alpha = tmp_node$dp[3], nu = tmp_node$dp[4] ),
         from = 0, to = x_lim,
         n = 1e4, col = "green",
         xlab = '', ylab = '', 
         main = paste( colnames( divtimes )[i], " = ",
                       "ST(", paste0( round(tmp_node$dp, 2),
                                      collapse = "," ),
                       ")", sep = "" ) )
  curve( dst( x, xi = tmp_node_OLD$dp[1], omega = tmp_node_OLD$dp[2],
              alpha = tmp_node_OLD$dp[3], nu = tmp_node_OLD$dp[4] ),
         from = 0, to = x_lim,
         n = 1e4, add = TRUE, col = "red" )
  # 2.3. Legend 
  info.legend <- c( "NEW ST",
                    "OLD ST" )
  col.legend  <- c( "green", "red" )
  legend( "topright", legend = info.legend, col = col.legend,
          lty = 1, bty = 'n', cex = 0.8 )
  # 2.3. Close graph
  dev.off()
  
}


# Plot analytical ST VS prior with analytical ST
png( filename = paste( "02_MCMCtree_prior_ST_new/plots/Compare_OLDvsNEW_ST.png", sep = "" ),
     width = 1024, height = 768 )
par( mfrow = c(8,9), mai=c(0,0,0,0))
for( i in 1:length( colnames( divtimes.prior.ST ) ) ){
  
  # 1. Get a tmp_node with ST distributions
  tmp_node     <- ST.fitted.objects[[ i ]]
  tmp_node_OLD <- ST.fitted.objects.old[[ i ]]
  # 2. Plot fitted ST for each node previously 
  #    computed using the values sampled during the MCMC
  
  # 2.1. Find limit axis
  max_x_st  <- round( max( density( rst( n = 1000, dp = tmp_node$dp.complete ) )$x ) + 0.5 )
  max_x_OLD <- round( max( density( rst( n = 1000, dp = tmp_node_OLD$dp.complete ) )$x ) + 0.5 )
  x_lim       <- max( max_x_st, max_x_OLD )
  max_y_st    <- round( max( density( rst( n = 1000, dp = tmp_node$dp.complete ) )$y ) + 0.5 )
  max_y_OLD   <- round( max( density( rst( n = 1000, dp = tmp_node_OLD$dp.complete ) )$y ) + 0.5 )
  y_lim       <- max( max_y_st, max_y_OLD )
  # 2.2. Plot
  curve( dst( x, xi = tmp_node$dp[1], omega = tmp_node$dp[2],
                    alpha = tmp_node$dp[3], nu = tmp_node$dp[4] ),
               from = 0, to = x_lim,
               n = 1e4, col = "green" )
  curve( dst( x, xi = tmp_node_OLD$dp[1], omega = tmp_node_OLD$dp[2],
              alpha = tmp_node_OLD$dp[3], nu = tmp_node_OLD$dp[4] ),
         from = 0, to = x_lim,
         n = 1e4, add = TRUE, col = "red" )
  
}

plot( x = 1, y = 1, col = "white" )
info.legend <- c( "NEW",
                  "OLD" )
col.legend  <- c( "green", "red" )
#coords.plot <- locator()
legend( x = 0.611875, y = 1.099012, legend = info.legend, col = col.legend,
        lty = 1, bty = 'n', cex = 0.8 )

dev.off()