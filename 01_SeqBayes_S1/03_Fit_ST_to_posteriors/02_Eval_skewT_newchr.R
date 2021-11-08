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
load( "00_fitST_newchr/Rdata/ST.fitted.objects.RData" )
load( "00_fitST_newchr/Rdata/ST.fitted.dists.RData" )
load( "00_fitST_newchr/Rdata/post.divtimes.RData" )
divtimes_UPDGC <- divtimes 
rm( divtimes )
load( "00_fitST/Rdata/post.divtimes.RData" )

#-------#
# PLOTS #
#-------#
# Plot analytical ST VS Posterior updated geochr VS Posterior with calibs before sep2021
if( ! dir.exists( "02_MCMCtree_prior_ST_newchr" ) ){
  dir.create( "02_MCMCtree_prior_ST_newchr" )
}
if( ! dir.exists( "02_MCMCtree_prior_ST_newchr/plots" ) ){
  dir.create( "02_MCMCtree_prior_ST_newchr/plots" )
}

png( filename = paste( "02_MCMCtree_prior_ST_newchr/plots/00_Check_analitST_priorST_posterior_SMALL.png", sep = "" ),
     width = 1024, height = 768 )
par( mfrow = c(8,9), mai=c(0.1,0.1,0.1,0.1))
for( i in 1:length( colnames( divtimes ) ) ){
  
  # 1. Get a tmp_node with ST distributions
  tmp_node <- ST.fitted.objects[[ i ]]
  tracemem( tmp_node )
  # 2. Plot fitted ST for each node previously 
  #    computed using the values sampled during the MCMC
  
  # 2.1. Find limit axis
  max_x_st            <- round( max( density( rst( n = 1000, dp = tmp_node$dp.complete ) )$x ) + 0.5 )
  max_x_chain_UPDGC   <- round( max( density( divtimes_UPDGC[,i] )$x ) + 0.5 )
  max_x_chain         <- round( max( density( divtimes[,i] )$x ) + 0.5 )
  x_lim       <- max( max_x_chain_UPDGC, max_x_st, max_x_chain )
  max_y_st            <- round( max( density( rst( n = 1000, dp = tmp_node$dp.complete ) )$y ) + 0.5 )
  max_y_chain_UPDGC   <- round( max( density( divtimes_UPDGC[,i] )$y ) + 0.5 )
  max_y_chain         <- round( max( density( divtimes[,i] )$y ) + 0.5 )
  y_lim       <- max( max_y_chain_UPDGC, max_y_st, max_y_chain )
  # 2.2. Plot
  plot( density( divtimes_UPDGC[,i], adj = 1 ),
        xlim = c( 0, x_lim ), ylim = c( 0, y_lim ), #main = colnames( divtimes_UPDGC )[i], cex.main = 0.6,
        main=NULL, xlab = '', ylab = '', cex.axis = 0.7, mgp = c(2.5,0,0),
        #xaxt = "n", yaxt = "n", 
        tck = 0.01, col.main = "white" )
  title( colnames( divtimes_UPDGC )[i], line = -2, cex.main = 0.8, adj = 0.9 )
  lines( density( divtimes[,i], adj = 1 ),
         col = "green" )
  curve( dst( x, xi = tmp_node$dp[1], omega = tmp_node$dp[2],
              alpha = tmp_node$dp[3], nu = tmp_node$dp[4] ),
         from = 0, to = x_lim,
         n = 1e4, add = TRUE, col = "red" )

}

plot( x = 1, y = 1, col = "white", xaxt = "n", yaxt = "n", xlab = '', ylab = '' )
info.legend <- c( "Post. ages new geochr",
                  "Post. ages old geochr",
                  "ST new geochr" )
col.legend  <- c( "black", "green", "red" )
#coords.plot <- locator()
legend( x = 0.611875, y = 1.099012, legend = info.legend, col = col.legend,
        lty = 1, bty = 'n', cex = 0.7 )

dev.off()

# Repeat the same but extract individual plots too
if( ! dir.exists( "02_MCMCtree_prior_ST_newchr/plots/individual_plots" ) ){
  dir.create( "02_MCMCtree_prior_ST_newchr/plots/individual_plots" )
}
for( i in 1:length( colnames( divtimes_UPDGC ) ) ){
  
  # 1. Get a tmp_node with ST distributions
  tmp_node <- ST.fitted.objects[[ i ]]
  
  # 2. Plot fitted ST for each node previously 
  #    computed using the values sampled during the MCMC
  png( filename = paste( "02_MCMCtree_prior_ST_newchr/plots/individual_plots/Compare_newchr.vs.STnewchr_", colnames( divtimes )[i], ".png", sep = "" ),
       width = 1024, height = 768 )
  # 2.1. Find limit axis
  max_x_st    <- round( max( density( rst( n = 1000, dp = tmp_node$dp.complete ) )$x ) + 0.5 )
  max_x_chain <- round( max( density( divtimes_UPDGC[,i] )$x ) + 0.5 )
  x_lim       <- max( max_x_st, max_x_chain )
  max_y_st    <- round( max( density( rst( n = 1000, dp = tmp_node$dp.complete ) )$y ) + 0.5 )
  max_y_chain <- round( max( density( divtimes_UPDGC[,i] )$y ) + 0.5 )
  y_lim       <- max( max_y_st, max_y_chain )
  # 2.2. Plot
  plot( density( divtimes_UPDGC[,i], adj = 1 ),
        xlim = c( 0, x_lim ), ylim = c( 0, y_lim ), #main = colnames( divtimes_UPDGC )[i],
        xlab = '', ylab = '', col = "green",
        main = paste( colnames( divtimes_UPDGC )[i], " = ",
                      "ST(", paste0( round(tmp_node$dp, 2),
                                     collapse = "," ),
                      ")", sep = "" ) )
  curve( dst( x, xi = tmp_node$dp[1], omega = tmp_node$dp[2],
              alpha = tmp_node$dp[3], nu = tmp_node$dp[4] ),
         from = 0, to = x_lim,
         n = 1e4, add = TRUE, col = "red" )
  # 2.3. Legend 
  info.legend <- c( "Post. new geochr",
                    "ST new geochr" )
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
load( "00_fitST_newchr/Rdata/ST.fitted.objects.RData" )
load( "00_fitST_newchr/Rdata/ST.fitted.dists.RData" )
if( ! dir.exists( "02_MCMCtree_prior_ST_newchr/plots/ST_comparison" ) ){
  dir.create( "02_MCMCtree_prior_ST_newchr/plots/ST_comparison" )
}
for( i in 1:length( colnames( divtimes ) ) ){
  
  # 1. Get a tmp_node with ST distributions
  tmp_node     <- ST.fitted.objects[[ i ]]
  tmp_node_OLD <- ST.fitted.objects.old[[ i ]]
  
  # 2. Plot fitted ST for each node previously 
  #    computed using the values sampled during the MCMC
  png( filename = paste( "02_MCMCtree_prior_ST_newchr/plots/ST_comparison/Compare_OLDvsNEW_ST_", colnames( divtimes )[i], ".png", sep = "" ),
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
png( filename = paste( "02_MCMCtree_prior_ST_newchr/plots/Compare_OLDvsNEW_ST_SMALL.png", sep = "" ),
     width = 1024, height = 768 )
#par( mfrow = c(8,9), mai=c(0,0,0,0))
par( mfrow = c(8,9), mai=c(0.1,0.1,0.1,0.1))
for( i in 1:length( colnames( divtimes ) ) ){
  
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
