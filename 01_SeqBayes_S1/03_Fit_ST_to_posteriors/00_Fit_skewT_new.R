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
wd_mcmctree <- gsub( pattern = "03_Fit_ST_to_posteriors/", replacement = "", x = wd )

#--------------#
# LOAD LIBRARY #
#--------------#
library( sn )

#---------------------------#
# LOAD AND PARSE DATA - NEW #
#---------------------------#
## POSTERIOR - GBM
# 1. Load files and get parameters
run3    <- read.table( paste( wd_mcmctree,
                                  "02_MCMCtree/00_MCMCtree_analyses/00_main_tree_T2/01_MCMCtree_posterior/02_atlantogenata_tarver2016/run3/mcmctree_GBM/mcmc.txt",
                                  sep = "" ), header = T, sep = "\t" )
run4    <- read.table( paste( wd_mcmctree,
                                  "02_MCMCtree/00_MCMCtree_analyses/00_main_tree_T2/01_MCMCtree_posterior/02_atlantogenata_tarver2016/run4/mcmctree_GBM/mcmc.txt",
                                  sep = "" ), header = T, sep = "\t" )
run7    <- read.table( paste( wd_mcmctree,
                                  "02_MCMCtree/00_MCMCtree_analyses/00_main_tree_T2/01_MCMCtree_posterior/02_atlantogenata_tarver2016/run7/mcmctree_GBM/mcmc.txt",
                                  sep = "" ), header = T, sep = "\t" )
run8    <- read.table( paste( wd_mcmctree,
                                  "02_MCMCtree/00_MCMCtree_analyses/00_main_tree_T2/01_MCMCtree_posterior/02_atlantogenata_tarver2016/run8/mcmctree_GBM/mcmc.txt",
                                  sep = "" ), header = T, sep = "\t" )
run9    <- read.table( paste( wd_mcmctree,
                                  "02_MCMCtree/00_MCMCtree_analyses/00_main_tree_T2/01_MCMCtree_posterior/02_atlantogenata_tarver2016/run9/mcmctree_GBM/mcmc.txt",
                                  sep = "" ), header = T, sep = "\t" )
run11    <- read.table( paste( wd_mcmctree,
                                   "02_MCMCtree/00_MCMCtree_analyses/00_main_tree_T2/01_MCMCtree_posterior/02_atlantogenata_tarver2016/run11/mcmctree_GBM/mcmc.txt",
                                   sep = "" ), header = T, sep = "\t" )
run13   <- read.table( paste( wd_mcmctree,
                                  "02_MCMCtree/00_MCMCtree_analyses/00_main_tree_T2/01_MCMCtree_posterior/02_atlantogenata_tarver2016/run13/mcmctree_GBM/mcmc.txt",
                                  sep = "" ), header = T, sep = "\t" )
run16   <- read.table( paste( wd_mcmctree,
                                  "02_MCMCtree/00_MCMCtree_analyses/00_main_tree_T2/01_MCMCtree_posterior/02_atlantogenata_tarver2016/run16/mcmctree_GBM/mcmc.txt",
                                  sep = "" ), header = T, sep = "\t" )

# 2. Summarise parameters for all runs 
divtimes3  <- run3[,-c(1, (dim(run3)[2]-8):(dim(run3)[2]) )]
divtimes4  <- run4[,-c(1, (dim(run4)[2]-8):(dim(run4)[2]) )]
divtimes7  <- run7[,-c(1, (dim(run7)[2]-8):(dim(run7)[2]) )]
divtimes8  <- run8[,-c(1, (dim(run8)[2]-8):(dim(run8)[2]) )]
divtimes9  <- run9[,-c(1, (dim(run9)[2]-8):(dim(run9)[2]) )]
divtimes11 <- run11[,-c(1, (dim(run11)[2]-8):(dim(run11)[2]) )]
divtimes13 <- run13[,-c(1, (dim(run13)[2]-8):(dim(run13)[2]) )]
divtimes16 <- run16[,-c(1, (dim(run16)[2]-8):(dim(run16)[2]) )]
divtimes   <- rbind( divtimes3, divtimes4, divtimes7, divtimes8,
                     divtimes9, divtimes11, divtimes13, divtimes16 )

mean_est_divt  <- apply( X = divtimes, MARGIN = 2, FUN = mean )
quant_est_divt <- apply( X = divtimes, MARGIN = 2, FUN = quantile, probs = c(0.025,0.975) )
quant_est_divt <- t( quant_est_divt )
all.equal( quant_est_divt[1,], quantile( divtimes[,1], probs = c( 0.025, 0.975 ) ) )

# Save object with post.divtimes to later be used 
save( divtimes, file = "00_fitST_new/Rdata/post.divtimes.RData" )

#---------------------------#
# FUNCTION TO FIND ST-FITS  #
#---------------------------#
# 1. Create list to store fitted ST-distributions
ST.fitted.dists          <- vector( mode = "list", length = length( colnames( divtimes ) ) )
names( ST.fitted.dists ) <- colnames( divtimes )
ST.fitted.objects          <- vector( mode = "list", length = length( colnames( divtimes ) ) )
names( ST.fitted.objects ) <- colnames( divtimes )

# 2. Set seed and start loop 
set.seed( 12345 )
for ( i in 1:length( colnames( divtimes ) ) ){
  
  # 1. Fit a ST distribution to each node. Then save in 
  #    the lists previously created both the object
  #    output by sn::st.mple and only the "dp" pars
  cat( "Working with node", colnames( divtimes )[i], "...\n\n" )
  write( paste( "Working with node ", colnames( divtimes )[i], "...\n\n", sep = "" ),
         file = "00_fitST_new/logs/log_file_convergence_BFGS.txt", sep = "\n", append = TRUE )
  tmp_node <- sn::st.mple( y = divtimes[,i], opt.method = "BFGS" )
  
  # 2. Check for convergence, otherwise keep trying
  count_tries_conv <- 1
  while( tmp_node$opt.method$convergence != 0 ){
    
    count_tries_conv <- count_tries_conv + 1
    cat( "Convergence has not been reached with node", colnames( divtimes )[i],
         "...\nSEARCH NUMBER", count_tries_conv, "...\n",
         "The parameters found in the previous search:\n",
         tmp_node$dp[1], "|", tmp_node$dp[2], "|",tmp_node$dp[3], "|",tmp_node$dp[4],
         "\nare now used\as starting values now\n\n" )
    write( paste( "Convergence has not been reached with node ", colnames( divtimes )[i],
                  "...\nSEARCH NUMBER ", count_tries_conv, "...\n",
                  "The parameters found in the previous search are used\n",
                  "as starting values now:\n",
                  tmp_node$dp[1], "|", tmp_node$dp[2], "|",tmp_node$dp[3], "|",tmp_node$dp[4],
                  "\n\n", sep = "" ),
           file = "00_fitST_new/logs/log_file_convergence_BFGS.txt", sep = "\n", append = TRUE )
    
    tmp_node <- sn::st.mple( y = divtimes[,i], opt.method = "BFGS",
                             dp = tmp_node$dp.complete )
    
    if( tmp_node$opt.method$convergence == 0 ){
      cat( "Convergenced reached now!\n" )
      cat( "Final parameters for node", colnames( divtimes )[i], "are:\n",
           tmp_node$dp[1], "|", tmp_node$dp[2], "|",tmp_node$dp[3], "|",tmp_node$dp[4],"\n\n" )
      write( paste( "Convergenced reached now!\n\n", 
                    "Final parameters for node ", colnames( divtimes )[i], "are:\n",
                    tmp_node$dp[1], "|", tmp_node$dp[2], "|",tmp_node$dp[3], "|",
                    tmp_node$dp[4],"\n\n", sep = "" ),
             file = "00_fitST_new/logs/log_file_convergence_BFGS.txt", sep = "\n", append = TRUE )
      
    }
    
    if( count_tries_conv == 50 ){
      
      cat( "You have tried 50 times\n",
           "We will try the optimizing approach within this function...\n" )
      write( paste( "You have tried 50 times\n",
                    "We will try the optimizing approach within this function...\n", sep = "" ),
             file = "00_fitST_new/logs/log_file_convergence_BFGS.txt", sep = "\n", append = TRUE )
      count_tries_conv_FUN <- 0
      
      while( tmp_node$opt.method$convergence != 0 ){
        
        count_tries_conv_FUN <- count_tries_conv_FUN + 1
        tmp_node <- sn::st.mple( y = divtimes[,i], dp = tmp_node$dp.complete )
        if( tmp_node$opt.method$convergence == 0 ){
          cat( "Convergenced reached with their method now!\n" )
          cat( "Final parameters for node", colnames( divtimes )[i], "are :\n",
               tmp_node$dp[1], "|", tmp_node$dp[2], "|",tmp_node$dp[3], "|",
               tmp_node$dp[4],"\n\n" )
          write( paste( "Convergenced reached with their method now!\n\n", 
                        "Final parameters for node ", colnames( divtimes )[i], "are:\n",
                        tmp_node$dp[1], "|", tmp_node$dp[2], "|",tmp_node$dp[3], "|",
                        tmp_node$dp[4],"\n\n", sep = "" ),
                 file = "00_fitST_new/logs/log_file_convergence_BFGS.txt", sep = "\n", append = TRUE )
        }
        if( count_tries_conv_FUN == 50 ){
          cat( "You have tried 50 times with their approach,",
               "this is going to be killed\n" )
          write( paste( "You have tried 50 times with their approach,",
                        "this is going to be killed\n", sep = "" ),
                 file = "00_fitST_new/logs/log_file_convergence_BFGS.txt", sep = "\n", append = TRUE )
          break
        }
        
      }
      
    }
    
  }
  
  # Get data 
  ST.fitted.objects[[ i ]] <- tmp_node
  ST.fitted.dists[[ i ]]   <- tmp_node$dp
  
  # 3. Plot fitted ST for each node previously 
  #    computed using the values sampled during the MCMC
  png( filename = paste( "00_fitST_new/plots/Fit_ST_", colnames( divtimes )[i], ".png", sep = "" ),
       width = 1024, height = 768 )
  # 3.1. Find limit axis
  max_x_st    <- round( max( density( rst( n = 1000, dp = tmp_node$dp.complete ) )$x ) + 0.5 )
  max_x_chain <- round( max( density( divtimes[,i] )$x ) + 0.5 )
  x_lim       <- max( max_x_chain, max_x_st )
  max_y_st    <- round( max( density( rst( n = 1000, dp = tmp_node$dp.complete ) )$y ) )
  max_y_chain <- round( max( density( divtimes[,i] )$y ) + 0.5 )
  y_lim       <- max( max_y_chain, max_y_st )
  write( paste( colnames( divtimes )[i], max_x_chain, max_y_chain, sep = "\t" ),
         file = "00_fitST_new/logs/log_limaxis.txt", sep = "\n", append = TRUE )
  # 3.2. Plot
  plot( density( divtimes[,i], adj = 1 ),
        xlim = c( 0, x_lim ), ylim = c( 0, y_lim ), 
        main = paste( colnames( divtimes )[i], " = ",
                      "ST(", paste0( round(tmp_node$dp, 2),
                                     collapse = "," ),
                      ")", sep = "" ) )
  curve( dst( x, xi = tmp_node$dp[1], omega = tmp_node$dp[2],
              alpha = tmp_node$dp[3], nu = tmp_node$dp[4] ),
         from = 0, to = x_lim,
         n = 1e4, add = TRUE, col = "red" )
  dev.off()
  
}

# 3. Save objects so they can later used for extra plots
save( ST.fitted.objects, file = "00_fitST_new/Rdata/ST.fitted.objects.RData" )
save( ST.fitted.dists, file = "00_fitST_new/Rdata/ST.fitted.dists.RData" )

# 4. Transform list into a matrix
mat_ST <- matrix( 0, nrow = length( ST.fitted.dists ), ncol = 6 )
#names( ST.fitted.dists[[1]] )
colnames( mat_ST ) <- c( "xi-location", "omega-scale", "alpha-shape", "nu-df",
                         "MCMCtree-calib", "MCMCtree-calib-rounded" )
rownames( mat_ST ) <- names( ST.fitted.dists )
for ( i in 1:length( ST.fitted.dists ) ){
  
  mat_ST[i,1] <- ST.fitted.dists[[ i ]][1]
  mat_ST[i,2] <- ST.fitted.dists[[ i ]][2]
  mat_ST[i,3] <- ST.fitted.dists[[ i ]][3]
  mat_ST[i,4] <- ST.fitted.dists[[ i ]][4]
  
  tmp.ST.fitted.dists.rounded <- round( ST.fitted.dists[[ i ]], 3 )
  mat_ST[i,5] <- paste( "ST(", ST.fitted.dists[[ i ]][1], ",", ST.fitted.dists[[ i ]][2],
                        ",", ST.fitted.dists[[ i ]][3], ",", ST.fitted.dists[[ i ]][4], ")",
                        sep = "" )
  mat_ST[i,6] <- paste( "ST(", tmp.ST.fitted.dists.rounded[1], ",",
                        tmp.ST.fitted.dists.rounded[2], ",", tmp.ST.fitted.dists.rounded[3],
                        ",", tmp.ST.fitted.dists.rounded[4], ")", sep = "" )
  
}

write.table( mat_ST, file = "00_fitST_new/Rout/ST.fitted.dists.G2.40.tsv", sep = "\t",
             quote = F )

##===============================================##
## STOP THIS SCRIPT HERE. NOW GO TO SCRIPT       ##
## "01_Add_STcalibs_to_tree.R" SO THE TREES      ##
## CAN BE IMPLEMENTED WITH THE ST DISTRIBUTIONS  ##
## THAT HAVE BEEN FITTED TO EACH NODE            ##
##===============================================##
