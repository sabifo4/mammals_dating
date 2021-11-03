#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------#
# SET ENVIRONMENT #
#-----------------#
library( rstudioapi ) 
# Get the path to current open R script and find main dir "00_Gene_filtering"
path_to_file <- getActiveDocumentContext()$path
wd <- paste( dirname( path_to_file ), "/", sep = "" )
wd2 <- gsub( pattern = "01_ESS_and_chain_convergence/", replacement = "", x = wd )
setwd( wd )

#-------------#
# DEFINE FUNS #
#-------------#
# Function to convert the matrix with the sampled values
# for each parameter in the MCMC into a 3D array readable 
# by `rstan::monitor`.
# 
# Parameters:
# x  Data frame, there are many rows as iterations and as many 
#    columns as parameters that had to be inferred during the 
#    MCMC
#
mat2arr <- function( x ){
  
  # Get rid of "Gen" column:
  ## Not needed anymore
  #x <- x[,-c(1)]
  
  # Get array format
  arr <- array( 0, dim = c( dim( x )[1], 1, dim( x )[2] ),
                dimnames = list( paste( "Iter", 1:dim( x )[1], sep = "" ),
                                 "1 chain",
                                 #paste( "Param", 1:dim( x )[2], sep = "" ) )
                                 colnames( x ) )
  )
  
  for( p in 1:dim( x )[2] ){
    arr[,1,p] <- x[,p]
  }
  
  # Return object
  return( arr )
  
}

# Function to run `rstan::monitor` and export the 
# median, min, and max values for the tail-ESS and the 
# bulk-ESS
#
# Parameters:
# x        Data frame, there are many rows as iterations and as many 
#          columns as parameters that had to be inferred during the 
#          MCMC
#
# coda_fun Boolean, TRUE if you want to also compute the ESS
#          with the coda::effectiveSize function
sum_MCMC_ESS <- function( x, coda_fun = FALSE ){
  
  # Define empty vector to report stats
  if( coda_fun == FALSE ){
    out_stats <- matrix( 0, nrow = 3, ncol = 2 )
    colnames( out_stats ) <- c( "Tail-ESS", "Bulk-ESS" )
  } else{ 
    out_stats <- matrix( 0, nrow = 3, ncol = 3 )
    colnames( out_stats ) <- c( "Tail-ESS", "Bulk-ESS", "coda-ESS" )
  }
  rownames( out_stats ) <- c( "Med.", "Min.", "Max." )
  
  # Compute stats
  stats_mcmc     <- rstan::monitor( sims = mat2arr( x ) )
  out_stats[1,1] <- median( stats_mcmc$Tail_ESS )
  out_stats[2,1] <- min( stats_mcmc$Tail_ESS )
  out_stats[3,1] <- max( stats_mcmc$Tail_ESS )
  out_stats[1,2] <- median( stats_mcmc$Bulk_ESS )
  out_stats[2,2] <- min( stats_mcmc$Bulk_ESS )
  out_stats[3,2] <- max( stats_mcmc$Bulk_ESS )
  
  if( coda_fun == TRUE ){
    ESS_coda       <- coda::effectiveSize( x = x )
    out_stats[1,3] <- median( ESS_coda )
    out_stats[2,3] <- min( ESS_coda )
    out_stats[3,3] <- max( ESS_coda )
  }
  
  # Return stats 
  return( list( tab = out_stats, stats = stats_mcmc, stats_CODA = ESS_coda ) )
  
}


# Function to load the two half of the subtrees.
# Used to get mean estimates for convergence plots
# 
# Parameters:
# mcmc    Character, path to mcmc.txt file
# subt    Character, name of the subtree
# delcol  Numeric, number of columns that are to be deleted 
#         as they do not contain divtimes. Default = 10.
load_subt <- function( mcmc1, subt, delcol = 10, perc = 0.975 ){
  
  # 1. Load files and get parameters
  cat( "Load combined mcmc.txt from ", subt, " ... ...\n" )
  run1    <- read.table( mcmc1, header = TRUE, sep = "\t" )
  
  # 2. Summarise parameters
  cat( "Generate objects with summarised estimated divergence times... ...\n")
  dim.r1   <- dim(run1)[2]
  # divtimes <- run1[,-c(1, dim(run1)[2])]
  divtimes <- run1[,-c( 1, ( dim.r1-delcol ):dim.r1 )]
  
  mean_est_divt  <- apply( X = divtimes, MARGIN = 2, FUN = mean )
  quant_est_divt <- apply( X = divtimes, MARGIN = 2, FUN = quantile, probs = c( 1-perc,perc ) )
  quant_est_divt <- t( quant_est_divt )
  test_alleq     <- all.equal( quant_est_divt[1,], quantile( divtimes[,1], probs = c( 1-perc,perc ) ) )
  if( test_alleq != TRUE ){
    stop( "There was an issue calculating quantiles!" )
  }
  
  # 3. Return object 
  cat( "\nTasks done! Return objects\n\n")
  return( list( divt = divtimes, mean_divt = mean_est_divt, quant_divt = quant_est_divt ) )
  
}


# Function to find problematic runs.
# It prints out the mean div time, qlow, and qup 
# for each node and for each run
# 
# Parameters:
# num_dirs       Numeric. Default is c(3,4,7,8,9,11,13,16). Change if needed.
# delcol         Numeric. Default is 8 so only the samples for time estimates are taken.
# name_dir_subt  Character. Name of the directory where the anlayses for this subtree ran.
# num_divt       Numeric. Number of columns of the mcmc.txt that correspond to the samples 
#                collected for the times.
# node_calib     Character. CSV file with two columns, where the first column has the names 
#                of the calibrations and in the second the corresponding node. If not available,
#                set to NULL.
# path_72sp      Character. Path to the directory that contains the 10 runs.
# tree_hyp       Character. Name of the tree hypothesis used to name the directory.
# clock          Character. GBM or ILN, default is "GBM".
# perc           Numeric. Default 0.975, change if you want different quantiles.
# out            Character. Path to where you want the output files saved. Default is 
#                current directory, change if needed.
# maintt         Boolean. TRUE if the main tree is being analysed, which is in a different 
#                directory. FALSE otherwise for the rest of tree hypotheses.
find_prob_MCMC_72sp <- function ( num_dirs = c(3,4,7,8,9,11,13,16), delcol = 8, name_dir_subt, num_divt, node_calib,
                                  path_72sp, tree_hyp, clock = "GBM", out = NULL, perc = 0.975, maintt = TRUE ){
  
  # 0. Allow for numbers not using exponentials if too low !
  options( scipen = 999 )
  
  # 1. Create global vars
  total_runs         <- length( num_dirs )
  subtree_list       <- vector( mode = "list", length = total_runs )
  subtree_meandivt   <- matrix( 0, nrow = total_runs, ncol = num_divt )
  subtree_qup        <- matrix( 0, nrow = total_runs, ncol = num_divt )
  subtree_qlow       <- matrix( 0, nrow = total_runs, ncol = num_divt )
  names( subtree_list )        <- rownames( subtree_meandivt ) <-
    rownames( subtree_qup ) <- rownames( subtree_qlow ) <- paste( "run", num_dirs, sep = "" )
  
  # 2. Get summarised data
  count     <- 0
  count_dim <- vector( "numeric", total_runs )
  for( i in num_dirs ){
      count <- count + 1
      cat( "Loading run ", i, "\n" )
      
      if ( maintt == TRUE ){
        subtree_list[[ count ]]  <- load_subt( mcmc1 = paste( path_72sp, "run", i, "/mcmctree_", clock, "/mcmc.txt", sep = "" ),
                                               subt = paste( name_dir_subt, "_run", i, sep = "" ),
                                               delcol = delcol, perc = perc )
      }else if( maintt == FALSE ){
        subtree_list[[ count ]]  <- load_subt( mcmc1 = paste( path_72sp, "mcmc", i, "/mcmctree_", clock, "/", tree_hyp, "/mcmc.txt", sep = "" ),
                                               subt = paste( name_dir_subt, "_run", i, sep = "" ),
                                               delcol = delcol, perc = perc )
      }
      
      count_dim[count]             <- dim( subtree_list[[ count ]]$divt )[1]
      subtree_meandivt[count,] <- subtree_list[[ count ]]$mean_divt
      subtree_qlow[count,]     <- subtree_list[[ count ]]$quant_divt[,1]
      subtree_qup[count,]      <- subtree_list[[ count ]]$quant_divt[,2]
  }
  colnames( subtree_meandivt ) <- colnames( subtree_qup ) <-
    colnames( subtree_qlow ) <- rownames( subtree_list[[ 1 ]]$quant_divt )
  
  # 3. Get mean data across runs. Apparently it is more accurate 
  #    than actually using meandivt per run and it matches metrics 
  #    calculated by MCMCtree :|
  mcmc_all <- matrix( 0, nrow = sum( count_dim ), ncol = num_divt )
  colnames( mcmc_all ) <- rownames( subtree_list[[ 1 ]]$quant_divt )
  start <- 1
  stop  <- 0
  cat( "Calculating mean and quantiles for all samples... ...\n\n")
  for( i in 1:length(num_dirs) ){
    # cat( "Start: ", start, "\n")
    stop <- stop + count_dim[i]
    # cat( "Stop: ", stop, "\n" )
    mcmc_all[c(start:stop),] <- matrix( unlist( subtree_list[[ i ]]$divt ), nrow = num_divt, byrow = FALSE ) 
    start <- stop + 1
  }
  
  mean_divt    <- apply( X = mcmc_all, MARGIN = 2, FUN = mean )
  mean_quants  <- apply( X = mcmc_all, MARGIN = 2, FUN = quantile, probs = c( 1-perc,perc ) )
  mean_est         <- matrix( 0, nrow = length( mean_divt ), ncol = 3 )
  mean_est_priors  <- matrix( 0, nrow = length( mean_divt ), ncol = 4 )
  mean_est[,1]     <- mean_est_priors[,1] <- round( mean_divt*100, digits = 3 )
  mean_est[,2]     <- mean_est_priors[,2] <- round( as.numeric( format( mean_quants[1,], scientific = FALSE ) )*100, digits = 2 )
  mean_est[,3]     <- mean_est_priors[,3] <- round( as.numeric( format( mean_quants[2,], scientific = FALSE ) )*100, digits = 2 )
  colnames( mean_est )        <- c( "Mean_time", "Mean_qlow", "Mean_qup" )
  colnames( mean_est_priors ) <- c( "Mean_time", "Mean_qlow", "Mean_qup", "Priors" )
  test_names <- all.equal( names( mean_divt ), colnames( mean_quants ) )
  if( test_names != TRUE ){
    stop( "Issue with names for mean divt, mean qup, and mean q low!" )
  }
  rownames( mean_est ) <- rownames( mean_est_priors ) <- names( mean_divt )
  
  # 4. Get matching nodes with calib names
  if( length( node_calib ) > 0 ) {
    match_csv <- read.table( node_calib, header = TRUE, sep = ";", stringsAsFactors = FALSE )
    ind_match <- which( as.numeric( gsub( pattern = "t_n", replacement = "", x = rownames( mean_est_priors ) ) )
                        %in% match_csv[,2] )
    for( i in ind_match ){
      tmp_ind <- which( paste( "t_n", match_csv[,2], sep = "" ) == rownames( mean_est_priors )[i] )
      rownames( mean_est_priors )[i] <- paste( "t_n", match_csv[tmp_ind,2], "_", match_csv[tmp_ind,1], sep = "" )
      mean_est_priors[i,4] <- match_csv[tmp_ind,3]
    }
  }
  
  # 5. Write separate output files for mean times, mean qup, and mean qlow for each 
  #    run
  write.table( x = subtree_meandivt, file = paste( out, "mean_divt_72sp.tsv", sep = "" ), sep = "\t", quote = FALSE )
  write.table( x = subtree_qup, file = paste( out, "mean_qup_7sp.tsv", sep = "" ), sep = "\t", quote = FALSE )
  write.table( x = subtree_qlow, file = paste( out, "mean_qlow_72sp.tsv", sep = "" ), sep = "\t", quote = FALSE )
  write.table( x = mean_est_priors, file = paste( out, "all_mean_est_72sp.tsv", sep = "" ), sep = "\t", quote = FALSE )
  
  cat( "Output files available! Check returned list with all objects generated too :) \n\n")
  return( list( tt_all = subtree_list, mean = subtree_meandivt, qup = subtree_qup, qdown = subtree_qlow,
                all_mean_est = mean_est, all_mcmc = mcmc_all ) )
  
}

# Plot convergence plot 
# Function to plot a convergence plot 
# Arguments:
#
# name_dir_subt  Character. Name of the directory where the anlayses for this subtree ran.
# mean_divt1     Character. Name of the object that contains the min divtimes for the
#                first half of runs.
# mean_divt2     Character. Name of the object that contains the min divtimes for the
#                second half of runs.
# num_runs       Integer. Number of runs.
plot_convergence <- function ( name_dir_subt, mean_divt1, mean_divt2, num_runs ){
  half <- round( num_runs/2, digits = 0 )
  tmp <- bquote( paste(  hat( italic(t) ), " | MCMC - run 1-", .(half), sep = "" ) )
  tmp2 <- bquote( paste( hat( italic(t) ), " | MCMC - run ", .(half+1), "-", .(num_runs), sep = "" ) )
  plot( x = mean_divt1, y = mean_divt2,
        main = paste( "Convergence plot - ", name_dir_subt, sep = "" ),
        xlab = tmp, ylab = tmp2 )
  abline( lm( mean_divt2~0 + mean_divt1 ),
          col="red", lty = 2 )
}


#-----------#
# LOAD DATA #
#-----------#
#-- MAIN TREE (T2) --#
# 72sp tree - NEW
if( ! dir.exists( paste( wd, "out_data/00_post_72sp_new", sep = "" ) ) ){
  dir.create( paste( wd, "out_data/00_post_72sp_new", sep = "" ) )
}
softbounds_72sp_sum  <- find_prob_MCMC_72sp( num_dirs = c(3,4,7,8,9,11,13,16), delcol = 8, name_dir_subt = "72sp",
                                             num_divt = 71, node_calib = "Calibs_nodes_72sp.csv", 
                                             tree_hyp = "02_atlantogenata_tarver2016", maintt = TRUE,
                                             clock = "GBM", out = "out_data/00_post_72sp_new/",
                                             path_72sp = paste( wd2, "00_MCMCtree_analyses/00_main_tree_T2/01_MCMCtree_posterior/02_atlantogenata_tarver2016/",
                                                                sep = "" ),
                                             perc = 0.975 )
#> CHECK: Issues with quantiles
for ( j in 2:8 ){
  tmp.qup     <- softbounds_72sp_sum$qup[1,] - softbounds_72sp_sum$qup[j,]
  tmp.ind.qup <- which( abs( tmp.qup ) > 0.1 )
  if( length( tmp.ind.qup ) > 0 ){
    cat( "q97.5%: Check the following nodes in chain ", j, ":", names( tmp.qup[tmp.ind.qup] ), "\n" )
    cat(  "Difference: ", tmp.qup[tmp.ind.qup], "\n")
  }
  tmp.qdown    <- softbounds_72sp_sum$qdown[1,] - softbounds_72sp_sum$qdown[j,]
  tmp.ind.qdown<- which( abs( tmp.qdown ) > 0.1 )
  if( length( tmp.ind.qdown ) > 0 ){
    cat( "q2.5%: Check the following nodes in chain ", j, ":", names( tmp.qdown[tmp.ind.qdown] ), "\n" )
    cat(  "Difference: ", tmp.qdown[tmp.ind.qdown], "\n")
  }
}
#>END CHECK -- No issues!
sp72_filt_half1 <- apply( X = softbounds_72sp_sum$mean[1:4,], MARGIN = 2, FUN = mean )
sp72_filt_half2 <- apply( X = softbounds_72sp_sum$mean[5:8,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_72sp_NEW.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "72sp", mean_divt1 = sp72_filt_half1,
                  mean_divt2 = sp72_filt_half2, num_runs = 8 )
dev.off()

# 72sp tree - OLD
if( ! dir.exists( paste( wd, "out_data/00_post_72sp_old", sep = "" ) ) ){
  dir.create( paste( wd, "out_data/00_post_72sp_old", sep = "" ) )
}
softbounds_72sp_OLD_sum  <- find_prob_MCMC_72sp( num_dirs = c(1,2,3,5,6,9,10), delcol = 8, name_dir_subt = "72sp",
                                                 num_divt = 71, node_calib = "Calibs_nodes_72sp.csv", 
                                                 tree_hyp = "02_atlantogenata_tarver2016", maintt = FALSE,
                                                 clock = "GBM", out = "out_data/00_post_72sp_old/",
                                                 path_72sp = paste( wd2, "00_MCMCtree_analyses/00_main_tree_T2/01_MCMCtree_posterior_old/", sep = "" ),
                                                 perc = 0.975 )
#> CHECK: Issues with quantiles
for ( j in 2:7 ){
  tmp.qup     <- softbounds_72sp_OLD_sum$qup[1,] - softbounds_72sp_OLD_sum$qup[j,]
  tmp.ind.qup <- which( abs( tmp.qup ) > 0.1 )
  if( length( tmp.ind.qup ) > 0 ){
    cat( "q97.5%: Check the following nodes in chain ", j, ":", names( tmp.qup[tmp.ind.qup] ), "\n" )
    cat(  "Difference: ", tmp.qup[tmp.ind.qup], "\n")
  }
  tmp.qdown    <- softbounds_72sp_OLD_sum$qdown[1,] - softbounds_72sp_OLD_sum$qdown[j,]
  tmp.ind.qdown<- which( abs( tmp.qdown ) > 0.1 )
  if( length( tmp.ind.qdown ) > 0 ){
    cat( "q2.5%: Check the following nodes in chain ", j, ":", names( tmp.qdown[tmp.ind.qdown] ), "\n" )
    cat(  "Difference: ", tmp.qdown[tmp.ind.qdown], "\n")
  }
}
#>END CHECK -- No issues!
sp72_OLD_filt_half1 <- apply( X = softbounds_72sp_OLD_sum$mean[1:4,], MARGIN = 2, FUN = mean )
sp72_OLD_filt_half2 <- apply( X = softbounds_72sp_OLD_sum$mean[5:7,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_72sp_OLD.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "72sp", mean_divt1 = sp72_OLD_filt_half1,
                  mean_divt2 = sp72_OLD_filt_half2, num_runs = 7 )
dev.off()

# 72sp tree - Updated geochronology (Sep 2021)
if( ! dir.exists( paste( wd, "out_data/00_post_72sp_updated_geochron", sep = "" ) ) ){
  dir.create( paste( wd, "out_data/00_post_72sp_updated_geochron", sep = "" ) )
}
softbounds_72sp_updgeochron_sum  <- find_prob_MCMC_72sp( num_dirs = c(1:4), delcol = 8, name_dir_subt = "72sp",
                                                         num_divt = 71, node_calib = "Calibs_nodes_72sp_UPGC.csv", 
                                                         tree_hyp = "02_atlantogenata_tarver2016", maintt = TRUE,
                                                         clock = "GBM", out = "out_data/00_post_72sp_updated_geochron/",
                                                         path_72sp = paste( wd2, "00_MCMCtree_analyses/00_main_tree_T2/02_MCMCtree_posterior_newchrono/",
                                                                    sep = "" ) )
#> CHECK: Issues with quantiles
for ( j in 2:4 ){
    tmp.qup     <- softbounds_72sp_updgeochron_sum$qup[1,] - softbounds_72sp_updgeochron_sum$qup[j,]
    tmp.ind.qup <- which( abs( tmp.qup ) > 0.1 )
    if( length( tmp.ind.qup ) > 0 ){
      cat( "q97.5%: Check the following nodes in chain ", j, ":", names( tmp.qup[tmp.ind.qup] ), "\n" )
      cat(  "Difference: ", tmp.qup[tmp.ind.qup], "\n")
    }
    tmp.qdown    <- softbounds_72sp_updgeochron_sum$qdown[1,] - softbounds_72sp_updgeochron_sum$qdown[j,]
    tmp.ind.qdown<- which( abs( tmp.qdown ) > 0.1 )
    if( length( tmp.ind.qdown ) > 0 ){
      cat( "q2.5%: Check the following nodes in chain ", j, ":", names( tmp.qdown[tmp.ind.qdown] ), "\n" )
      cat(  "Difference: ", tmp.qdown[tmp.ind.qdown], "\n")
    }
}
#>END CHECK -- Issues with chain 2, so we do not use it!
sp72_UPDGC_filt_half1 <- apply( X = softbounds_72sp_updgeochron_sum$mean[c(1,3),], MARGIN = 2, FUN = mean )
sp72_UPDGC_filt_half2 <- softbounds_72sp_updgeochron_sum$mean[4,]
pdf( "plots/Convergence_plot_72sp_UPDGC.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "72sp", mean_divt1 = sp72_UPDGC_filt_half1,
                  mean_divt2 = sp72_UPDGC_filt_half2, num_runs = 3 )
dev.off()


#-- OTHER 6 TREE HYPOTHESES FOR 72sp --#
## T1
T1_72sp_sum  <- find_prob_MCMC_72sp( num_dirs = 1:10, delcol = 8, name_dir_subt = "72sp_tree1",
                                     num_divt = 71, node_calib = NULL, 
                                     path_72sp = paste( wd2, "00_MCMCtree_analyses/01_alternative_tree_hypotheses/00_MCMCtree/", sep = "" ),
                                     maintt = FALSE,
                                     tree_hyp = "01_atlantogenata_scandentia_primates_tarver2016",
                                     clock = "GBM", out = "out_data/01_post_alt_treehyp/tree_1/tree1_", perc = 0.975 )
#> CHECK: Issues with quantiles
for ( j in 2:10 ){
  tmp.qup     <- T1_72sp_sum$qup[1,] - T1_72sp_sum$qup[j,]
  tmp.ind.qup <- which( abs( tmp.qup ) > 0.1 )
  if( length( tmp.ind.qup ) > 0 ){
    cat( "q97.5%: Check the following nodes in chain ", j, ":", names( tmp.qup[tmp.ind.qup] ), "\n" )
    cat(  "Difference: ", tmp.qup[tmp.ind.qup], "\n")
  }
  tmp.qdown    <- T1_72sp_sum$qdown[1,] - T1_72sp_sum$qdown[j,]
  tmp.ind.qdown<- which( abs( tmp.qdown ) > 0.1 )
  if( length( tmp.ind.qdown ) > 0 ){
    cat( "q2.5%: Check the following nodes in chain ", j, ":", names( tmp.qdown[tmp.ind.qdown] ), "\n" )
    cat(  "Difference: ", tmp.qdown[tmp.ind.qdown], "\n")
  }
}
#>END CHECK -- Issues with chain 2
T1_filt_half1 <- apply( X = T1_72sp_sum$mean[1:5,], MARGIN = 2, FUN = mean )
T1_filt_half2 <- apply( X = T1_72sp_sum$mean[6:10,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_T1.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "T1-72sp", mean_divt1 = T1_filt_half1,
                  mean_divt2 = T1_filt_half2, num_runs = 10 )
dev.off()
T1_filt_half1_filt <- apply( X = T1_72sp_sum$mean[c(1,3:6),], MARGIN = 2, FUN = mean )
T1_filt_half2_filt <- apply( X = T1_72sp_sum$mean[7:10,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_T1_filtered.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "T1-72sp", mean_divt1 = T1_filt_half1_filt,
                  mean_divt2 = T1_filt_half2_filt, num_runs = 9 )
dev.off()
T1_72sp_FILT_sum  <- find_prob_MCMC_72sp( num_dirs = c(1,3:10), delcol = 8, name_dir_subt = "72sp_tree1",
                                          num_divt = 71, node_calib = NULL, 
                                          path_72sp = paste( wd2, "00_MCMCtree_analyses/01_alternative_tree_hypotheses/00_MCMCtree/", sep = "" ),
                                          maintt = FALSE,
                                          tree_hyp = "01_atlantogenata_scandentia_primates_tarver2016",
                                          clock = "GBM", out = "out_data/01_post_alt_treehyp/tree_1/tree1_filt_", perc = 0.975 )

## T3
T3_72sp_sum  <- find_prob_MCMC_72sp( num_dirs = 1:10, delcol = 8, name_dir_subt = "72sp_tree3",
                                     num_divt = 71, node_calib = NULL, 
                                     path_72sp = paste( wd2, "00_MCMCtree_analyses/01_alternative_tree_hypotheses/00_MCMCtree/", sep = "" ),
                                     maintt = FALSE,
                                     tree_hyp = "03_atlantogenata_tarver2016_laurasiatheria_dR2012",
                                     clock = "GBM", out = "out_data/01_post_alt_treehyp/tree_3/tree3_", perc = 0.975 )
#> CHECK: Issues with quantiles
for ( j in 2:10 ){
  tmp.qup     <- T3_72sp_sum$qup[1,] - T3_72sp_sum$qup[j,]
  tmp.ind.qup <- which( abs( tmp.qup ) > 0.1 )
  if( length( tmp.ind.qup ) > 0 ){
    cat( "q97.5%: Check the following nodes in chain ", j, ":", names( tmp.qup[tmp.ind.qup] ), "\n" )
    cat(  "Difference: ", tmp.qup[tmp.ind.qup], "\n")
  }
  tmp.qdown    <- T3_72sp_sum$qdown[1,] - T3_72sp_sum$qdown[j,]
  tmp.ind.qdown<- which( abs( tmp.qdown ) > 0.1 )
  if( length( tmp.ind.qdown ) > 0 ){
    cat( "q2.5%: Check the following nodes in chain ", j, ":", names( tmp.qdown[tmp.ind.qdown] ), "\n" )
    cat(  "Difference: ", tmp.qdown[tmp.ind.qdown], "\n")
  }
}
#>END CHECK -- No issues!
T3_filt_half1 <- apply( X = T3_72sp_sum$mean[1:5,], MARGIN = 2, FUN = mean )
T3_filt_half2 <- apply( X = T3_72sp_sum$mean[6:10,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_T3.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "T3-72sp", mean_divt1 = T3_filt_half1,
                  mean_divt2 = T3_filt_half2, num_runs = 10 )
dev.off()

## T4
T4_72sp_sum  <- find_prob_MCMC_72sp( num_dirs = 1:10, delcol = 8, name_dir_subt = "72sp_tree4",
                                     num_divt = 71, node_calib = NULL, 
                                     path_72sp = paste( wd2, "00_MCMCtree_analyses/01_alternative_tree_hypotheses/00_MCMCtree/", sep = "" ),
                                     maintt = FALSE,
                                     tree_hyp = "04_atlantogenata_tarver2016_laurasiatheriaensembl",
                                     clock = "GBM", out = "out_data/01_post_alt_treehyp/tree_4/tree4_", perc = 0.975 )
#> CHECK: Issues with quantiles
for ( j in 2:10 ){
  tmp.qup     <- T4_72sp_sum$qup[1,] - T4_72sp_sum$qup[j,]
  tmp.ind.qup <- which( abs( tmp.qup ) > 0.1 )
  if( length( tmp.ind.qup ) > 0 ){
    cat( "q97.5%: Check the following nodes in chain ", j, ":", names( tmp.qup[tmp.ind.qup] ), "\n" )
    cat(  "Difference: ", tmp.qup[tmp.ind.qup], "\n")
  }
  tmp.qdown    <- T4_72sp_sum$qdown[1,] - T4_72sp_sum$qdown[j,]
  tmp.ind.qdown<- which( abs( tmp.qdown ) > 0.1 )
  if( length( tmp.ind.qdown ) > 0 ){
    cat( "q2.5%: Check the following nodes in chain ", j, ":", names( tmp.qdown[tmp.ind.qdown] ), "\n" )
    cat(  "Difference: ", tmp.qdown[tmp.ind.qdown], "\n")
  }
}
#>END CHECK -- No issues!
T4_filt_half1 <- apply( X = T4_72sp_sum$mean[1:5,], MARGIN = 2, FUN = mean )
T4_filt_half2 <- apply( X = T4_72sp_sum$mean[6:10,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_T4.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "T4-72sp", mean_divt1 = T4_filt_half1,
                  mean_divt2 = T4_filt_half2, num_runs = 10 )
dev.off()

## T5
T5_72sp_sum  <- find_prob_MCMC_72sp( num_dirs = 1:10, delcol = 8, name_dir_subt = "72sp_tree5",
                                     num_divt = 71, node_calib = NULL, 
                                     path_72sp = paste( wd2, "00_MCMCtree_analyses/01_alternative_tree_hypotheses/00_MCMCtree/", sep = "" ),
                                     maintt = FALSE,
                                     tree_hyp = "05_ensembl",
                                     clock = "GBM", out = "out_data/01_post_alt_treehyp/tree_5/tree5_", perc = 0.975 )
#> CHECK: Issues with quantiles
for ( j in 2:10 ){
  tmp.qup     <- T5_72sp_sum$qup[1,] - T5_72sp_sum$qup[j,]
  tmp.ind.qup <- which( abs( tmp.qup ) > 0.1 )
  if( length( tmp.ind.qup ) > 0 ){
    cat( "q97.5%: Check the following nodes in chain ", j, ":", names( tmp.qup[tmp.ind.qup] ), "\n" )
    cat(  "Difference: ", tmp.qup[tmp.ind.qup], "\n")
  }
  tmp.qdown    <- T5_72sp_sum$qdown[1,] - T5_72sp_sum$qdown[j,]
  tmp.ind.qdown<- which( abs( tmp.qdown ) > 0.1 )
  if( length( tmp.ind.qdown ) > 0 ){
    cat( "q2.5%: Check the following nodes in chain ", j, ":", names( tmp.qdown[tmp.ind.qdown] ), "\n" )
    cat(  "Difference: ", tmp.qdown[tmp.ind.qdown], "\n")
  }
}
#>END CHECK -- Issues with chain 10!
T5_filt_half1 <- apply( X = T5_72sp_sum$mean[1:5,], MARGIN = 2, FUN = mean )
T5_filt_half2 <- apply( X = T5_72sp_sum$mean[6:10,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_T5.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "T5-72sp", mean_divt1 = T5_filt_half1,
                  mean_divt2 = T5_filt_half2, num_runs = 10 )
dev.off()
T5_filt_half1_filt <- apply( X = T5_72sp_sum$mean[1:5,], MARGIN = 2, FUN = mean )
T5_filt_half2_filt <- apply( X = T5_72sp_sum$mean[6:9,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_T5_filtered.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "T5-72sp", mean_divt1 = T5_filt_half1_filt,
                  mean_divt2 = T5_filt_half2_filt, num_runs = 9 )
dev.off()
T5_72sp_FILT_sum  <- find_prob_MCMC_72sp( num_dirs = 1:9, delcol = 8, name_dir_subt = "72sp_tree5",
                                          num_divt = 71, node_calib = NULL, 
                                          path_72sp = paste( wd2, "00_MCMCtree_analyses/01_alternative_tree_hypotheses/00_MCMCtree/", sep = "" ),
                                          maintt = FALSE,
                                          tree_hyp = "05_ensembl",
                                          clock = "GBM", out = "out_data/01_post_alt_treehyp/tree_5/tree5_filt_", perc = 0.975 )

## T6
T6_72sp_sum  <- find_prob_MCMC_72sp( num_dirs = 1:10, delcol = 8, name_dir_subt = "72sp_tree6",
                                     num_divt = 71, node_calib = NULL, 
                                     path_72sp = paste( wd2, "00_MCMCtree_analyses/01_alternative_tree_hypotheses/00_MCMCtree/", sep = "" ),
                                     maintt = FALSE,
                                     tree_hyp = "06_xenarthra_tarver2016",
                                     clock = "GBM", out = "out_data/01_post_alt_treehyp/tree_6/tree6_", perc = 0.975 )
#> CHECK: Issues with quantiles
for ( j in 2:10 ){
  tmp.qup     <- T6_72sp_sum$qup[1,] - T6_72sp_sum$qup[j,]
  tmp.ind.qup <- which( abs( tmp.qup ) > 0.1 )
  if( length( tmp.ind.qup ) > 0 ){
    cat( "q97.5%: Check the following nodes in chain ", j, ":", names( tmp.qup[tmp.ind.qup] ), "\n" )
    cat(  "Difference: ", tmp.qup[tmp.ind.qup], "\n")
  }
  tmp.qdown    <- T6_72sp_sum$qdown[1,] - T6_72sp_sum$qdown[j,]
  tmp.ind.qdown<- which( abs( tmp.qdown ) > 0.1 )
  if( length( tmp.ind.qdown ) > 0 ){
    cat( "q2.5%: Check the following nodes in chain ", j, ":", names( tmp.qdown[tmp.ind.qdown] ), "\n" )
    cat(  "Difference: ", tmp.qdown[tmp.ind.qdown], "\n")
  }
}
#>END CHECK -- No issues!
T6_filt_half1 <- apply( X = T6_72sp_sum$mean[1:5,], MARGIN = 2, FUN = mean )
T6_filt_half2 <- apply( X = T6_72sp_sum$mean[6:10,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_T6.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "T6-72sp", mean_divt1 = T6_filt_half1,
                  mean_divt2 = T6_filt_half2, num_runs = 10 )
dev.off()

## T7
T7_72sp_sum  <- find_prob_MCMC_72sp( num_dirs = 1:10, delcol = 8, name_dir_subt = "72sp_tree7",
                                     num_divt = 71, node_calib = NULL, 
                                     path_72sp = paste( wd2, "00_MCMCtree_analyses/01_alternative_tree_hypotheses/00_MCMCtree/", sep = "" ),
                                     maintt = FALSE,
                                     tree_hyp = "07_afrotheria_tarver2016",
                                     clock = "GBM", out = "out_data/01_post_alt_treehyp/tree_7/tree7_", perc = 0.975 )
#> CHECK: Issues with quantiles
for ( j in c(1,3:10) ){
  tmp.qup     <- T7_72sp_sum$qup[2,] - T7_72sp_sum$qup[j,]
  tmp.ind.qup <- which( abs( tmp.qup ) > 0.1 )
  if( length( tmp.ind.qup ) > 0 ){
    cat( "q97.5%: Check the following nodes in chain ", j, ":", names( tmp.qup[tmp.ind.qup] ), "\n" )
    cat(  "Difference: ", tmp.qup[tmp.ind.qup], "\n")
  }
  tmp.qdown    <- T7_72sp_sum$qdown[2,] - T7_72sp_sum$qdown[j,]
  tmp.ind.qdown<- which( abs( tmp.qdown ) > 0.1 )
  if( length( tmp.ind.qdown ) > 0 ){
    cat( "q2.5%: Check the following nodes in chain ", j, ":", names( tmp.qdown[tmp.ind.qdown] ), "\n" )
    cat(  "Difference: ", tmp.qdown[tmp.ind.qdown], "\n")
  }
}
#>END CHECK -- Issues with chains 1 and 4!
T7_filt_half1 <- apply( X = T7_72sp_sum$mean[1:5,], MARGIN = 2, FUN = mean )
T7_filt_half2 <- apply( X = T7_72sp_sum$mean[6:10,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_T7.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "T7-72sp", mean_divt1 = T7_filt_half1,
                  mean_divt2 = T7_filt_half2, num_runs = 10 )
dev.off()
T7_filt_half1_filt <- apply( X = T7_72sp_sum$mean[c(2:3,5:6),], MARGIN = 2, FUN = mean )
T7_filt_half2_filt <- apply( X = T7_72sp_sum$mean[7:10,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_T7_filtered.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "T7-72sp", mean_divt1 = T7_filt_half1_filt,
                  mean_divt2 = T7_filt_half2_filt, num_runs = 8 )
dev.off()
T7_72sp_FILT_sum  <- find_prob_MCMC_72sp( num_dirs = c(2:3,5:10), delcol = 8, name_dir_subt = "72sp_tree7",
                                     num_divt = 71, node_calib = NULL, 
                                     path_72sp = paste( wd2, "00_MCMCtree_analyses/01_alternative_tree_hypotheses/00_MCMCtree/", sep = "" ),
                                     maintt = FALSE,
                                     tree_hyp = "07_afrotheria_tarver2016",
                                     clock = "GBM", out = "out_data/01_post_alt_treehyp/tree_7/tree7_filt_", perc = 0.975 )

## ALL 72SP HYPOTHESES !
pdf( "plots/Convergence_plot_all_72sp_hypotheses.pdf", paper = "a4" )
par( mfrow = c(4,2) )
plot_convergence( name_dir_subt = "main_72sp", mean_divt1 = sp72_filt_half1,
                  mean_divt2 = sp72_filt_half2, num_runs = 8 )
# plot_convergence( name_dir_subt = "T1_72sp", mean_divt1 = T1_filt_half1,
#                   mean_divt2 = T1_filt_half2, num_runs = 10 )
plot_convergence( name_dir_subt = "T1_72sp", mean_divt1 = T1_filt_half1_filt,
                  mean_divt2 = T1_filt_half2_filt, num_runs = 9 )
plot_convergence( name_dir_subt = "T3_72sp", mean_divt1 = T3_filt_half1,
                  mean_divt2 = T3_filt_half2, num_runs = 10 )
plot_convergence( name_dir_subt = "T4_72sp", mean_divt1 = T4_filt_half1,
                  mean_divt2 = T4_filt_half2, num_runs = 10 )
# plot_convergence( name_dir_subt = "T5_72sp", mean_divt1 = T5_filt_half1,
#                   mean_divt2 = T5_filt_half2, num_runs = 10 )
plot_convergence( name_dir_subt = "T5_72sp", mean_divt1 = T5_filt_half1_filt,
                  mean_divt2 = T5_filt_half2_filt, num_runs = 9 )
plot_convergence( name_dir_subt = "T6_72sp", mean_divt1 = T6_filt_half1,
                  mean_divt2 = T6_filt_half2, num_runs = 10 )
# plot_convergence( name_dir_subt = "T7_72sp", mean_divt1 = T7_filt_half1,
#                   mean_divt2 = T7_filt_half2, num_runs = 10 )
plot_convergence( name_dir_subt = "T7_72sp", mean_divt1 = T7_filt_half1_filt,
                  mean_divt2 = T7_filt_half2_filt, num_runs = 8 )
dev.off()


pdf( "plots/Convergence_plot_all_72sp_hypotheses_UPDGC.pdf", paper = "a4" )
par( mfrow = c(4,2) )
plot_convergence( name_dir_subt = "main_72sp_UPDGC", mean_divt1 = sp72_UPDGC_filt_half1,
                  mean_divt2 = sp72_UPDGC_filt_half2, num_runs = 3 )
# plot_convergence( name_dir_subt = "T1_72sp", mean_divt1 = T1_filt_half1,
#                   mean_divt2 = T1_filt_half2, num_runs = 10 )
plot_convergence( name_dir_subt = "T1_72sp", mean_divt1 = T1_filt_half1_filt,
                  mean_divt2 = T1_filt_half2_filt, num_runs = 9 )
plot_convergence( name_dir_subt = "T3_72sp", mean_divt1 = T3_filt_half1,
                  mean_divt2 = T3_filt_half2, num_runs = 10 )
plot_convergence( name_dir_subt = "T4_72sp", mean_divt1 = T4_filt_half1,
                  mean_divt2 = T4_filt_half2, num_runs = 10 )
# plot_convergence( name_dir_subt = "T5_72sp", mean_divt1 = T5_filt_half1,
#                   mean_divt2 = T5_filt_half2, num_runs = 10 )
plot_convergence( name_dir_subt = "T5_72sp", mean_divt1 = T5_filt_half1_filt,
                  mean_divt2 = T5_filt_half2_filt, num_runs = 9 )
plot_convergence( name_dir_subt = "T6_72sp", mean_divt1 = T6_filt_half1,
                  mean_divt2 = T6_filt_half2, num_runs = 10 )
# plot_convergence( name_dir_subt = "T7_72sp", mean_divt1 = T7_filt_half1,
#                   mean_divt2 = T7_filt_half2, num_runs = 10 )
plot_convergence( name_dir_subt = "T7_72sp", mean_divt1 = T7_filt_half1_filt,
                  mean_divt2 = T7_filt_half2_filt, num_runs = 8 )
dev.off()



#---------------#
# CALCULATE ESS #
#---------------#
#> ESS with RStan
# Each column is assumed to be an MCMC. Rows are iterations for parameter X
# Source explaining why it is preferable than the function in coda:
# https://nature.berkeley.edu/~pdevalpine/MCMC_comparisons/nimble_MCMC_comparisons.html

#\\\\\\\\\#
# MAIN T2 #
#---------#
ESS_s72sp        <- sum_MCMC_ESS( x = softbounds_72sp_sum$all_mcmc, coda_fun = TRUE )
dim(softbounds_72sp_sum$all_mcmc)[1]
# length = 160008
min(ESS_s72sp$stats$Rhat)
# min(Rhat) 1.000024
ESS_s72sp$tab
#      Tail-ESS Bulk-ESS   coda-ESS
# Med.      929      374 1255.2725
# Min.      216       56  179.1359
# Max.    11669     8358 9867.9955

ESS_s72sp_UPDGC  <- sum_MCMC_ESS( x = softbounds_72sp_updgeochron_sum$all_mcmc, coda_fun = TRUE )
dim(softbounds_72sp_updgeochron_sum$all_mcmc)[1]
# length = 80004
min(ESS_s72sp_UPDGC$stats$Rhat)
# min(Rhat) 0.9999882
ESS_s72sp_UPDGC$tab
# Tail-ESS Bulk-ESS   coda-ESS
# Med.      271      210  508.18333
# Min.       91       31   81.36824
# Max.     5377     3801 7693.93093

ESS_s72sp_OLD    <- sum_MCMC_ESS( x = softbounds_72sp_OLD_sum$all_mcmc, coda_fun = TRUE )
# length = 140007
# min(Rhat) 0.999986
# Tail-ESS Bulk-ESS   coda-ESS
# Med.     1727      796  1898.5583
# Min.      223      130   268.2216
# Max.    17543    11285 24970.0482

#\\\\#
# T1 #
#----#
ESS_T1_72sp      <- sum_MCMC_ESS( x = T1_72sp_sum$all_mcmc, coda_fun = TRUE )
dim( T1_72sp_sum$all_mcmc )[1]
# length = 200010
min( ESS_T1_72sp$stats$Rhat )
# min( Rhat ) = 0.9999984
ESS_T1_72sp$tab
# Tail-ESS Bulk-ESS   coda-ESS
# Med.     2243     1149  2410.3653
# Min.      530      160   249.1687
# Max.    25359    18028 34733.4022

ESS_T1_filt_72sp <- sum_MCMC_ESS( x = T1_72sp_FILT_sum$all_mcmc, coda_fun = TRUE )
dim( T1_72sp_FILT_sum$all_mcmc )[1]
# length = 180009
min( ESS_T1_filt_72sp$stats$Rhat )
# min( Rhat ) = 0.9999901
ESS_T1_filt_72sp$tab
# Tail-ESS Bulk-ESS   coda-ESS
# Med.     2074     1062  2494.8053
# Min.      507      142   284.6633
# Max.    22442    16211 31304.7031

#\\\\#
# T3 #
#----#
ESS_T3_72sp      <- sum_MCMC_ESS( x = T3_72sp_sum$all_mcmc, coda_fun = TRUE )
dim( T3_72sp_sum$all_mcmc )[1]
# length = 200010
min( ESS_T3_72sp$stats$Rhat )
# min( Rhat ) = 0.999998
ESS_T3_72sp$tab
# Tail-ESS Bulk-ESS   coda-ESS
# Med.     1839      949  2141.2830
# Min.      503       98   352.6019
# Max.    28525    18054 34458.9801

#\\\\#
# T4 #
#----#
ESS_T4_72sp      <- sum_MCMC_ESS( x = T4_72sp_sum$all_mcmc, coda_fun = TRUE )
dim( T4_72sp_sum$all_mcmc )[1]
# length = 200010
min( ESS_T4_72sp$stats$Rhat )
# min( Rhat ) = 0.9999954
ESS_T4_72sp$tab
# Tail-ESS Bulk-ESS  coda-ESS
# Med.     2007      968  2301.386
# Min.      451      177   373.214
# Max.    27303    19469 39685.681

#\\\\#
# T5 #
#----#
ESS_T5_72sp      <- sum_MCMC_ESS( x = T5_72sp_sum$all_mcmc, coda_fun = TRUE )
dim( T5_72sp_sum$all_mcmc )[1]
# length = 200010
min( ESS_T5_72sp$stats$Rhat )
# min( Rhat ) = 1.000006
ESS_T5_72sp$tab
# Tail-ESS Bulk-ESS  coda-ESS
# Med.     1085      745 1857.6749
# Min.      330      142  367.5558
# Max.    10119     6628 9738.5197

ESS_T5_filt_72sp <- sum_MCMC_ESS( x = T5_72sp_FILT_sum$all_mcmc, coda_fun = TRUE )
dim( T5_72sp_FILT_sum$all_mcmc )[1]
# length = 180009
min( ESS_T5_filt_72sp$stats$Rhat )
# min( Rhat ) = 0.9999916
ESS_T5_filt_72sp$tab
# Tail-ESS Bulk-ESS  coda-ESS
# Med.     1118      709  1765.0908
# Min.      361      125   330.1042
# Max.    15411    11241 21647.0629

#\\\\#
# T6 #
#----#
ESS_T6_72sp      <- sum_MCMC_ESS( x = T6_72sp_sum$all_mcmc, coda_fun = TRUE )
dim( T6_72sp_sum$all_mcmc )[1]
# length = 200010
min( ESS_T6_72sp$stats$Rhat )
# min( Rhat ) = 0.9999919
ESS_T6_72sp$tab
# Tail-ESS Bulk-ESS   coda-ESS
# Med.     2721     1295  2633.6739
# Min.      442      193   354.5854
# Max.    29174    19324 16135.7356

#\\\\#
# T7 #
#----#
ESS_T7_72sp      <- sum_MCMC_ESS( x = T7_72sp_sum$all_mcmc, coda_fun = TRUE )
dim( T7_72sp_sum$all_mcmc )[1]
# length = 200010
min( ESS_T7_72sp$stats$Rhat )
# min( Rhat ) = 1.000016
ESS_T7_72sp$tab
# Tail-ESS Bulk-ESS   coda-ESS
# Med.     1725      962  2505.2485
# Min.      259      132   313.6458
# Max.    27271    18096 10359.6552

ESS_T7_filt_72sp <- sum_MCMC_ESS( x = T7_72sp_FILT_sum$all_mcmc, coda_fun = TRUE )
dim( T7_72sp_FILT_sum$all_mcmc )[1]
# length = 160008
min( ESS_T7_filt_72sp$stats$Rhat )
# min( Rhat ) = 1.000035
ESS_T7_filt_72sp$tab
# Tail-ESS Bulk-ESS   coda-ESS
# Med.     1362      751  2147.1282
# Min.      196      103   265.7601
# Max.    21129    14015 28987.7681
