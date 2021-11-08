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
wd2 <- gsub( pattern = "00_ESS_and_chain_convergence/", replacement = "", x = wd )
wd3 <- gsub( pattern = "02_MCMCtree_updcrh/", replacement = "02_MCMCtree/", x = wd2 )
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
# clean   Boolean, TRUE if the `mcmc.txt` file does not have a header,
#         FALSE otherwise.
load_subt <- function( mcmc1, subt, delcol = 10, perc = 0.975  ){
  
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
# num_dirs 
# num_chains     Numeric. Default is 2 chains per run. Change if needed.
# delcol         Numeric. Default is 10 so only the samples for time estimates are taken.
# name_dir_subt  Character. Name of the directory/Path tho the directory where the anlayses
#                for this subtree ran.
# num_divt       Numeric. Number of columns of the mcmc.txt that correspond to the samples 
#                collected for the times.
# node_calib     Character. CSV file with two columns, where the first column has the names 
#                of the calibrations and in the second the corresponding node.
# perc           Numeric. Percentile to calculate the quantiles. Default: 0.975.
# clean          Boolean. FALSE if a "clean" file removing incomplete lines in the `mcmc.txt`
#                has not been created. TRUE otherwise.
# new            Boolean. FALSE if the first type of file architecture is used. TRUE otherwise.
#                Note that if new = TRUE, directories `mcmc1` and `mcmc2` do not appear, hence 
#                the path changes.
# subtree        Character, name of the subtree.
# out_dat        Character, path to the directory to save output data. Do not add the last "/".
find_prob_MCMC <- function ( num_dirs = 6, num_chains = 2, delcol = 10, name_dir_subt, num_divt,
                             subtree, node_calib, perc = 0.975, clean = FALSE, new = FALSE,
                             new_rev = FALSE, out_dat ){
  
  # 0. Allow for numbers not using exponentials if too low !
  options( scipen = 999 )
  
  # 1. Create global vars
  if( length( num_dirs ) == 1 ){
    run_dirs   <- num_dirs           # If all dirs are there
    seq_dirs   <- c(1:run_dirs)
  }else{
    run_dirs   <- length( num_dirs ) # If some dirs have been deleted
    seq_dirs   <- num_dirs
  }
  total_runs         <- run_dirs*num_chains
  subtree_list       <- vector( mode = "list", length = total_runs )
  subtree_meandivt   <- matrix( 0, nrow = total_runs, ncol = num_divt )
  subtree_qup        <- matrix( 0, nrow = total_runs, ncol = num_divt )
  subtree_qlow       <- matrix( 0, nrow = total_runs, ncol = num_divt )
  names( subtree_list )        <- rownames( subtree_meandivt ) <-
    rownames( subtree_qup ) <- rownames( subtree_qlow ) <- paste( "run", 1:total_runs, sep = "" )
  
  # 2. Get summarised data
  count <- 0
  count_dim <- vector( "numeric", total_runs )
  #for( i in 1:run_dirs ){
  for( i in seq_dirs ){
    for ( j in 1:num_chains ){
      count <- count + 1
      cat( "Loading run ", i, " and chain ", j, "\n" )
      if( new_rev == TRUE ){
        subtree_list[[ count ]]  <- load_subt( mcmc1 = paste( name_dir_subt, "/run", i, "/mcmctree_GBM/mcmc_clean.txt", sep = "" ),
                                               subt = paste( name_dir_subt, " - run", i, sep = "" ),
                                               delcol = delcol, perc = perc )
      }else{
        if( new == TRUE ){
          subtree_list[[ count ]]  <- load_subt( mcmc1 = paste( name_dir_subt, "/7/run", i, "/mcmctree_GBM/mcmc_clean.txt", sep = "" ),
                                                 subt = paste( name_dir_subt, " - run", i, sep = "" ),
                                                 delcol = delcol, perc = perc )
        }else if( new == FALSE ){
          if( clean == FALSE ){
            subtree_list[[ count ]]  <- load_subt( mcmc1 = paste( name_dir_subt, "/7/run", i, "/mcmc", j, "/mcmctree_GBM/mcmc.txt", sep = "" ),
                                                   subt = paste( name_dir_subt, "-", i, "_run", j, sep = "" ),
                                                   delcol = delcol, perc = perc )
          }else if( clean == TRUE ){
            subtree_list[[ count ]]  <- load_subt( mcmc1 = paste( name_dir_subt, "/7/run", i, "/mcmc", j, "/mcmctree_GBM/mcmc_clean.txt", sep = "" ),
                                                   subt = paste( name_dir_subt, "-", i, "_run", j, sep = "" ),
                                                   delcol = delcol, perc = perc )
          }
        }
      }
      
      count_dim[count]         <- dim( subtree_list[[ count ]]$divt )[1]
      subtree_meandivt[count,] <- subtree_list[[ count ]]$mean_divt
      subtree_qlow[count,]     <- subtree_list[[ count ]]$quant_divt[,1]
      subtree_qup[count,]      <- subtree_list[[ count ]]$quant_divt[,2]
    }
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
  for( i in 1:total_runs ){
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
  match_csv <- read.table( node_calib, header = TRUE, sep = ";", stringsAsFactors = FALSE )
  ind_match <- which( as.numeric( gsub( pattern = "t_n", replacement = "", x = rownames( mean_est_priors ) ) )
                      %in% match_csv[,2] )
  for( i in ind_match ){
    tmp_ind <- which( paste( "t_n", match_csv[,2], sep = "" ) == rownames( mean_est_priors )[i] )
    rownames( mean_est_priors )[i] <- paste( "t_n", match_csv[tmp_ind,2], "_", match_csv[tmp_ind,1], sep = "" )
    mean_est_priors[i,4] <- match_csv[tmp_ind,3]
  }
  
  # 5. Write separate output files for mean times, mean qup, and mean qlow for each 
  #    run
  if ( ! dir.exists( paste( out_dat, "/", subtree, sep = "" ) ) ){
    dir.create( paste( out_dat, "/", subtree, sep = "" ) )
  }
  write.table( x = subtree_meandivt, file = paste( out_dat, "/", subtree, "/mean_divt.tsv", sep = "" ),
               sep = "\t", quote = FALSE )
  write.table( x = subtree_qup, file = paste( out_dat, "/", subtree, "/mean_qup.tsv", sep = "" ), sep = "\t",
               quote = FALSE )
  write.table( x = subtree_qlow, file = paste( out_dat, "/", subtree, "/mean_qlow.tsv", sep = "" ), sep = "\t",
               quote = FALSE )
  if( new_rev == TRUE ){
    write.table( x = mean_est_priors, file = paste( out_dat, "/", subtree, "/all_mean_est.tsv", sep = "" ),
                 sep = "\t", quote = FALSE )
  }else{
    write.table( x = mean_est_priors,
                 file = paste( out_dat, "/", subtree, "/", subtree, "_all_mean_est.tsv", sep = "" ),
                 sep = "\t", quote = FALSE )
  }
  
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
# Subtrees

#\\\\\\\\\\\\#
# AFROTHERIA #
#------------#
Afro_sum        <- find_prob_MCMC( num_dirs = 32, num_chains = 1, delcol = 10,
                                   name_dir_subt = paste( wd2, "Afrotheria", sep = "" ),
                                   num_divt = 59, subtree = "Afrotheria",
                                   node_calib = paste( wd2, "Afrotheria/Calibs_nodes_Afrotheria.csv", sep = "" ),
                                   perc = 0.975, new = TRUE, clean = TRUE,
                                   out_dat = paste( wd, "out_data", sep = "" ) )

#> CHECK: Issues with quantiles
for ( j in 2:32 ){
  tmp.qup     <- Afro_sum$qup[1,] - Afro_sum$qup[j,]
  tmp.ind.qup <- which( abs( tmp.qup ) > 0.1 )
  if( length( tmp.ind.qup ) > 0 ){
    cat( "q97.5%: Check the following nodes in chain ", j, ":", names( tmp.qup[tmp.ind.qup] ), "\n" )
    cat(  "Difference: ", tmp.qup[tmp.ind.qup], "\n")
  }
  tmp.qdown    <- Afro_sum$qdown[1,] - Afro_sum$qdown[j,]
  tmp.ind.qdown<- which( abs( tmp.qdown ) > 0.1 )
  if( length( tmp.ind.qdown ) > 0 ){
    cat( "q2.5%: Check the following nodes in chain ", j, ":", names( tmp.qdown[tmp.ind.qdown] ), "\n" )
    cat(  "Difference: ", tmp.qdown[tmp.ind.qdown], "\n")
  }
}
#>END CHECK -- No issues!
Afro_filt_half1 <- apply( X = Afro_sum$mean[1:17,], MARGIN = 2, FUN = mean )
Afro_filt_half2 <- apply( X = Afro_sum$mean[17:32,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_Afrotheria.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "Afrotheria", mean_divt1 = Afro_filt_half1,
                  mean_divt2 = Afro_filt_half2, num_runs = 32 )
dev.off() 

#\\\\\\\\\\\#
# XENARTHRA #
#-----------#
Xen_sum        <- find_prob_MCMC( num_dirs = 32, num_chains = 1, delcol = 10,
                                  name_dir_subt = paste( wd2, "Xenarthra", sep = "" ),
                                  num_divt = 32, subtree = "Xenarthra",
                                  node_calib = paste( wd2, "Xenarthra/Calibs_nodes_Xenarthra.csv", sep = "" ),
                                  perc = 0.975, clean = TRUE, new = TRUE,
                                  out_dat = paste( wd, "out_data", sep = "" ) )

#> CHECK: Issues with quantiles
for ( j in 2:32 ){
  tmp.qup     <- Xen_sum$qup[1,] - Xen_sum$qup[j,]
  tmp.ind.qup <- which( abs( tmp.qup ) > 0.1 )
  if( length( tmp.ind.qup ) > 0 ){
    cat( "q97.5%: Check the following nodes in chain ", j, ":", names( tmp.qup[tmp.ind.qup] ), "\n" )
    cat(  "Difference: ", tmp.qup[tmp.ind.qup], "\n")
  }
  tmp.qdown    <- Xen_sum$qdown[1,] - Xen_sum$qdown[j,]
  tmp.ind.qdown<- which( abs( tmp.qdown ) > 0.1 )
  if( length( tmp.ind.qdown ) > 0 ){
    cat( "q2.5%: Check the following nodes in chain ", j, ":", names( tmp.qdown[tmp.ind.qdown] ), "\n" )
    cat(  "Difference: ", tmp.qdown[tmp.ind.qdown], "\n")
  }
}
#>END CHECK -- No issues!
Xen_filt_half1 <- apply( X = Xen_sum$mean[1:16,], MARGIN = 2, FUN = mean )
Xen_filt_half2 <- apply( X = Xen_sum$mean[17:32,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_Xenarthra.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "Xenarthra", mean_divt1 = Xen_filt_half1,
                  mean_divt2 = Xen_filt_half2, num_runs = 32 )
dev.off()

#\\\\\\\\\\\\#
# LAGOMORPHA #
#------------#
Lag_sum        <- find_prob_MCMC( num_dirs = 32, num_chains = 1, delcol = 10,
                                  name_dir_subt = paste( wd2, "Lagomorpha", sep = "" ),
                                  num_divt = 87, subtree = "Lagomorpha",
                                  node_calib = paste( wd2, "Lagomorpha/Calibs_nodes_Lagomorpha.csv", sep = "" ),
                                  perc = 0.975, clean = TRUE, new = TRUE,
                                  out_dat = paste( wd, "out_data", sep = "" ) )

#> CHECK: Issues with quantiles
for ( j in 1:32 ){
  tmp.qup     <- Lag_sum$qup[1,] - Lag_sum$qup[j,]
  tmp.ind.qup <- which( abs( tmp.qup ) > 0.1 )
  if( length( tmp.ind.qup ) > 0 ){
    cat( "q97.5%: Check the following nodes in chain ", j, ":", names( tmp.qup[tmp.ind.qup] ), "\n" )
    cat(  "Difference: ", tmp.qup[tmp.ind.qup], "\n")
  }
  tmp.qdown    <- Lag_sum$qdown[1,] - Lag_sum$qdown[j,]
  tmp.ind.qdown<- which( abs( tmp.qdown ) > 0.1 )
  if( length( tmp.ind.qdown ) > 0 ){
    cat( "q2.5%: Check the following nodes in chain ", j, ":", names( tmp.qdown[tmp.ind.qdown] ), "\n" )
    cat(  "Difference: ", tmp.qdown[tmp.ind.qdown], "\n")
  }
}
#>END CHECK -- Check node 174 (q2.5%) - run 3 | also run 17
# Lag_sum$qdown[,86:87]
# We decide to delete runs 3 and 17 !
Lag_filt_half1 <- apply( X = Lag_sum$mean[1:17,], MARGIN = 2, FUN = mean )
Lag_filt_half2 <- apply( X = Lag_sum$mean[17:32,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_Lagomorpha.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "Lagomorpha", mean_divt1 = Lag_filt_half1,
                  mean_divt2 = Lag_filt_half2, num_runs = 32 )
dev.off()
Lag_FILT_sum  <- find_prob_MCMC( num_dirs = c(1:2,4:16,18:32), num_chains = 1, delcol = 10,
                                 name_dir_subt = paste( wd2, "Lagomorpha", sep = "" ),
                                 num_divt = 87, subtree = "Lagomorpha_filt",
                                 node_calib = paste( wd2, "Lagomorpha/Calibs_nodes_Lagomorpha.csv", sep = "" ),
                                 perc = 0.975, clean = TRUE, new = TRUE,
                                 out_dat = paste( wd, "out_data", sep = "" ) )
Lag_filt_half1_filt <- apply( X = Lag_FILT_sum$mean[c(1:15),], MARGIN = 2, FUN = mean )
Lag_filt_half2_filt <- apply( X = Lag_FILT_sum$mean[c(16:30),], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_Lagomorpha_filt.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "Lagomorpha", mean_divt1 = Lag_filt_half1_filt,
                  mean_divt2 = Lag_filt_half2_filt, num_runs = 30 )
dev.off() 



#\\\\\\\\\\\\\\\\\\#
# R. CTENOHYSTRICA #
#------------------#
Rcte_sum        <- find_prob_MCMC( num_dirs = 32, num_chains = 1, delcol = 10,
                                   name_dir_subt = paste( wd2, "Ctenohystrica", sep = "" ),
                                   num_divt = 209, subtree = "Ctenohystrica",
                                   node_calib = paste( wd2, "Ctenohystrica/Calibs_nodes_Rctenohystrica.csv", sep = "" ),
                                   perc = 0.975, clean = TRUE, new = TRUE,
                                   out_dat = paste( wd, "out_data", sep = "" ) )

#> CHECK: Issues with quantiles
for ( j in 2:32 ){
  tmp.qup     <- Rcte_sum$qup[1,] - Rcte_sum$qup[j,]
  tmp.ind.qup <- which( abs( tmp.qup ) > 0.1 )
  if( length( tmp.ind.qup ) > 0 ){
    cat( "q97.5%: Check the following nodes in chain ", j, ":", names( tmp.qup[tmp.ind.qup] ), "\n" )
    cat(  "Difference: ", tmp.qup[tmp.ind.qup], "\n")
  }
  tmp.qdown    <- Rcte_sum$qdown[1,] - Rcte_sum$qdown[j,]
  tmp.ind.qdown<- which( abs( tmp.qdown ) > 0.1 )
  if( length( tmp.ind.qdown ) > 0 ){
    cat( "q2.5%: Check the following nodes in chain ", j, ":", names( tmp.qdown[tmp.ind.qdown] ), "\n" )
    cat(  "Difference: ", tmp.qdown[tmp.ind.qdown], "\n")
  }
}
#>END CHECK -- No issues!
Rcte_filt_half1 <- apply( X = Rcte_sum$mean[1:16,], MARGIN = 2, FUN = mean )
Rcte_filt_half2 <- apply( X = Rcte_sum$mean[17:32,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_Ctenohystrica.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "Ctenohystrica", mean_divt1 = Rcte_filt_half1,
                  mean_divt2 = Rcte_filt_half2, num_runs = 32 )
dev.off()


#\\\\\\\\\\\\\#
## EUARCHONTA #
#-------------#
Eua_sum        <- find_prob_MCMC( num_dirs = 32, num_chains = 1, delcol = 10,
                                   name_dir_subt = paste( wd2, "Euarchonta", sep = "" ),
                                   num_divt = 485, subtree = "Euarchonta",
                                   node_calib = paste( wd2, "Euarchonta/Calibs_nodes_Euarchonta.csv", sep = "" ),
                                   perc = 0.975, clean = TRUE, new = TRUE,
                                   out_dat = paste( wd, "out_data", sep = "" ) )

#> CHECK: Issues with quantiles
for ( j in c(1,3:32) ){
  tmp.qup     <- Eua_sum$qup[2,] - Eua_sum$qup[j,]
  tmp.ind.qup <- which( abs( tmp.qup ) > 0.1 )
  if( length( tmp.ind.qup ) > 0 ){
    cat( "q97.5%: Check the following nodes in chain ", j, ":", names( tmp.qup[tmp.ind.qup] ), "\n" )
    cat(  "Difference: ", tmp.qup[tmp.ind.qup], "\n")
  }
  tmp.qdown    <- Eua_sum$qdown[2,] - Eua_sum$qdown[j,]
  tmp.ind.qdown<- which( abs( tmp.qdown ) > 0.1 )
  if( length( tmp.ind.qdown ) > 0 ){
    cat( "q2.5%: Check the following nodes in chain ", j, ":", names( tmp.qdown[tmp.ind.qdown] ), "\n" )
    cat(  "Difference: ", tmp.qdown[tmp.ind.qdown], "\n")
  }
}
#>END CHECK -- Check nodes 487 (ch24), 970 and 971 
#>>Eua_sum$qup[,c(1,484:485)]
# t_n487    t_n970    t_n971
# run1  1.938179 0.9511099 0.2441901 # FLAG
# run2  1.948146 1.2247042 0.4134621
# run3  1.971643 1.0245978 0.2933448 # FLAG
# run4  1.931044 1.1986960 0.4803852
# run5  2.003563 1.1248824 0.4423592
# run6  1.949487 1.2884518 0.5546632 # FLAG
# run7  1.960210 1.2454225 0.4964323
# run8  1.949048 1.1945955 0.4456011
# run9  1.936518 1.1110881 0.3137833 # FLAG
# run10 1.957106 1.1298622 0.3657234
# run11 1.941477 1.2569711 0.4473436
# run12 1.920345 1.1995027 0.4017404
# run13 1.914819 1.2465957 0.4629562
# run14 1.960764 1.2633181 0.4727171
# run15 1.987226 1.1399787 0.3942400
# run16 1.919351 1.2120873 0.4306810
# run17 1.985282 1.1993951 0.4755227
# run18 1.942054 1.2041228 0.4638106
# run19 1.958328 1.1455336 0.3494939
# run20 1.899131 1.2554744 0.5313648 # FLAG
# run21 1.956280 1.0675936 0.4013122 # FLAG
# run22 1.958010 1.2630192 0.4052652
# run23 1.893222 1.2173361 0.4651070
# run24 2.129136 1.0017620 0.2711866 # FLAG
# run25 1.945283 1.2402958 0.5176750 # FLAG
# run26 1.980307 1.1773204 0.3404507
# run27 1.965868 1.1871419 0.4336397
# run28 1.981234 1.2591370 0.4821802
# run29 1.970305 1.1370791 0.3117859 # FLAG
# run30 1.944837 1.2596137 0.5027089
# run31 1.923018 1.3071386 0.5304448 # FLAG 
# run32 1.921758 1.3113743 0.6009471 # FAG 
#>>mean( Eua_sum$qup[,484] )  # 1.188913
#>>mean( Eua_sum$qup[,485] )  # 0.4294537
#>>mean( Eua_sum$qup[1] )     # 1.938179
#
# It seems we should remove chains 1,3,6,9,20,21,24,25,29,31,32
# as there are some issues in q97.5% for these nodes
Eua_filt_half1 <- apply( X = Eua_sum$mean[1:16,], MARGIN = 2, FUN = mean )
Eua_filt_half2 <- apply( X = Eua_sum$mean[17:32,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_Euarchonta.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "Euarchonta", mean_divt1 = Eua_filt_half1,
                  mean_divt2 = Eua_filt_half2, num_runs = 32 )
dev.off()

# Get filtered object
Eua_FILT_sum  <- find_prob_MCMC( num_dirs = c(2,4:5,7:8,10:19,22:23,26:28,30), num_chains = 1, delcol = 10,
                                 name_dir_subt = paste( wd2, "Euarchonta", sep = "" ),
                                 num_divt = 485, subtree = "Euarchonta_filt",
                                 node_calib = paste( wd2, "Euarchonta/Calibs_nodes_Euarchonta.csv", sep = "" ),
                                 perc = 0.975, clean = TRUE, new = TRUE,
                                 out_dat = paste( wd, "out_data", sep = "" ) )
Eua_filt_half1_filt <- apply( X = Eua_FILT_sum$mean[1:11,], MARGIN = 2, FUN = mean )
Eua_filt_half2_filt <- apply( X = Eua_FILT_sum$mean[12:21,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_Euarchonta_filt.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "Euarchonta", mean_divt1 = Eua_filt_half1_filt,
                  mean_divt2 = Eua_filt_half2_filt, num_runs = 21 )
dev.off()


#\\\\\\\\\\\\\\\\\\\\#
# L. CETARTIODACTYLA #
#--------------------#
Lcet_sum      <- find_prob_MCMC( num_dirs = 32, num_chains = 1, delcol = 10,
                                 name_dir_subt = paste( wd2, "Laurasiatheria_cetartiodactyla", sep = "" ),
                                 num_divt = 430, subtree = "Laurasiatheria_cetartiodactyla",
                                 node_calib = paste( wd2, "Laurasiatheria_cetartiodactyla/Calibs_nodes_Lcetartiodactyla.csv", sep = "" ),
                                 perc = 0.975, clean = TRUE, new = TRUE,
                                 out_dat = paste( wd, "out_data", sep = "" ) )

#> CHECK: Issues with quantiles
for ( j in 2:32 ){
  tmp.qup     <- Lcet_sum$qup[1,] - Lcet_sum$qup[j,]
  tmp.ind.qup <- which( abs( tmp.qup ) > 0.1 )
  if( length( tmp.ind.qup ) > 0 ){
    cat( "q97.5%: Check the following nodes in chain ", j, ":", names( tmp.qup[tmp.ind.qup] ), "\n" )
    cat(  "Difference: ", tmp.qup[tmp.ind.qup], "\n")
  }
  tmp.qdown    <- Lcet_sum$qdown[1,] - Lcet_sum$qdown[j,]
  tmp.ind.qdown<- which( abs( tmp.qdown ) > 0.1 )
  if( length( tmp.ind.qdown ) > 0 ){
    cat( "q2.5%: Check the following nodes in chain ", j, ":", names( tmp.qdown[tmp.ind.qdown] ), "\n" )
    cat(  "Difference: ", tmp.qdown[tmp.ind.qdown], "\n")
  }
}
#>END CHECK -- Check 432, 860 and 861, chains 6, 9, and 14
#>>Lcet_sum$qup[,c(1,429:430)]
# t_n432    t_n860     t_n861
# run1  2.171405 0.4245035 0.03334582
# run2  2.123790 0.4329004 0.03382007
# run3  2.189501 0.4476814 0.03422314
# run4  2.138390 0.4293888 0.03300362
# run5  2.113950 0.4352294 0.03330243
# run6  2.411562 0.4979102 0.07265557 #FLAG
# run7  2.105690 0.4195928 0.03334724
# run8  2.100637 0.4412562 0.03318118
# run9  2.164359 0.9062172 0.19050481 #FLAG
# run10 2.173062 0.4364323 0.03352854
# run11 2.107687 0.4474975 0.03336505
# run12 2.203032 0.4435409 0.03347107
# run13 2.085876 0.4322146 0.03344317
# run14 2.178462 0.8968227 0.17164374 #FLAG
# run15 2.228708 0.4302702 0.03293270
# run16 2.160245 0.4380021 0.03319492
# run17 2.116666 0.4248868 0.03306136
# run18 2.102717 0.4164012 0.03336293
# run19 2.108087 0.4399108 0.03399278
# run20 2.187483 0.4488942 0.03429424
# run21 2.140220 0.4352882 0.03319815
# run22 2.145495 0.4307841 0.03344220
# run23 2.160099 0.4460811 0.03356035
# run24 2.169507 0.4351895 0.03379331
# run25 2.141054 0.4248787 0.03334901
# run26 2.145054 0.4341742 0.03348767
# run27 2.125519 0.4416872 0.03349883
# run28 2.153073 0.4217558 0.03333228
# run29 2.168479 0.4244151 0.03296793
# run30 2.118906 0.4133485 0.03347198
# run31 2.202001 0.4244427 0.03369588
# run32 2.140581 0.4227526 0.03333586
# 
# Get rid of chains 6, 9, and 14
Lcet_filt_half1 <- apply( X = Lcet_sum$mean[1:16,], MARGIN = 2, FUN = mean )
Lcet_filt_half2 <- apply( X = Lcet_sum$mean[17:32,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_Lcetartiodactyla.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "L.cetartiodactyla", mean_divt1 = Lcet_filt_half1,
                  mean_divt2 = Lcet_filt_half2, num_runs = 32 )
dev.off()

Lcet_FILT_sum      <- find_prob_MCMC( num_dirs = c(1:5,7:8,10:13,15:32), num_chains = 1, delcol = 10,
                                 name_dir_subt = paste( wd2, "Laurasiatheria_cetartiodactyla", sep = "" ),
                                 num_divt = 430, subtree = "Laurasiatheria_cetartiodactyla_filt",
                                 node_calib = paste( wd2, "Laurasiatheria_cetartiodactyla/Calibs_nodes_Lcetartiodactyla.csv", sep = "" ),
                                 perc = 0.975, clean = TRUE, new = TRUE,
                                 out_dat = paste( wd, "out_data", sep = "" ) )
Lcet_filt_half1_filt <- apply( X = Lcet_FILT_sum$mean[1:15,], MARGIN = 2, FUN = mean )
Lcet_filt_half2_filt <- apply( X = Lcet_FILT_sum$mean[16:29,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_Lcetartiodactyla_filt.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "Artiodactyla", mean_divt1 = Lcet_filt_half1_filt,
                  mean_divt2 = Lcet_filt_half2_filt, num_runs = 29 )
dev.off()



#\\\\\\\\\\\\\\\\\\\\\#
# L. CHIROPTERA SUBT1 #
#---------------------#
Lchir1_sum      <- find_prob_MCMC( num_dirs = 32, num_chains = 1, delcol = 10,
                                   name_dir_subt = paste( wd2, "Laurasiatheria_chiroptera_subt1", sep = "" ),
                                   num_divt = 255, subtree = "Laurasiatheria_chiroptera_subt1",
                                   node_calib = paste( wd2, "Laurasiatheria_chiroptera_subt1/Calibs_nodes_Lchirosubt1.csv", sep = "" ),
                                   perc = 0.975, clean = TRUE, new = TRUE,
                                   out_dat = paste( wd, "out_data", sep = "" ) )

#> CHECK: Issues with quantiles
for ( j in 2:32 ){
  tmp.qup     <- Lchir1_sum$qup[1,] - Lchir1_sum$qup[j,]
  tmp.ind.qup <- which( abs( tmp.qup ) > 0.1 )
  if( length( tmp.ind.qup ) > 0 ){
    cat( "q97.5%: Check the following nodes in chain ", j, ":", names( tmp.qup[tmp.ind.qup] ), "\n" )
    cat(  "Difference: ", tmp.qup[tmp.ind.qup], "\n")
  }
  tmp.qdown    <- Lchir1_sum$qdown[1,] - Lchir1_sum$qdown[j,]
  tmp.ind.qdown<- which( abs( tmp.qdown ) > 0.1 )
  if( length( tmp.ind.qdown ) > 0 ){
    cat( "q2.5%: Check the following nodes in chain ", j, ":", names( tmp.qdown[tmp.ind.qdown] ), "\n" )
    cat(  "Difference: ", tmp.qdown[tmp.ind.qdown], "\n")
  }
}
#>END CHECK -- Check nodes 510 and 511 (Chains 2, 8, 14, 15, 17, 29, 30)
#>>Lchir1_sum$qdown[,254:255]
# t_n510     t_n511
# run1  0.2102302 0.03928404
# run2  0.4337148 0.14172330 # FLAG
# run3  0.2159084 0.04385610
# run4  0.2377703 0.06007100
# run5  0.2221567 0.04974500
# run6  0.2271077 0.05140910
# run7  0.2338950 0.05499221
# run8  0.4990360 0.17529588 # FLAG
# run9  0.2201967 0.04792140
# run10 0.2119941 0.04135309
# run11 0.2169209 0.04673320
# run12 0.2102672 0.03991320
# run13 0.2097102 0.03881766
# run14 0.3987790 0.13283317 # FLAG
# run15 0.4617068 0.15200290 # FLAG
# run16 0.2324512 0.05666854
# run17 0.4935974 0.17367410 # FLAG
# run18 0.2200883 0.04829630
# run19 0.2093086 0.03932080
# run20 0.2127781 0.04093770
# run21 0.2148765 0.04412680
# run22 0.2083565 0.03887900
# run23 0.2798823 0.09349880
# run24 0.2126695 0.04074754
# run25 0.2137631 0.04293567
# run26 0.2110519 0.04027395
# run27 0.2180419 0.04467210
# run28 0.2226099 0.04783174
# run29 0.4889852 0.17086900 # FLAG
# run30 0.4989271 0.17374380 # FLAG
# run31 0.2243496 0.05014722
# run32 0.2211193 0.04750311
#
# It seems we should remove chains 2, 8, 14, 15, 17, 29, 30
Lchir1_filt_half1 <- apply( X = Lchir1_sum$mean[1:16,], MARGIN = 2, FUN = mean )
Lchir1_filt_half2 <- apply( X = Lchir1_sum$mean[17:32,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_Lchirsubt1.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "L.chirosubt1", mean_divt1 = Lchir1_filt_half1,
                  mean_divt2 = Lchir1_filt_half2, num_runs = 32 )
dev.off() 

Lchir1_FILT_sum  <- find_prob_MCMC( num_dirs = c(1,3:7,9:13,16,18:28,31:32), num_chains = 1, delcol = 10,
                                    name_dir_subt = paste( wd2, "Laurasiatheria_chiroptera_subt1", sep = "" ),
                                    num_divt = 255, subtree = "Laurasiatheria_chiroptera_subt1_filt",
                                    node_calib = paste( wd2, "Laurasiatheria_chiroptera_subt1/Calibs_nodes_Lchirosubt1.csv", sep = "" ),
                                    perc = 0.975, clean = TRUE, new = TRUE,
                                    out_dat = paste( wd, "out_data", sep = "" ) )
Lchir1_filt_half1_filt <- apply( X = Lchir1_FILT_sum$mean[1:13,], MARGIN = 2, FUN = mean )
Lchir1_filt_half2_filt <- apply( X = Lchir1_FILT_sum$mean[14:25,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_Lchirsubt1_filt.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "L.chirosubt1", mean_divt1 = Lchir1_filt_half1_filt,
                  mean_divt2 = Lchir1_filt_half2_filt, num_runs = 25 )
dev.off() 


#\\\\\\\\\\\\\\\\\\\\\#
# L. CHIROPTERA SUBT2 #
#---------------------#
Lchir2_sum        <- find_prob_MCMC( num_dirs = 32, num_chains = 1, delcol = 10,
                                     name_dir_subt = paste( wd2, "Laurasiatheria_chiroptera_subt2", sep = "" ),
                                     num_divt = 633, subtree = "Laurasiatheria_chiroptera_subt2",
                                     node_calib = paste( wd2, "Laurasiatheria_chiroptera_subt2/Calibs_nodes_Lchirosubt2.csv", sep = "" ),
                                     perc = 0.975, clean = TRUE, new = TRUE,
                                     out_dat = paste( wd, "out_data", sep = "" ) )

#> CHECK: Issues with quantiles
for ( j in c(1,3:32) ){
  tmp.qup     <- Lchir2_sum$qup[2,] - Lchir2_sum$qup[j,]
  tmp.ind.qup <- which( abs( tmp.qup ) > 0.1 )
  if( length( tmp.ind.qup ) > 0 ){
    cat( "q97.5%: Check the following nodes in chain ", j, ":", names( tmp.qup[tmp.ind.qup] ), "\n" )
    cat(  "Difference: ", tmp.qup[tmp.ind.qup], "\n")
  }
  tmp.qdown    <- Lchir2_sum$qdown[2,] - Lchir2_sum$qdown[j,]
  tmp.ind.qdown<- which( abs( tmp.qdown ) > 0.1 )
  if( length( tmp.ind.qdown ) > 0 ){
    cat( "q2.5%: Check the following nodes in chain ", j, ":", names( tmp.qdown[tmp.ind.qdown] ), "\n" )
    cat(  "Difference: ", tmp.qdown[tmp.ind.qdown], "\n")
  }
}
#>END CHECK -- Check nodes 635 (chains 1,3,8,16)
#>>Lchir2_sum$qup[,1]
#     run1     run2     run3     run4     run5     run6     run7     run8     run9    run10    run11    run12    run13 
# 1.924783 2.048294 2.251309 2.109359 2.035339 2.058640 1.960736 1.948202 2.015431 2.030492 2.026201 2.104283 2.024430 
# run14    run15    run16    run17    run18    run19    run20    run21    run22    run23    run24    run25    run26 
# 2.012282 1.997706 2.186448 2.052785 1.999196 2.103862 2.010054 1.969700 2.012531 2.092393 2.111929 1.979143 2.023125 
# run27    run28    run29    run30    run31    run32 
# 2.034702 2.078196 2.033466 2.113430 2.065556 1.959560 
#
# Remove chains 1, 3, 8, 16
Lchir2_filt_half1 <- apply( X = Lchir2_sum$mean[1:16,], MARGIN = 2, FUN = mean )
Lchir2_filt_half2 <- apply( X = Lchir2_sum$mean[17:32,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_Lchirsubt2.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "L.chirosubt2", mean_divt1 = Lchir2_filt_half1,
                  mean_divt2 = Lchir2_filt_half2, num_runs = 32 )
dev.off() 

Lchir2_FILT_sum  <- find_prob_MCMC( num_dirs = c(2,4:7,9:15,17:32), num_chains = 1, delcol = 10,
                                    name_dir_subt = paste( wd2, "Laurasiatheria_chiroptera_subt2", sep = "" ),
                                    num_divt = 633, subtree = "Laurasiatheria_chiroptera_subt2_filt",
                                    node_calib = paste( wd2, "Laurasiatheria_chiroptera_subt2/Calibs_nodes_Lchirosubt2.csv", sep = "" ),
                                    perc = 0.975, clean = TRUE, new = TRUE,
                                    out_dat = paste( wd, "out_data", sep = "" ) )
Lchir2_filt_half1_filt <- apply( X = Lchir2_FILT_sum$mean[1:14,], MARGIN = 2, FUN = mean )
Lchir2_filt_half2_filt <- apply( X = Lchir2_FILT_sum$mean[15:28,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_Lchirsubt2_filt.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "L.chirosubt2", mean_divt1 = Lchir2_filt_half1_filt,
                  mean_divt2 = Lchir2_filt_half2_filt, num_runs = 28 )
dev.off() 


#\\\\\\\\\\\\\#
# L. THE REST #
#-------------#
Ltherest_sum <- find_prob_MCMC( num_dirs = 32, num_chains = 1, delcol = 10,
                                name_dir_subt = paste( wd2, "Laurasiatheria_therest", sep = "" ),
                                num_divt = 658, subtree = "Laurasiatheria_therest",
                                node_calib = paste( wd2, "Laurasiatheria_therest/Calibs_nodes_Ltherest.csv", sep = "" ),
                                perc = 0.975, clean = TRUE, new = TRUE,
                                out_dat = paste( wd, "out_data", sep = "" ) )

#> CHECK: Issues with quantiles
for ( j in 2:32 ){
  tmp.qup     <- Ltherest_sum$qup[1,] - Ltherest_sum$qup[j,]
  tmp.ind.qup <- which( abs( tmp.qup ) > 0.1 )
  if( length( tmp.ind.qup ) > 0 ){
    cat( "q97.5%: Check the following nodes in chain ", j, ":", names( tmp.qup[tmp.ind.qup] ), "\n" )
    cat(  "Difference: ", tmp.qup[tmp.ind.qup], "\n")
  }
  tmp.qdown    <- Ltherest_sum$qdown[1,] - Ltherest_sum$qdown[j,]
  tmp.ind.qdown<- which( abs( tmp.qdown ) > 0.1 )
  if( length( tmp.ind.qdown ) > 0 ){
    cat( "q2.5%: Check the following nodes in chain ", j, ":", names( tmp.qdown[tmp.ind.qdown] ), "\n" )
    cat(  "Difference: ", tmp.qdown[tmp.ind.qdown], "\n")
  }
}
#>END CHECK -- -- Check nodes 660, 1316
#>>mean(Ltherest_sum$qup[,1])   # 2.186838
#>>mean(Ltherest_sum$qup[,657]) # 0.437407
#>>Ltherest_sum$qup[,c(1,657)]
# t_n660   t_n1316
# run1  2.225667 0.4307261
# run2  2.263907 0.3955360
# run3  2.183523 0.3893073
# run4  2.141647 0.6785645 #FLAG
# run5  2.131324 0.4112076
# run6  2.096538 0.3848733 #minor FLAG
# run7  1.986810 0.3880225 #FLAG
# run8  2.312993 0.5704267 #FLAG
# run9  2.100390 0.4024216 #minorFLAG
# run10 2.073261 0.7870835 #FLAG
# run11 2.185124 0.3961284
# run12 2.176709 0.4192948
# run13 2.337228 0.4193877 #FLAG
# run14 2.089918 0.4013724
# run15 2.160952 0.3961232
# run16 2.403804 0.4238602 #FLAG
# run17 2.178006 0.4086121
# run18 2.263541 0.4559471
# run19 2.167255 0.4222983
# run20 2.153214 0.4538322
# run21 2.246202 0.4015642
# run22 2.345491 0.4167805 #FLAG
# run23 2.182697 0.4044650
# run24 2.192947 0.4436916
# run25 2.105792 0.4063519 #minorFLAG
# run26 2.471961 0.4183441 #FLAG
# run27 2.315073 0.4077386 #FLAG
# run28 2.089845 0.4187536 #minorFLAG
# run29 2.030113 0.4260218 #FLAG
# run30 2.170818 0.3971042
# run31 2.094361 0.4096207 #minorFLAG
# run32 2.101702 0.4115612 #minorFLAG
# 
# We should remove 29,27,26,22,16,13,10,8,7,4
Ltherest_filt_half1 <- apply( X = Ltherest_sum$mean[1:16,], MARGIN = 2, FUN = mean )
Ltherest_filt_half2 <- apply( X = Ltherest_sum$mean[17:32,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_Ltherest.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "Ltherest_subt2", mean_divt1 = Ltherest_filt_half1,
                  mean_divt2 = Ltherest_filt_half2, num_runs = 32 )
dev.off() 

Ltherest_FILT_sum <- find_prob_MCMC( num_dirs = c(1:3,5:6,9,11:12,14:15,17:21,23:25,28,31:32), num_chains = 1, delcol = 10,
                                     name_dir_subt = paste( wd2, "Laurasiatheria_therest", sep = "" ),
                                     num_divt = 658, subtree = "Laurasiatheria_therest_filt",
                                     node_calib = paste( wd2, "Laurasiatheria_therest/Calibs_nodes_Ltherest.csv", sep = "" ),
                                     perc = 0.975, clean = TRUE, new = TRUE,
                                     out_dat = paste( wd, "out_data", sep = "" ) )
Ltherest_filt_half1_filt <- apply( X = Ltherest_FILT_sum$mean[1:11,], MARGIN = 2, FUN = mean )
Ltherest_filt_half2_filt <- apply( X = Ltherest_FILT_sum$mean[12:21,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_Ltherest_filt.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "Ltherest_subt2", mean_divt1 = Ltherest_filt_half1_filt,
                  mean_divt2 = Ltherest_filt_half2_filt, num_runs = 21 )
dev.off() 


#\\\\\\\\\\\\\#
# MARSUPIALIA #
#-------------#
Mar_sum <- find_prob_MCMC( num_dirs = 32, num_chains = 1, delcol = 10,
                           name_dir_subt = paste( wd2, "Marsupialia", sep = "" ),
                           num_divt = 306, subtree = "Marsupialia",
                           node_calib = paste( wd2, "Marsupialia/Calibs_nodes_Marsupialia.csv", sep = "" ),
                           perc = 0.975, clean = TRUE, new = TRUE,
                           out_dat = paste( wd, "out_data", sep = "" ) )

#> CHECK: Issues with quantiles
for ( j in c(1:2, 4:32) ){
  tmp.qup     <- Mar_sum$qup[3,] - Mar_sum$qup[j,]
  tmp.ind.qup <- which( abs( tmp.qup ) > 0.1 )
  if( length( tmp.ind.qup ) > 0 ){
    cat( "q97.5%: Check the following nodes in chain ", j, ":", names( tmp.qup[tmp.ind.qup] ), "\n" )
    cat(  "Difference: ", tmp.qup[tmp.ind.qup], "\n")
  }
  tmp.qdown    <- Mar_sum$qdown[3,] - Mar_sum$qdown[j,]
  tmp.ind.qdown<- which( abs( tmp.qdown ) > 0.1 )
  if( length( tmp.ind.qdown ) > 0 ){
    cat( "q2.5%: Check the following nodes in chain ", j, ":", names( tmp.qdown[tmp.ind.qdown] ), "\n" )
    cat(  "Difference: ", tmp.qdown[tmp.ind.qdown], "\n")
  }
}
#>END CHECK -- Check node 612
#>> Mar_sum$qup[,305]
#>> mean( Mar_sum$qup[,305] ) # 1.187398
# run1     run2     run3     run4     run5     run6     run7     run8     run9    run10    run11    run12    run13 
# 1.095507 1.275830 1.229224 1.214567 1.175542 1.197641 1.221704 1.176550 1.188157 1.132783 1.218467 1.128640 1.240312 
# run14    run15    run16    run17    run18    run19    run20    run21    run22    run23    run24    run25    run26 
# 1.223386 1.165080 1.186361 1.183855 1.211758 1.134489 1.158973 1.271129 1.242588 1.157243 1.162130 1.176051 1.179549 
# run27    run28    run29    run30    run31    run32 
# 1.142137 1.190902 1.174010 1.220614 1.207928 1.113616 
# Minor flags, no chain removed
Mar_filt_half1 <- apply( X = Mar_sum$mean[1:16,], MARGIN = 2, FUN = mean )
Mar_filt_half2 <- apply( X = Mar_sum$mean[17:32,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_Marsupialia.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "Marsupialia", mean_divt1 = Mar_filt_half1,
                  mean_divt2 = Mar_filt_half2, num_runs = 32 )
dev.off()


#\\\\\\\\\\\\\\\#
# ROD. SQUIRREL #
#---------------#
Rsq_sum <- find_prob_MCMC( num_dirs = 32, num_chains = 1, delcol = 10,
                           name_dir_subt = paste( wd2, "Rodentia_squirrel", sep = "" ),
                           num_divt = 266, subtree = "Rodentia_squirrel",
                           node_calib = paste( wd2, "Rodentia_squirrel/Calibs_nodes_Rsquirrel.csv", sep = "" ),
                           perc = 0.975, clean = TRUE, new = TRUE,
                           out_dat = paste( wd, "out_data", sep = "" ) )

#> CHECK: Issues with quantiles
for ( j in 2:32 ){
  tmp.qup     <- Rsq_sum$qup[1,] - Rsq_sum$qup[j,]
  tmp.ind.qup <- which( abs( tmp.qup ) > 0.1 )
  if( length( tmp.ind.qup ) > 0 ){
    cat( "q97.5%: Check the following nodes in chain ", j, ":", names( tmp.qup[tmp.ind.qup] ), "\n" )
    cat(  "Difference: ", tmp.qup[tmp.ind.qup], "\n")
  }
  tmp.qdown    <- Rsq_sum$qdown[1,] - Rsq_sum$qdown[j,]
  tmp.ind.qdown<- which( abs( tmp.qdown ) > 0.1 )
  if( length( tmp.ind.qdown ) > 0 ){
    cat( "q2.5%: Check the following nodes in chain ", j, ":", names( tmp.qdown[tmp.ind.qdown] ), "\n" )
    cat(  "Difference: ", tmp.qdown[tmp.ind.qdown], "\n")
  }
}
#>END CHECK -- No issues!
Rsq_filt_half1 <- apply( X = Rsq_sum$mean[1:16,], MARGIN = 2, FUN = mean )
Rsq_filt_half2 <- apply( X = Rsq_sum$mean[17:32,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_Rsquirrel.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "R.squirrel", mean_divt1 = Rsq_filt_half1,
                  mean_divt2 = Rsq_filt_half2, num_runs = 32 )
dev.off()



#\\\\\\\\\\\\\\\\#
# RODENTIA SUBT1 #
#----------------#
Rod1_sum        <- find_prob_MCMC( num_dirs = 32, num_chains = 1, delcol = 10,
                                   name_dir_subt = paste( wd2, "Rodentia_subtree1", sep = "" ),
                                   num_divt = 629, subtree = "Rodentia_subtree1",
                                   node_calib = paste( wd2, "Rodentia_subtree1/Calibs_nodes_Rtherestsubt1.csv", sep = "" ),
                                   perc = 0.975, clean = TRUE, new = TRUE,
                                   out_dat = paste( wd, "out_data", sep = "" ) )

#> CHECK: Issues with quantiles
for ( j in 2:32 ){
  tmp.qup     <- Rod1_sum$qup[1,] - Rod1_sum$qup[j,]
  tmp.ind.qup <- which( abs( tmp.qup ) > 0.1 )
  if( length( tmp.ind.qup ) > 0 ){
    cat( "q97.5%: Check the following nodes in chain ", j, ":", names( tmp.qup[tmp.ind.qup] ), "\n" )
    cat(  "Difference: ", tmp.qup[tmp.ind.qup], "\n")
  }
  tmp.qdown    <- Rod1_sum$qdown[1,] - Rod1_sum$qdown[j,]
  tmp.ind.qdown<- which( abs( tmp.qdown ) > 0.1 )
  if( length( tmp.ind.qdown ) > 0 ){
    cat( "q2.5%: Check the following nodes in chain ", j, ":", names( tmp.qdown[tmp.ind.qdown] ), "\n" )
    cat(  "Difference: ", tmp.qdown[tmp.ind.qdown], "\n")
  }
}
#>END CHECK -- Check nodes 1258, 1259, and and 631 (check chains 6, 8, 9, 11, 24, 26)
#>>Rod1_sum$qup[,c(1,628:629)]
#>>mean(Rod1_sum$qup[,1]) # 2.002539
# t_n631   t_n1258    t_n1259
# run1  1.993509 0.5089442 0.06889407
# run2  1.986677 0.4745428 0.06342703
# run3  1.965037 0.4806053 0.06698000
# run4  1.928109 0.4800277 0.05650141
# run5  2.035802 0.4424235 0.06098179
# run6  2.104896 0.4800854 0.06477308 #minorFLAG
# run7  2.010490 0.4662140 0.06532791
# run8  2.192369 0.5135035 0.06404995 #FLAG
# run9  1.994642 0.7987016 0.24205050 #FLAG
# run10 2.084739 0.5069069 0.06546588
# run11 2.009673 0.6440756 0.06432784 #FLAG
# run12 1.968975 0.5202663 0.06654988
# run13 2.018705 0.5206900 0.06428834
# run14 2.024224 0.4978588 0.06636380
# run15 2.008385 0.4877355 0.06057077
# run16 1.916874 0.4867055 0.05980248
# run17 1.982954 0.5260911 0.06246453
# run18 1.939622 0.4821591 0.06373859
# run19 2.018655 0.4323503 0.06315555
# run20 1.959108 0.4515949 0.06162420
# run21 2.004343 0.4430664 0.05803009
# run22 1.973040 0.5232832 0.06258525
# run23 1.972366 0.5612533 0.06170675
# run24 2.183535 0.4604043 0.06234533 #FLAG
# run25 1.972228 0.5081107 0.05679595
# run26 2.163885 0.4599946 0.06458364 #FLAG
# run27 2.060494 0.4439007 0.06430816
# run28 1.963925 0.4853832 0.06237924
# run29 2.024140 0.4752104 0.05952564
# run30 1.962845 0.4676016 0.06309357
# run31 1.955132 0.5073369 0.07006046
# run32 2.011702 0.4316273 0.06142246
#
# We decide to remove chains 8,9,11,24,26
# seems to be convergence issues for the highlighted nodes
Rod1_filt_half1 <- apply( X = Rod1_sum$mean[1:16,], MARGIN = 2, FUN = mean )
Rod1_filt_half2 <- apply( X = Rod1_sum$mean[17:32,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_Rodsubt1.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "Rod_subt1", mean_divt1 = Rod1_filt_half1,
                  mean_divt2 = Rod1_filt_half2, num_runs = 32 )
dev.off() 

Rod1_FILT_sum        <- find_prob_MCMC( num_dirs = c(1:7,10,12:23,25,27:32), num_chains = 1, delcol = 10,
                                        name_dir_subt = paste( wd2, "Rodentia_subtree1", sep = "" ),
                                        num_divt = 629, subtree = "Rodentia_subtree1_filt",
                                        node_calib = paste( wd2, "Rodentia_subtree1/Calibs_nodes_Rtherestsubt1.csv", sep = "" ),
                                        perc = 0.975, clean = TRUE, new = TRUE,
                                        out_dat = paste( wd, "out_data", sep = "" ) )
Rod1_filt_half1_filt <- apply( X = Rod1_FILT_sum$mean[1:14,], MARGIN = 2, FUN = mean )
Rod1_filt_half2_filt <- apply( X = Rod1_FILT_sum$mean[15:27,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_Rodsubt1_filt.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "Rod_subt1", mean_divt1 = Rod1_filt_half1_filt,
                  mean_divt2 = Rod1_filt_half2_filt, num_runs = 27 )
dev.off() 

#\\\\\\\\\\\\\\\\#
# RODENTIA SUBT2 #
#----------------#
Rod2_sum        <- find_prob_MCMC( num_dirs = 32, num_chains = 1, delcol = 10,
                                   name_dir_subt = paste( wd2, "Rodentia_subtree2", sep = "" ),
                                   num_divt = 690, subtree = "Rodentia_subtree2",
                                   node_calib = paste( wd2, "Rodentia_subtree2/Calibs_nodes_Rtherestsubt2.csv", sep = "" ),
                                   perc = 0.975, clean = TRUE, new = TRUE,
                                   out_dat = paste( wd, "out_data", sep = "" ) )

#> CHECK: Issues with quantiles
for ( j in 2:32 ){
  tmp.qup     <- Rod2_sum$qup[1,] - Rod2_sum$qup[j,]
  tmp.ind.qup <- which( abs( tmp.qup ) > 0.1 )
  if( length( tmp.ind.qup ) > 0 ){
    cat( "q97.5%: Check the following nodes in chain ", j, ":", names( tmp.qup[tmp.ind.qup] ), "\n" )
    cat(  "Difference: ", tmp.qup[tmp.ind.qup], "\n")
  }
  tmp.qdown    <- Rod2_sum$qdown[1,] - Rod2_sum$qdown[j,]
  tmp.ind.qdown<- which( abs( tmp.qdown ) > 0.1 )
  if( length( tmp.ind.qdown ) > 0 ){
    cat( "q2.5%: Check the following nodes in chain ", j, ":", names( tmp.qdown[tmp.ind.qdown] ), "\n" )
    cat(  "Difference: ", tmp.qdown[tmp.ind.qdown], "\n")
  }
}
#>END CHECK -- Check nodes 692, 1380, 1381 (31,30,29,28,)
#>>mean(Rod2_sum$qup[,1])   # 1.999782
#>>mean(Rod2_sum$qup[,689]) # 1.174721
#>>mean(Rod2_sum$qup[,690]) # 0.3957367
#>>Rod2_sum$qup[,c(1,689:690)]
# t_n692  t_n1380   t_n1381
# run1  1.965029 0.855048 0.1994138 #FLAG
# run2  1.977064 1.144278 0.2633044 #FLAG
# run3  2.065265 1.233980 0.5067273 #FLAG
# run4  2.010376 1.067954 0.2382681 #FLAG
# run5  2.033788 1.199060 0.3530274
# run6  1.965708 1.238936 0.3799263
# run7  2.041261 1.298408 0.5734956 #FALG
# run8  1.960196 1.198066 0.4595468
# run9  2.045990 1.212336 0.4268513
# run10 2.000020 1.200428 0.2113327 #FLAG
# run11 1.969235 1.080392 0.3911291 
# run12 1.946113 1.233165 0.4582122 
# run13 2.016982 0.981314 0.2306939 #FLAG
# run14 2.205339 1.315827 0.5842551 #FLAG
# run15 1.898208 1.254307 0.3654984 #MINFLAG
# run16 1.988952 1.205946 0.4067286
# run17 1.893988 1.173317 0.4757131
# run18 2.019466 1.106654 0.2460795 #FLAG
# run19 1.975426 1.127866 0.3427578
# run20 2.045458 1.246087 0.5086278 #FLAG
# run21 2.002541 1.039606 0.3572143 #FLAG
# run22 2.090224 1.214881 0.4208013
# run23 2.002488 1.052082 0.2734093 #FLAG
# run24 1.962167 1.172585 0.4362453
# run25 1.968436 1.279485 0.6605014 #FLAG
# run26 2.007749 1.315582 0.5884110 #FLAG
# run27 1.941942 1.315407 0.4404022 #FLAG
# run28 2.000316 1.372240 0.6933800 #FLAG
# run29 2.075566 1.056583 0.2476310 #FLAG
# run30 1.967542 1.101254 0.2701617 #FLAG
# run31 2.011743 1.292423 0.4471973 #FLAG
# run32 1.938455 1.005567 0.2066288 #FLAG
#
# keep> 5,6,8,9,11,12,16,17,19,22,24
Rod2_filt_half1 <- apply( X = Rod2_sum$mean[1:16,], MARGIN = 2, FUN = mean )
Rod2_filt_half2 <- apply( X = Rod2_sum$mean[17:32,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_Rodsubt2.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "Rod_subt2", mean_divt1 = Rod2_filt_half1,
                  mean_divt2 = Rod2_filt_half2, num_runs = 32 )
dev.off() 

Rod2_FILT_sum  <- find_prob_MCMC( num_dirs = c(5,6,8,9,11,12,16,17,19,22,24), num_chains = 1, delcol = 10,
                                  name_dir_subt = paste( wd2, "Rodentia_subtree2", sep = "" ),
                                  num_divt = 690, subtree = "Rodentia_subtree2_filt",
                                  node_calib = paste( wd2, "Rodentia_subtree2/Calibs_nodes_Rtherestsubt2.csv", sep = "" ),
                                  perc = 0.975, clean = TRUE, new = TRUE,
                                  out_dat = paste( wd, "out_data", sep = "" ) )
Rod2_filt_half1_filt <- apply( X = Rod2_FILT_sum$mean[1:6,], MARGIN = 2, FUN = mean )
Rod2_filt_half2_filt <- apply( X = Rod2_FILT_sum$mean[7:11,], MARGIN = 2, FUN = mean )
pdf( "plots/Convergence_plot_Rodsubt2_filt.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "Rod_subt2", mean_divt1 = Rod2_filt_half1_filt,
                  mean_divt2 = Rod2_filt_half2_filt, num_runs = 11 )
dev.off() 


#----#

pdf( "plots/Convergence_plot_all_subtrees_filt.pdf" )
# bottom, left, top, right
par( mfrow = c(7,2), mai=c(0.5,0.6,0.4,0.3) )
plot_convergence( name_dir_subt = "Afrotheria", mean_divt1 = Afro_filt_half1,
                  mean_divt2 = Afro_filt_half2, num_runs = 32 )
plot_convergence( name_dir_subt = "Artiodactyla", mean_divt1 = Lcet_filt_half1_filt,
                  mean_divt2 = Lcet_filt_half2_filt, num_runs = 29 )
plot_convergence( name_dir_subt = "Chiroptera-I", mean_divt1 = Lchir1_filt_half1_filt,
                  mean_divt2 = Lchir1_filt_half2_filt, num_runs = 25 )
plot_convergence( name_dir_subt = "Chiroptera-II", mean_divt1 = Lchir2_filt_half1_filt,
                  mean_divt2 = Lchir2_filt_half2_filt, num_runs = 28 )
plot_convergence( name_dir_subt = "Ctenohystrica", mean_divt1 = Rcte_filt_half1,
                  mean_divt2 = Rcte_filt_half2, num_runs = 32 )
plot_convergence( name_dir_subt = "Euarchonta", mean_divt1 = Eua_filt_half1_filt,
                  mean_divt2 = Eua_filt_half2_filt, num_runs = 21 )
plot_convergence( name_dir_subt = "Lagomorpha", mean_divt1 = Lag_filt_half1_filt,
                  mean_divt2 = Lag_filt_half2_filt, num_runs = 30 )
plot_convergence( name_dir_subt = "Marsupialia", mean_divt1 = Mar_filt_half1,
                  mean_divt2 = Mar_filt_half2, num_runs = 32 )
plot_convergence( name_dir_subt = "Rest of Laurasiatheria", mean_divt1 = Ltherest_filt_half1_filt,
                  mean_divt2 = Ltherest_filt_half2_filt, num_runs = 21 )
plot_convergence( name_dir_subt = "Rest of Rodentia - I", mean_divt1 = Rod1_filt_half1_filt,
                  mean_divt2 = Rod1_filt_half2_filt, num_runs = 27 )
plot_convergence( name_dir_subt = "Rest of Rodentia - II", mean_divt1 = Rod2_filt_half1_filt,
                  mean_divt2 = Rod2_filt_half2_filt, num_runs = 11 )
plot_convergence( name_dir_subt = "Sciuridae", mean_divt1 = Rsq_filt_half1,
                  mean_divt2 = Rsq_filt_half2, num_runs = 32 )
plot_convergence( name_dir_subt = "Xenarthra", mean_divt1 = Xen_filt_half1,
                  mean_divt2 = Xen_filt_half2, num_runs = 32 )
dev.off()

#---------------#
# CALCULATE ESS #
#---------------#
#> ESS with RStan
# Each column is assumed to be an MCMC. Rows are iterations for parameter X
# Source explaining why it is preferable than the function in coda:
# https://nature.berkeley.edu/~pdevalpine/MCMC_comparisons/nimble_MCMC_comparisons.html

#\\\\\\\\\\\\#
# AFROTHERIA #
#------------#
ESS_Afro     <- sum_MCMC_ESS( x = Afro_sum$all_mcmc, coda_fun = TRUE )
ESS_Afro$tab
# Tail-ESS Bulk-ESS  coda-ESS
# Med.   130906    68429 138178.62
# Min.    27957    16209  29675.23
# Max.   278303   244573 519552.17
min( ESS_Afro$stats$Rhat )
#[1] 0.9999971
dim( Afro_sum$all_mcmc )
#[1] 640032     59

#\\\\\\\\\\\\#
# EUARCHONTA #
#------------#
ESS_Eua_FILT      <- sum_MCMC_ESS( x = Eua_FILT_sum$all_mcmc, coda_fun = TRUE )
ESS_Eua_FILT$tab
# Tail-ESS Bulk-ESS   coda-ESS
# Med.     4651     1527  3397.0280
# Min.      305      167   475.9892
# Max.    17309    14074 27387.7487
min( ESS_Eua_FILT$stats$Rhat )
#[1] 0.9999437
dim( Eua_FILT_sum$all_mcmc )
#[1] 35507   485

#\\\\\\\\\\\\#
# LAGOMORPHA #
#------------#
ESS_Lag_FILT      <- sum_MCMC_ESS( x = Lag_FILT_sum$all_mcmc, coda_fun = TRUE )
ESS_Lag_FILT$tab
# Tail-ESS Bulk-ESS  coda-ESS
# Med.    36837    11006  22876.28
# Min.     9338     4679  10150.63
# Max.   251588   136426 263722.45
min( ESS_Lag_FILT$stats$Rhat )
#[1] 0.9999968
dim( Lag_FILT_sum$all_mcmc )
#[1] 600030     87


#\\\\\\\\\\\\\\#
# ARTIODACTYLA #
#--------------#
ESS_Lcet     <- sum_MCMC_ESS( x = Lcet_FILT_sum$all_mcmc, coda_fun = TRUE )
ESS_Lcet$tab
# Tail-ESS Bulk-ESS  coda-ESS
# Med.  14400.5     6102 12622.535
# Min.   1713.0     1066  2397.277
# Max.  51125.0    47870 96251.224
min( ESS_Lcet$stats$Rhat )
#[1] 0.999981
dim( Lcet_FILT_sum$all_mcmc )
#[1] 104802    430

#\\\\\\\\\\\\\\\\#
# CHIROPTERA - I #
#----------------#
ESS_FILT_Lchiro1 <- sum_MCMC_ESS( x = Lchir1_FILT_sum$all_mcmc, coda_fun = TRUE )
ESS_FILT_Lchiro1$tab
# Tail-ESS Bulk-ESS    coda-ESS
# Med.    65955    26290  59813.6131
# Min.      461      329    408.7703
# Max.   211161   172817 343651.1053
min( ESS_FILT_Lchiro1$stats$Rhat )
#[1] 0.9999965
dim( Lchir1_FILT_sum$all_mcmc )
#[1] 493192    255

#\\\\\\\\\\\\\\\\\#
# CHIROPTERA - II #
#-----------------#
ESS_filt_Lchiro2  <- sum_MCMC_ESS( x = Lchir2_FILT_sum$all_mcmc, coda_fun = TRUE )
ESS_filt_Lchiro2$tab
# Tail-ESS Bulk-ESS    coda-ESS
# Med.     1958      419  1622.24135
# Min.      112       19    82.29594
# Max.    10379     9045 17424.76113
min( ESS_filt_Lchiro2$stats$Rhat )
#[1] 0.9999122
dim( Lchir2_FILT_sum$all_mcmc )
#[1] 21725      633

#\\\\\\\\\\\\\\\\\\\\\\\\#
# REST OF LAURASIATHERIA #
#------------------------#
ESS_filt_Ltherest <- sum_MCMC_ESS( x = Ltherest_FILT_sum$all_mcmc, coda_fun = TRUE )
ESS_filt_Ltherest$tab
# Tail-ESS Bulk-ESS   coda-ESS
# Med.     1573    522.5  1222.5574
# Min.      171     68.0   163.8509
# Max.     6221   5870.0 11987.5936
min( ESS_filt_Ltherest$stats$Rhat )
#[1] 0.9998421
dim( Ltherest_FILT_sum$all_mcmc )
#[1] 12572      658


#\\\\\\\\\\\\\#
# MARSUPIALIA #
#-------------#
ESS_Mar      <- sum_MCMC_ESS( x = Mar_sum$all_mcmc, coda_fun = TRUE )
ESS_Mar$tab
# Tail-ESS Bulk-ESS   coda-ESS
# Med.    48352  20549.5  43469.937
# Min.     4047   2824.0   7616.737
# Max.   135255 123593.0 245808.216
min( ESS_Mar$stats$Rhat )
#[1] 0.999993
dim( Mar_sum$all_mcmc )
#[1] 284546 306


#\\\\\\\\\\\\\\\#
# CTENOHYSTRICA #
#---------------#
ESS_Rcte     <- sum_MCMC_ESS( x = Rcte_sum$all_mcmc, coda_fun = TRUE )
ESS_Rcte$tab
# Tail-ESS Bulk-ESS  coda-ESS
# Med.    94398    37585  75039.09
# Min.    13121     5727  11372.26
# Max.   315674   308325 634321.30
min( ESS_Rcte$stats$Rhat )
#[1] 0.9999969
dim( Rcte_sum$all_mcmc )
#[1] 640032    209

#\\\\\\\\\\\\\\\\\\\\\\\#
# SCIURIDAE AND RELATED #
#-----------------------#
ESS_Rsq      <- sum_MCMC_ESS( x = Rsq_sum$all_mcmc, coda_fun = TRUE )
ESS_Rsq$tab
# Tail-ESS Bulk-ESS  coda-ESS
# Med.  74819.5  29476.5  61316.96
# Min.   7537.0   7785.0  13228.96
# Max. 260577.0 201905.0 354876.51
min( ESS_Rsq$stats$Rhat )
#[1] 0.9999969
dim( Rsq_sum$all_mcmc )
#[1] 639984     266

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
# RODENTIA THE REST - SUBT 1 #
#----------------------------#
ESS_FILT_Rsubt1   <- sum_MCMC_ESS( x = Rod1_FILT_sum$all_mcmc, coda_fun = TRUE )
ESS_FILT_Rsubt1$tab
# Tail-ESS Bulk-ESS   coda-ESS
# Med.     1615      487  1145.0901
# Min.      126       74   144.5766
# Max.     7475     6336 13930.6301   
min( ESS_FILT_Rsubt1$stats$Rhat )
#[1] 0.9999267
dim( Rod1_FILT_sum$all_mcmc )
#[1] 26231   629

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
# RODENTIA THE REST - SUBT 2 #
#----------------------------#
ESS_filt_Rsubt2   <- sum_MCMC_ESS( x = Rod2_FILT_sum$all_mcmc, coda_fun = TRUE )
ESS_filt_Rsubt2$tab
# Tail-ESS Bulk-ESS   coda-ESS
# Med.    568.5     94.5  328.90016
# Min.     40.0     14.0   29.06115
# Max.   2471.0   2502.0 4951.25049    
min( ESS_filt_Rsubt2$stats$Rhat )
#[1] 0.999821
dim( Rod2_FILT_sum$all_mcmc )
#[1] 11096   690


#\\\\\\\\\\\#
# XENARTHRA #
#-----------#
ESS_Xen      <- sum_MCMC_ESS( x = Xen_sum$all_mcmc, coda_fun = TRUE )
ESS_Xen$tab
# Tail-ESS Bulk-ESS  coda-ESS
# Med. 149874.5  75307.5 150844.28
# Min.  28838.0  13837.0  30199.57
# Max. 285304.0 226337.0 439261.63
min( ESS_Xen$stats$Rhat )
#[1] 0.999997
dim( Xen_sum$all_mcmc )
#[1] 640032     32

