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
# num_dirs 
# num_chains     Numeric. Default is 2 chains per run. Change if needed.
# delcol         Numeric. Default is 10 so only the samples for time estimates are taken.
# name_dir_subt  Character. Name of the directory where the anlayses for this subtree ran.
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
find_prob_MCMC <- function ( num_dirs = 6, num_chains = 2, delcol = 10, name_dir_subt, num_divt,
                             node_calib, perc = 0.975, clean = FALSE, new = FALSE, new_rev = FALSE ){
  
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
                                                 subt = paste( name_dir_subt, i, "_run", i, sep = "" ),
                                                 delcol = delcol, perc = perc )
        }else if( new == FALSE ){
          if( clean == FALSE ){
            subtree_list[[ count ]]  <- load_subt( mcmc1 = paste( name_dir_subt, "/7/run", i, "/mcmc", j, "/mcmctree_GBM/mcmc.txt", sep = "" ),
                                                   subt = paste( name_dir_subt, i, "_run", j, sep = "" ),
                                                   delcol = delcol, perc = perc )
          }else if( clean == TRUE ){
            subtree_list[[ count ]]  <- load_subt( mcmc1 = paste( name_dir_subt, "/7/run", i, "/mcmc", j, "/mcmctree_GBM/mcmc_clean.txt", sep = "" ),
                                                   subt = paste( name_dir_subt, i, "_run", j, sep = "" ),
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
  write.table( x = subtree_meandivt, file = paste( name_dir_subt, "/mean_divt.tsv", sep = "" ), sep = "\t",
               quote = FALSE )
  write.table( x = subtree_qup, file = paste( name_dir_subt, "/mean_qup.tsv", sep = "" ), sep = "\t",
               quote = FALSE )
  write.table( x = subtree_qlow, file = paste( name_dir_subt, "/mean_qlow.tsv", sep = "" ), sep = "\t",
               quote = FALSE )
  if( new_rev == TRUE ){
    write.table( x = mean_est_priors, file = paste( name_dir_subt, "/all_mean_est.tsv", sep = "" ),
                 sep = "\t", quote = FALSE )
  }else{
    write.table( x = mean_est_priors, file = paste( name_dir_subt, "/", name_dir_subt, "_all_mean_est.tsv", sep = "" ),
               sep = "\t", quote = FALSE )
  }
  
  cat( "Output files available! Check returned list with all objects generated too :) \n\n")
  
  return( list( tt_all = subtree_list, mean = subtree_meandivt, qup = subtree_qup, qdown = subtree_qlow,
                all_mean_est = mean_est, all_mcmc = mcmc_all ) )
  
}

# Function to find problematic runs.
# It prints out the mean div time, qlow, and qup 
# for each node and for each run
# 
# Parameters:
# num_dirs       Numeric. Default is c(1,2,3,5,6,9,10) for first run. Now,
#                default is c(3,4,7,8,9,11,13,16). Change if needed.
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
#-- MAIN TREE (T2, data set 1) --#
# 72sp tree - Load it to generate the plots with the results obtained 
# with the individual data sets for the extra analyses
wd2 <- gsub( pattern = "03_Extra_analyses/..*",
             replacement = "01_SeqBayes_S1/02_MCMCtree/01_ESS_and_chain_convergence/",
             x = wd )
wd3 <- gsub( pattern = "01_ESS_and_chain_convergence/", replacement = "", x = wd2 )
setwd( wd2 )
softbounds_72sp_sum  <- find_prob_MCMC_72sp( num_dirs = c(3,4,7,8,9,11,13,16), delcol = 8, name_dir_subt = "72sp",
                                             num_divt = 71, node_calib = "Calibs_nodes_72sp.csv", 
                                             tree_hyp = "02_atlantogenata_tarver2016", maintt = TRUE,
                                             clock = "GBM", out = "out_data/00_post_72sp_new/",
                                             path_72sp = paste( wd3, "00_MCMCtree_analyses/00_main_tree_T2/01_MCMCtree_posterior/02_atlantogenata_tarver2016/",
                                                                sep = "" ),
                                             perc = 0.975 )

#-- DATA SUBSETS (data set 2) --#
wd3 <- gsub( pattern = "03_Extra_analyses/..*",
             replacement = "02_SeqBayes_S2/02_MCMCtree/",
             x = wd )
setwd( wd3 )

## FUNCTION 
# This function generates a data.frame with matching nodes 
# between data 1 and data 2
#
# Arguments:
# inp_mat     Matrix/Data.frame. This object is passed to the function, which 
#             keeps getting updated after having gone through all data subsets
# tmp_nodes   Numeric, temporary nodes for each data subset (2)
# tmp_pos     Numeric, temporary positions to replace in rows in `inp_mat`.
# sum_list    List, list generated for each data subset
# tmp_nodes72 Numeric, temporary nodes for data 72sp (1)
# sum_list72  List, list generated for data 72sp
sum_tableD1D2 <- function( inp_mat, tmp_nodes, tmp_pos, sum_list, tmp_nodes72, sum_list72, subtree ){
  
  # Check nodes in data 2
  for ( i in 1:length( tmp_nodes ) ){
    
    cat( "4.7K sp:Evaluating node ", tmp_nodes[i], "...\n" ) 
    tmp_nums <- gsub( pattern = "t_n", replacement = "", x = colnames( sum_list$mean) )
    tmp_ind  <- which( tmp_nums %in% tmp_nodes[i] )
    colnames( sum_list$mean )[tmp_ind] 
    inp_mat[tmp_pos[i],4] <- mean( sum_list$mean[,tmp_ind])
    inp_mat[tmp_pos[i],5] <- mean( sum_list$qdown[,tmp_ind])
    inp_mat[tmp_pos[i],6] <- mean( sum_list$qup[,tmp_ind])
    inp_mat[tmp_pos[i],7] <- paste( colnames( sum_list$mean )[tmp_ind], subtree, sep = "|" ) 
    cat( "Mean: ", mean( sum_list$mean[,tmp_ind]), "\n",
         "Q-2.5%:", mean( sum_list$qdown[,tmp_ind]), "\n",
         "Q-97.5%", mean( sum_list$qup[,tmp_ind]), "\n\n" )
    
    # Data 1
    cat( "72 sp:Evaluating node ", tmp_nodes72[i], "...\n" ) 
    tmp_nums72 <- gsub( pattern = "t_n", replacement = "", x = colnames( sum_list72$mean) )
    tmp_ind72  <- which( tmp_nums72 %in% tmp_nodes72[i] )
    colnames( sum_list72$mean )[tmp_ind72] 
    inp_mat[tmp_pos[i],1] <- mean( sum_list72$mean[,tmp_ind72])
    inp_mat[tmp_pos[i],2] <- mean( sum_list72$qdown[,tmp_ind72])
    inp_mat[tmp_pos[i],3] <- mean( sum_list72$qup[,tmp_ind72])
    cat( "Mean: ", mean( sum_list72$mean[,tmp_ind72]), "\n",
         "Q-2.5%:", mean( sum_list72$qdown[,tmp_ind72]), "\n",
         "Q-97.5%", mean( sum_list72$qup[,tmp_ind72]), "\n\n" )
    
    # Add flag
    inp_mat[tmp_pos[i],8] <- 0
    
  }
  
  # Return object
  return( inp_mat )
  
}

# 0. Create lists/vectors needed to find matching nodes
matching_nodes_d1VSd2 <- matrix( 0, ncol = 8, nrow = length( colnames(softbounds_72sp_sum$mean) )+1 )
colnames( matching_nodes_d1VSd2 ) <- c( "mean.data1", "qdown.data1", "qup.data1", 
                                        "mean.data2", "qdown.data1", "qup.data2",
                                        "data2-nodenames", "only_bb" )
names_nodes <- c( "MAMMALIA", "THERIA", "DIDELPHIMORPHIA-AUSTRALIDELPHIA",
                  "EOMETATHERIA", "PLACENTALIA", "XENARTHRA-AFROTHERIA (ATLANTOGENATA)", "XENARTHRA",
                  "AFROTHERIA", "PAENUNGULATA", "BOREOEUTHERIA", "LAURASIATHERIA (Lipotyphla-Carnivora)", 
                  "ERINACEIDAE-SORICIDAE", "SCROTIFERA", "CHIROPTERA", "FEREUNGULATA", "CARNIVORA", 
                  "FELIDAE", "PANTHERINAE", "CANIFORMIA", "ARCTOIDEA", "EUUNGULATA", "ARTIODACTYLA", "ARTIOFABULA",
                  "CETRUMINANTIA", "BOVIDAE", "ovis_aries-capra_hircus", "EUARCHONTOGLIRES", "Tupaia-Glires", 
                  "GLIRES", "LAGOMORPHA", "RODENTIA", "ROD-NOSQUIRREL", "CAVIOMORPHA-PHIOMORPHA", "PHIOMORPHA",
                  "CAVIOMORPHA", "Cavia_porcellus-Cavia_aperea", "Chinchilla_lanigera-Octodon_degus", 
                  "dipodomys_ordi-MYOMORPHA", "MYOMORPHA", "nannospalax_galili-MURIDAE", "MURIDAE", 
                  "CRICETIDAE", "Mesocricetus_auratus-Cricetulus_griseus", "Peromyscus_maniculatus-Microtus_ochrogaster",
                  "MURINAE", "Mus_pahari-rest", "Mus_caroli-rest", "Mus_spretus-Mus_musculus", "PRIMATES", 
                  "STREPSIRHINI", "within_Lemuroidea (propithecus-microcebus)", "HAPLORHINI (Simiiformes-Tarsiiformes)", 
                  "ANTHROPOIDEA/SIMIIFORMES", "AOTIDAE-CALLITRICHIDAE", "aotus_nancymaae - callithrix_jacchus", "CEBIDAE", 
                  "CATARRHINI", "CERCOPITHECOIDEA", "CERCOPITHECINAE", "PAPIONINI", "Papio-Mandrillus",
                  "cercocebus_atys-mandrillus_leucophaeus", "GENUS MACACA", "macaca_fascicularis-macaca_mulatta", "COLOBINAE",
                  "GENUS RHINOPITHECUS", "HOMINOIDEA", "HOMINIDAE", "HOMININAE", "HOMININI", "pan_paniscus-pan_troglodites" )

rownames( matching_nodes_d1VSd2 ) <- c( paste( colnames( softbounds_72sp_sum$mean )[1:3], names_nodes[1:3], sep = " || " ),
                                        "NA    || MARSUPIALIA ",
                                        paste( colnames( softbounds_72sp_sum$mean )[4:length(colnames( softbounds_72sp_sum$mean ))],
                                               names_nodes[4:length(colnames( softbounds_72sp_sum$mean ))], sep = " || " ) )
names_nodes_2match <- c( colnames( softbounds_72sp_sum$mean )[1:3], "NA",
                         colnames( softbounds_72sp_sum$mean )[4:length(colnames( softbounds_72sp_sum$mean ))] )

# rownames( matching_nodes_d1VSd2 ) <- paste( colnames( softbounds_72sp_sum$mean ), names_nodes, sep = " || " )
# rownames( matching_nodes_d1VSd2 ) <- paste( "t_n", c(75, 76, 79, 80, 81, 84, 86, 88, 91,
#                                                      93, 94, 96, 97, 99, 102, 103, 110,
#                                                      117, 121, 122, 125, 129, 131, 132,
#                                                      140, 141, 142 ))
# rownames( matching_nodes_d1VSd2 ) <- colnames( softbounds_72sp_sum$mean )

## 1. MARSUPIALIA
Mar_sum        <- find_prob_MCMC( num_dirs = 6, num_chains = 2, delcol = 10, name_dir_subt = "Marsupialia",
                                  num_divt = 306, node_calib = "Marsupialia/Calibs_nodes_Marsupialia.csv",
                                  perc = 0.975, clean = TRUE, new = FALSE )

# Nodes                            | 72sp   | Marsupialia
# Mammalia                      	 | t_n72  | t_n308
# DIDELPHIMORPHIA-AUSTRALIDELPHIA	 | t_n75  | t_n310
# Marsupialia                      | NA     | t_n309
# Eometatheria                     | t_n76  | t_n312

# Update inp_mat with
tmp_nodes72 <- c(75, 76)
tmp_pos     <- which( names_nodes_2match %in% paste( "t_n", tmp_nodes72, sep = "" ) )
matching_nodes_d1VSd2 <- sum_tableD1D2( inp_mat = matching_nodes_d1VSd2,
                                        tmp_nodes = c(310, 312), tmp_pos = tmp_pos,
                                        sum_list = Mar_sum,
                                        tmp_nodes72 = tmp_nodes72, sum_list72 = softbounds_72sp_sum,
                                        subtree = "Marsupialia" )
# Add manually results for Marsupialia, which come only from 
# the subtree
matching_nodes_d1VSd2[4,1] <- matching_nodes_d1VSd2[4,4] <- apply( Mar_sum$mean, 2, mean )[2]
matching_nodes_d1VSd2[4,2] <- matching_nodes_d1VSd2[4,5] <- apply( Mar_sum$qdown, 2, mean )[2]
matching_nodes_d1VSd2[4,3] <- matching_nodes_d1VSd2[4,6] <- apply( Mar_sum$qup, 2, mean )[2]
matching_nodes_d1VSd2[4,7] <- "t_n309|Marsupialia"
matching_nodes_d1VSd2[4,8] <- 1

## 2. XENARTHRA
Xen_sum        <- find_prob_MCMC( num_dirs = 6, num_chains = 2, delcol = 10, name_dir_subt = "Xenarthra",
                                  num_divt = 32, node_calib = "Xenarthra/Calibs_nodes_Xenarthra.csv",
                                  perc = 0.975, clean = FALSE, new = FALSE )

# Nodes                                       | 72sp   | Rod.subt2
# Mammalia                      	            | t_n72  | t_n34
# Xenarthra                      	            | t_n79  | t_n35

# Update inp_mat
tmp_nodes72 <- 79
tmp_pos     <- which( names_nodes_2match %in% paste( "t_n", tmp_nodes72, sep = "" ) )
matching_nodes_d1VSd2 <- sum_tableD1D2( inp_mat = matching_nodes_d1VSd2,
                                        tmp_nodes = 35, tmp_pos = tmp_pos,
                                        sum_list = Xen_sum,
                                        tmp_nodes72 = tmp_nodes72, sum_list72 = softbounds_72sp_sum,
                                        subtree = "Xenarthra")

## 3. AFROTHERIA 
Afro_sum            <- find_prob_MCMC( num_dirs = 16, num_chains = 1, delcol = 10,
                                       name_dir_subt = "Afrotheria",
                                       num_divt = 59,
                                       node_calib = "Afrotheria/Calibs_nodes_Afrotheria.csv",
                                       perc = 0.975, clean = TRUE, new = TRUE )

# Nodes        | 72sp   | Afrotheria
# Mammalia     | t_n72  | t_n61
# Afrotheria   | t_n80  | t_n62
# Paenungulata | t_n81  | t_n110

# Update inp_mat 
tmp_nodes72 <-  c( 80, 81 )
tmp_pos     <- which( names_nodes_2match %in% paste( "t_n", tmp_nodes72, sep = "" ) )
matching_nodes_d1VSd2 <- sum_tableD1D2( inp_mat = matching_nodes_d1VSd2,
                                        tmp_nodes = c( 62, 110 ), tmp_pos = tmp_pos,
                                        sum_list = Afro_sum,
                                        tmp_nodes72 = tmp_nodes72, sum_list72 = softbounds_72sp_sum,
                                        subtree = "Afrotheria")
## 4. L. THE REST
Ltherest_sum        <- find_prob_MCMC( num_dirs = 32, num_chains = 1, delcol = 10, name_dir_subt = "Laurasiatheria_therest",
                                       num_divt = 658,
                                       node_calib = "Laurasiatheria_therest/Calibs_nodes_Ltherest.csv",
                                       perc = 0.975, clean = TRUE, new = TRUE )

# Nodes                                           | 72sp   | L.therest
# Mammalia                                        | t_n72  | t_n660
# LAURASIATHERIA (Lipotyphla-Carnivora)	          | t_n83  | t_n661
# ERINACEIDAE-SORICIDAE                           | t_n84  | t_n1005
# SCROTIFERA	                                    | t_n85  | t_n662
# CHIROPTERA	                                    | t_n86  | t_n1008 # node used from this tree
# FEREUNGULATA (Carnivora-Euungulata)	            | t_n87  | t_n664
# CARNIVORA	                                      | t_n88  | t_n665
# FELIDAE (Felis-Panthera)                       	| t_n89  | t_n917
# PANTHERINAE (P.trigis-P.pardus)	                | t_n90  | t_n955
# CANIFORMIA	                                    | t_n91  | t_n666
# ARCTOIDEA                                    	  | t_n92  | t_n667
# EUUNGULATA                                     	| t_n93  | t_n969
# ARTIODACTYLA	                                  | t_n94  | t_n1000 # From L.cetartiodactyla

# Update inp_mat  
# NOTE: Skip node 94 -- L.cetartiodactyla
tmp_nodes72 <-  c( c(83:93) )
tmp_pos     <- which( names_nodes_2match %in% paste( "t_n", tmp_nodes72, sep = "" ) )
matching_nodes_d1VSd2 <- sum_tableD1D2( inp_mat = matching_nodes_d1VSd2,
                                        tmp_nodes = c( 661, 1005, 662, 1008, 664, 665, 917, 
                                                       955, 666, 667, 969 ),
                                        tmp_pos = tmp_pos,
                                        sum_list = Ltherest_sum,
                                        tmp_nodes72 = tmp_nodes72,
                                        sum_list72 = softbounds_72sp_sum,
                                        subtree = "L.therest" )

## 5. L. CHIROPTERA SUBT1 -- NODE NOT USED IN MAIN TREE
Lchir1_sum        <- find_prob_MCMC( num_dirs = 32, num_chains = 1, delcol = 10, name_dir_subt = "Laurasiatheria_chiroptera_subt1",
                                     num_divt = 255,
                                     node_calib = "Laurasiatheria_chiroptera_subt1/Calibs_nodes_Lchirosubt1.csv",
                                     perc = 0.975, clean = TRUE, new = TRUE )

# Nodes         | 72sp   | L.chiroptera.subt1
# Mammalia      | t_n72  | t_n256
# Chiroptera    | t_n86  | t_n258

# Check nodes in data2 (use for loop to avoid issues with number ordering)
# tmp_nodes72 <-  86
# tmp_pos     <- which( names_nodes_2match %in% paste( "t_n", tmp_nodes72, sep = "" ) )
# matching_nodes_d1VSd2 <- sum_tableD1D2( inp_mat = matching_nodes_d1VSd2,
#                                         tmp_nodes = 258,
#                                         tmp_pos = tmp_pos,
#                                         sum_list = Lchir1_sum,
#                                         tmp_nodes72 = tmp_nodes72,
#                                         sum_list72 = softbounds_72sp_sum,
#                                         subtree = "Chiroptera")

## 5. L. CHIROPTERA SUBT2 -- NODE NOT USED IN MAIN TREE
Lchir2_sum        <- find_prob_MCMC( num_dirs = 32, num_chains = 1, delcol = 10, name_dir_subt = "Laurasiatheria_chiroptera_subt2",
                                     num_divt = 633,
                                     node_calib = "Laurasiatheria_chiroptera_subt2/Calibs_nodes_Lchirosubt2.csv",
                                     perc = 0.975, clean = TRUE, new = TRUE )

# Nodes         | 72sp   | L.chiroptera.subt2
# Mammalia      | t_n72  | t_n635
# Chiroptera    | t_n86  | t_n636

# Check nodes in data2 (use for loop to avoid issues with number ordering)
# tmp_nodes72 <-  86
# tmp_pos     <- which( colnames( softbounds_72sp_sum$mean ) %in% paste( "t_n", tmp_nodes72, sep = "" ) )
# matching_nodes_d1VSd2 <- sum_tableD1D2( inp_mat = matching_nodes_d1VSd2,
#                                         tmp_nodes = 636,
#                                         tmp_pos = tmp_pos,
#                                         sum_list = Lchir2_sum,
#                                         tmp_nodes72 = tmp_nodes72,
#                                         sum_list72 = softbounds_72sp_sum )

## 6. L. CETARTIODACTYLA
Lcet_sum        <- find_prob_MCMC( num_dirs = 6, num_chains = 2, delcol = 10, name_dir_subt = "Laurasiatheria_cetartiodactyla",
                                   num_divt = 430,
                                   node_calib = "Laurasiatheria_cetartiodactyla/Calibs_nodes_Lcetartiodactyla.csv",
                                   perc = 0.975, clean = TRUE, new = FALSE )

# Nodes                      | 72sp   | L.cetartiodactyla
# Mammalia                   | t_n72  | t_n432
# ARTIODACTYLA               | t_n94  | t_n433
# ARTIOFABULA                | t_n95  | t_n434 
# CETRUMINANTIA (WHIP-RUM)   | t_n96  | t_n435
# BOVIDAE                    | t_n97  | t_n440
# ovis_aries-capra_hircus    | t_n98  | t_n446

tmp_nodes72 <- 94:98
tmp_pos     <- which( names_nodes_2match %in% paste( "t_n", tmp_nodes72, sep = "" ) )
matching_nodes_d1VSd2 <- sum_tableD1D2( inp_mat = matching_nodes_d1VSd2,
                                        tmp_nodes = c(433,434,435,440,446),
                                        tmp_pos = tmp_pos,
                                        sum_list = Lcet_sum,
                                        tmp_nodes72 = tmp_nodes72,
                                        sum_list72 = softbounds_72sp_sum,
                                        subtree = "Artiodactyla" )

## 7. EUARCHONTA
Eua_sum        <- find_prob_MCMC( num_dirs = 6, num_chains = 2, delcol = 10, name_dir_subt = "Euarchonta",
                                  num_divt = 485, node_calib = "Euarchonta/Calibs_nodes_Euarchonta.csv",
                                  perc = 0.975, clean = TRUE, new = FALSE )

# Nodes                                        | 72sp    | Eaurchonta
# Mammalia                                     | t_n72   | t_n487
# Euarchontoglires                             | t_n99   | t_n488
# PRIMATES	                                   | t_n121  | t_n490
# STREPSIRHINI	                               | t_n122  | t_n825
# within_Lemuroidea (propithecus-microcebus)	 | t_n123  | t_n827
# HAPLORHINI (Simiiformes-Tarsiiformes)	       | t_n124  | t_n491
# ANTHROPOIDEA/SIMIIFORMES	                   | t_n125  | t_n492
# AOTIDAE-CALLITRICHIDAE	                     | t_n126  | t_n695
# CEBIDAE (Cebus-Saimiri)	                     | t_n128  | t_n731
# CATARRHINI	                                 | t_n129  | t_n493
# CERCOPITHECOIDEA (Old World Monkeys)         | t_n130  | t_n494
# CERCOPITHECINAE	                             | t_n131  | t_n495
# PAPIONINI		                                 | t_n132  | t_n552
# Papio-Mandrillus	                           | t_n133  | t_n575
# cercocebus_atys-mandrillus_leucophaeus	     | t_n134  | t_n585
# GENUS MACACA                                 | t_n135  | t_n555
# macaca_fascicularis-macaca_mulatta	         | t_n136  | t_n557
# COLOBINAE (Colobus-Rhinopithecus)	           | t_n137  | t_n592
# GENUS RHINOPITHECUS (R.roxellana-R.bieti)    | t_n138  | t_n624
# HOMINOIDEA	                                 | t_n139  | t_n657
# HOMINIDAE (great apes)	                     | t_n140  | t_n683
# HOMININAE	                                   | t_n141  | t_n684
# HOMININI 	                                   | t_n142  | t_n685
# pan_paniscus-pan_troglodites	               | t_n143  | t_n686

tmp_nodes72 <- c(99, 121:126, 128:143)
tmp_pos     <- which( names_nodes_2match %in% paste( "t_n", tmp_nodes72, sep = "" ) )
matching_nodes_d1VSd2 <- sum_tableD1D2( inp_mat = matching_nodes_d1VSd2,
                                        tmp_nodes = c(488, 490, 825, 827, 491, 492, 695, 
                                                      731, 493, 494, 495, 552, 575, 585, 
                                                      555, 557, 592, 624, 657, 683, 684, 
                                                      685, 686 ),
                                        tmp_pos = tmp_pos,
                                        sum_list = Eua_sum,
                                        tmp_nodes72 = tmp_nodes72,
                                        sum_list72 = softbounds_72sp_sum,
                                        subtree = "Euarchonta" )

## 8. LAGOMORPHA
Lag_sum        <- find_prob_MCMC( num_dirs = 1, num_chains = 2, delcol = 10, name_dir_subt = "Lagomorpha",
                                  num_divt = 87, node_calib = "Lagomorpha/Calibs_nodes_Lagomorpha.csv",
                                  perc = 0.975, clean = TRUE )

# Nodes       | 72sp    | Lagomorpha
# Mammalia    | t_n72   | t_n89
# Lagomorpha  | t_n102  | t_n90

tmp_nodes72 <- 102
tmp_pos     <- which( names_nodes_2match %in% paste( "t_n", tmp_nodes72, sep = "" ) )
matching_nodes_d1VSd2 <- sum_tableD1D2( inp_mat = matching_nodes_d1VSd2,
                                        tmp_nodes = 90,
                                        tmp_pos = tmp_pos,
                                        sum_list = Lag_sum,
                                        tmp_nodes72 = tmp_nodes72,
                                        sum_list72 = softbounds_72sp_sum,
                                        subtree = "Lagomorpha" )

## 9. R. CTENOHYSTRICA
Rcte_sum        <- find_prob_MCMC( num_dirs = 6, num_chains = 2, delcol = 10, name_dir_subt = "Rodentia_ctenohystrica",
                                   num_divt = 209,
                                   node_calib = "Rodentia_ctenohystrica/Calibs_nodes_Rctenohystrica.csv",
                                   perc = 0.975, clean = TRUE, new = FALSE )

# Nodes                             | 72sp   | R.ctenohystrica
# Mammalia                      	  | t_n72  | t_n211
# Rodentia                      	  | t_n103 | t_n212
# CAVIOMORPHA-PHIOMORPHA	          | t_n105 | t_n214
# BATHYERGIDAE | PHIOMORPHA         | t_n106 | t_n385
# CAVIOMORPHA                       | t_n107 | t_n215
# Cavia_porcellus-Cavia_aperea	    | t_n108 | t_n360
# Chinchilla_lanigera-Octodon_degus	| t_n109 | t_n3217

tmp_nodes72 <- c(103, 105:109)
tmp_pos     <- which( names_nodes_2match %in% paste( "t_n", tmp_nodes72, sep = "" ) )
matching_nodes_d1VSd2 <- sum_tableD1D2( inp_mat = matching_nodes_d1VSd2,
                                        tmp_nodes = c( 212, 214, 385, 215, 360, 217 ),
                                        tmp_pos = tmp_pos,
                                        sum_list = Rcte_sum,
                                        tmp_nodes72 = tmp_nodes72,
                                        sum_list72 = softbounds_72sp_sum,
                                        subtree = "Ctenohystrica")

## 10. RODENTIA SUBT1
Rod1_sum        <- find_prob_MCMC( num_dirs = 32, num_chains = 1, delcol = 10, name_dir_subt = "Rodentia_subtree1",
                                   num_divt = 629,
                                   node_calib = "Rodentia_subtree1/Calibs_nodes_Rtherestsubt1.csv",
                                   perc = 0.975, clean = TRUE, new = TRUE )

# Nodes                                       | 72sp   | Rod.subt1
# Mammalia                      	            | t_n72  | t_n631
# d.ordi-MYOMORPHA                            | t_n110 | t_n632
# MYOMORPHA	                                  | t_n111 | t_n633
# nannospalax_galili-MURIDAE	                | t_n112 | t_n635
# MURIDAE   	                                | t_n113 | t_n636
# Mesocricetus_auratus-Cricetulus_griseus	    | t_n115 | t_n938
# Peromyscus_maniculatus-Microtus_ochrogaster	| t_n116 | t_n638
# MURINAE	                                    | t_n117 | t_n1011 # from Rod. subt 2 !

tmp_nodes72 <- c(110:113, 115:116)
tmp_pos     <- which( names_nodes_2match %in% paste( "t_n", tmp_nodes72, sep = "" ) )
matching_nodes_d1VSd2 <- sum_tableD1D2( inp_mat = matching_nodes_d1VSd2,
                                        tmp_nodes = c( 632,633,635,636,938,638 ),
                                        tmp_pos = tmp_pos,
                                        sum_list = Rod1_sum,
                                        tmp_nodes72 = tmp_nodes72,
                                        sum_list72 = softbounds_72sp_sum,
                                        subtree = "Rodentia subt1" )

## 10. RODENTIA SUBT2 
Rod2_sum        <- find_prob_MCMC( num_dirs = 32, num_chains = 1, delcol = 10, name_dir_subt = "Rodentia_subtree2",
                                   num_divt = 690,
                                   node_calib = "Rodentia_subtree2/Calibs_nodes_Rtherestsubt2.csv",
                                   perc = 0.975, clean = TRUE, new = TRUE )

# Nodes                                       | 72sp   | Rod.subt2
# Mammalia                      	            | t_n72  | t_n692
# MURIDAE                                     | t_n113 | t_n693 # From Rod. subt 1 !
# CRICETIDAE                                  | t_n114 | t_n1379
# MURINAE (rattus-mus)                        | t_n117 | t_n698
# Mus_pahari-rest                             | t_n118 | t_n752
# Mus_caroli-rest                             | t_n119 | t_n770
# Mus_spretus-Mus_musculus                    | t_n120 | t_n771

tmp_nodes72 <- c(114,117:120)
tmp_pos     <- which( names_nodes_2match %in% paste( "t_n", tmp_nodes72, sep = "" ) )
matching_nodes_d1VSd2 <- sum_tableD1D2( inp_mat = matching_nodes_d1VSd2,
                                        tmp_nodes = c( 1379,698,752,770,771 ),
                                        tmp_pos = tmp_pos,
                                        sum_list = Rod2_sum,
                                        tmp_nodes72 = tmp_nodes72,
                                        sum_list72 = softbounds_72sp_sum,
                                        subtree = "Rodentia subt2" )

## 11. ROD. SQUIRREL -- NO MAIN NODES IN MAIN TREE
Rsq_sum        <- find_prob_MCMC( num_dirs = 6, num_chains = 2, delcol = 10, name_dir_subt = "Rodentia_squirrel",
                                  num_divt = 266,
                                  node_calib = "Rodentia_squirrel/Calibs_nodes_Rsquirrel.csv",
                                  perc = 0.975, clean = TRUE, new = FALSE )

# Nodes                             | 72sp   | Rod.squirrel
# Mammalia                      	  | t_n72  | t_n268

## 12. ADD THOSE NODES FOR WHICH THERE ARE NO EQUIVALENTS IN 4.7K DATA
##     IN ANY OF THE SUBTREES

# Add manually results for Mammalia
tmp_pos     <- which( names_nodes_2match %in% "t_n73" )
tmp_pos2    <- which( colnames( softbounds_72sp_sum$mean ) %in% "t_n73" )
matching_nodes_d1VSd2[tmp_pos,1] <- matching_nodes_d1VSd2[tmp_pos,4] <- apply( softbounds_72sp_sum$mean, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,2] <- matching_nodes_d1VSd2[tmp_pos,5] <- apply( softbounds_72sp_sum$qdown, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,3] <- matching_nodes_d1VSd2[tmp_pos,6] <- apply( softbounds_72sp_sum$qup, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,7] <- "t_n73|main-72sp"
matching_nodes_d1VSd2[tmp_pos,8] <- 1

# Add manually results for Theria
tmp_pos     <- which( names_nodes_2match %in% "t_n74" )
tmp_pos2    <- which( colnames( softbounds_72sp_sum$mean ) %in% "t_n74" )
matching_nodes_d1VSd2[tmp_pos,1] <- matching_nodes_d1VSd2[tmp_pos,4] <- apply( softbounds_72sp_sum$mean, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,2] <- matching_nodes_d1VSd2[tmp_pos,5] <- apply( softbounds_72sp_sum$qdown, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,3] <- matching_nodes_d1VSd2[tmp_pos,6] <- apply( softbounds_72sp_sum$qup, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,7] <- "t_n73|main-72sp"
matching_nodes_d1VSd2[tmp_pos,8] <- 1

# Add manually results for Placentalia
tmp_pos     <- which( names_nodes_2match %in% "t_n77" )
tmp_pos2    <- which( colnames( softbounds_72sp_sum$mean ) %in% "t_n77" )
matching_nodes_d1VSd2[tmp_pos,1] <- matching_nodes_d1VSd2[tmp_pos,4] <- apply( softbounds_72sp_sum$mean, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,2] <- matching_nodes_d1VSd2[tmp_pos,5] <- apply( softbounds_72sp_sum$qdown, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,3] <- matching_nodes_d1VSd2[tmp_pos,6] <- apply( softbounds_72sp_sum$qup, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,7] <- "t_n77|main-72sp"
matching_nodes_d1VSd2[tmp_pos,8] <- 1

# Add manually results for Atlantogenata
tmp_pos     <- which( names_nodes_2match %in% "t_n78" )
tmp_pos2    <- which( colnames( softbounds_72sp_sum$mean ) %in% "t_n78" )
matching_nodes_d1VSd2[tmp_pos,1] <- matching_nodes_d1VSd2[tmp_pos,4] <- apply( softbounds_72sp_sum$mean, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,2] <- matching_nodes_d1VSd2[tmp_pos,5] <- apply( softbounds_72sp_sum$qdown, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,3] <- matching_nodes_d1VSd2[tmp_pos,6] <- apply( softbounds_72sp_sum$qup, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,7] <- "t_n78|main-72sp"
matching_nodes_d1VSd2[tmp_pos,8] <- 1

# Add manually results for Boreoeutheria
tmp_pos     <- which( names_nodes_2match %in% "t_n82" )
tmp_pos2    <- which( colnames( softbounds_72sp_sum$mean ) %in% "t_n82" )
matching_nodes_d1VSd2[tmp_pos,1] <- matching_nodes_d1VSd2[tmp_pos,4] <- apply( softbounds_72sp_sum$mean, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,2] <- matching_nodes_d1VSd2[tmp_pos,5] <- apply( softbounds_72sp_sum$qdown, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,3] <- matching_nodes_d1VSd2[tmp_pos,6] <- apply( softbounds_72sp_sum$qup, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,7] <- "t_n82|main-72sp"
matching_nodes_d1VSd2[tmp_pos,8] <- 1

# Add manually results for Tupaia-Glires
tmp_pos     <- which( names_nodes_2match %in% "t_n100" )
tmp_pos2    <- which( colnames( softbounds_72sp_sum$mean ) %in% "t_n100" )
matching_nodes_d1VSd2[tmp_pos,1] <- matching_nodes_d1VSd2[tmp_pos,4] <- apply( softbounds_72sp_sum$mean, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,2] <- matching_nodes_d1VSd2[tmp_pos,5] <- apply( softbounds_72sp_sum$qdown, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,3] <- matching_nodes_d1VSd2[tmp_pos,6] <- apply( softbounds_72sp_sum$qup, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,7] <- "t_n100|main-72sp"
matching_nodes_d1VSd2[tmp_pos,8] <- 1

# Add manually results for Glires
tmp_pos     <- which( names_nodes_2match %in% "t_n101" )
tmp_pos2    <- which( colnames( softbounds_72sp_sum$mean ) %in% "t_n101" )
matching_nodes_d1VSd2[tmp_pos,1] <- matching_nodes_d1VSd2[tmp_pos,4] <- apply( softbounds_72sp_sum$mean, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,2] <- matching_nodes_d1VSd2[tmp_pos,5] <- apply( softbounds_72sp_sum$qdown, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,3] <- matching_nodes_d1VSd2[tmp_pos,6] <- apply( softbounds_72sp_sum$qup, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,7] <- "t_n101|main-72sp"
matching_nodes_d1VSd2[tmp_pos,8] <- 1

# Add manually results for Rod-nosquirrel
tmp_pos     <- which( names_nodes_2match %in% "t_n104" )
tmp_pos2    <- which( colnames( softbounds_72sp_sum$mean ) %in% "t_n104" )
matching_nodes_d1VSd2[tmp_pos,1] <- matching_nodes_d1VSd2[tmp_pos,4] <- apply( softbounds_72sp_sum$mean, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,2] <- matching_nodes_d1VSd2[tmp_pos,5] <- apply( softbounds_72sp_sum$qdown, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,3] <- matching_nodes_d1VSd2[tmp_pos,6] <- apply( softbounds_72sp_sum$qup, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,7] <- "t_n104|main-72sp"
matching_nodes_d1VSd2[tmp_pos,8] <- 1

# Add manually results for Rod-nosquirrel
tmp_pos     <- which( names_nodes_2match %in% "t_n104" )
tmp_pos2    <- which( colnames( softbounds_72sp_sum$mean ) %in% "t_n104" )
matching_nodes_d1VSd2[tmp_pos,1] <- matching_nodes_d1VSd2[tmp_pos,4] <- apply( softbounds_72sp_sum$mean, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,2] <- matching_nodes_d1VSd2[tmp_pos,5] <- apply( softbounds_72sp_sum$qdown, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,3] <- matching_nodes_d1VSd2[tmp_pos,6] <- apply( softbounds_72sp_sum$qup, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,7] <- "t_n104|main-72sp"
matching_nodes_d1VSd2[tmp_pos,8] <- 1

# Add manually results for aotus_nancymaae - callithrix_jacchus
tmp_pos     <- which( names_nodes_2match %in% "t_n127" )
tmp_pos2    <- which( colnames( softbounds_72sp_sum$mean ) %in% "t_n127" )
matching_nodes_d1VSd2[tmp_pos,1] <- matching_nodes_d1VSd2[tmp_pos,4] <- apply( softbounds_72sp_sum$mean, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,2] <- matching_nodes_d1VSd2[tmp_pos,5] <- apply( softbounds_72sp_sum$qdown, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,3] <- matching_nodes_d1VSd2[tmp_pos,6] <- apply( softbounds_72sp_sum$qup, 2, mean )[tmp_pos2]
matching_nodes_d1VSd2[tmp_pos,7] <- "t_n127|main-72sp"
matching_nodes_d1VSd2[tmp_pos,8] <- 1

## FINAL STEP: WRITE A TABLE WITH THE SUMMARY !
write.table( x = matching_nodes_d1VSd2, file = "../../03_Extra_analyses/03_Analyses/Matching_table_72sp_4705sp_v0.tsv",
             sep = "\t", quote = FALSE  )


## FINAL STEP: WRITE A TABLE WITH THE SUMMARY !
# Remove row 4 
matching_nodes_d1VSd2_2 <- matching_nodes_d1VSd2[-c(4),]
write.table( x = matching_nodes_d1VSd2_2, file = "../../03_Extra_analyses/03_Analyses/Matching_table_72sp_4705sp.tsv",
             sep = "\t", quote = FALSE  )

## --- RESULTS WITH EXTRA ANALYSES --- ##
# 0. Come back to main wd 
setwd( wd )
wd5 <- gsub( pattern = "03_Analyses..*", replacement = "02_MCMCtree/", x = wd )

# 1. Create lists/vectors needed to find matching nodes
matching_nodes_D1_nuc12mit12 <- matrix( 0, ncol = 12, nrow = length( colnames(softbounds_72sp_sum$mean) ) )
colnames( matching_nodes_D1_nuc12mit12 ) <- c( "mean.data1", "qdown.data1", "qup.data1", 
                                        "mean.nuc12", "qdown.nuc12", "qup.nuc12",
                                        "mean.mit12", "qdown.mit12", "qup.mit12",
                                        "mean.all3p", "qdown.all3p", "qup.all3p" )
names_nodes <- c( "MAMMALIA", "THERIA", "DIDELPHIMORPHIA-AUSTRALIDELPHIA",
                  "EOMETATHERIA", "PLACENTALIA", "XENARTHRA-AFROTHERIA (ATLANTOGENATA)", "XENARTHRA",
                  "AFROTHERIA", "PAENUNGULATA", "BOREOEUTHERIA", "LAURASIATHERIA (Lipotyphla-Carnivora)", 
                  "ERINACEIDAE-SORICIDAE", "SCROTIFERA", "CHIROPTERA", "FEREUNGULATA", "CARNIVORA", 
                  "FELIDAE", "PANTHERINAE", "CANIFORMIA", "ARCTOIDEA", "EUUNGULATA", "ARTIODACTYLA", "ARTIOFABULA",
                  "CETRUMINANTIA", "BOVIDAE", "ovis_aries-capra_hircus", "EUARCHONTOGLIRES", "Tupaia-Glires", 
                  "GLIRES", "LAGOMORPHA", "RODENTIA", "ROD-NOSQUIRREL", "CAVIOMORPHA-PHIOMORPHA", "PHIOMORPHA",
                  "CAVIOMORPHA", "Cavia_porcellus-Cavia_aperea", "Chinchilla_lanigera-Octodon_degus", 
                  "dipodomys_ordi-MYOMORPHA", "MYOMORPHA", "nannospalax_galili-MURIDAE", "MURIDAE", 
                  "CRICETIDAE", "Mesocricetus_auratus-Cricetulus_griseus", "Peromyscus_maniculatus-Microtus_ochrogaster",
                  "MURINAE", "Mus_pahari-rest", "Mus_caroli-rest", "Mus_spretus-Mus_musculus", "PRIMATES", 
                  "STREPSIRHINI", "within_Lemuroidea (propithecus-microcebus)", "HAPLORHINI (Simiiformes-Tarsiiformes)", 
                  "ANTHROPOIDEA/SIMIIFORMES", "AOTIDAE-CALLITRICHIDAE", "aotus_nancymaae - callithrix_jacchus", "CEBIDAE", 
                  "CATARRHINI", "CERCOPITHECOIDEA", "CERCOPITHECINAE", "PAPIONINI", "Papio-Mandrillus",
                  "cercocebus_atys-mandrillus_leucophaeus", "GENUS MACACA", "macaca_fascicularis-macaca_mulatta", "COLOBINAE",
                  "GENUS RHINOPITHECUS", "HOMINOIDEA", "HOMINIDAE", "HOMININAE", "HOMININI", "pan_paniscus-pan_troglodites" )

rownames( matching_nodes_D1_nuc12mit12 ) <- paste( colnames( softbounds_72sp_sum$mean ), names_nodes, sep = " || " )

## FUNCTION 
# This function generates a data.frame with matching nodes 
# between data 1, NUC12CP, and MIT12CP
#
# Arguments:
# inp_mat     Matrix/Data.frame. This object is passed to the function, which 
#             keeps getting updated after having gone through all data subsets
# tmp_nodes   Numeric, temporary nodes for each data subset (2).
# tmp_pos     Numeric, temporary positions to replace in rows in `inp_mat`.
# sum_list    List, list generated for each data subset.
# tmp_nodes72 Numeric, temporary nodes for data 72sp (1).
# sum_list72  List, list generated for data 72sp.
# col_pos     Numeric, vector with the 3 columns to where the info for `sum_list`
#             subtree is to be added.
# ready72     Boolean, TRUE if columns with info from 72s is already added.
#             Otherwise, use FALSE.
sum_tableD1NM <- function( inp_mat, tmp_nodes, tmp_pos, sum_list, tmp_nodes72, sum_list72, subtree,
                           ready72 = TRUE, col_pos = c(4,5,6) ){
  
  # Check if data 1 added
  if ( ready72 == FALSE ){
    cat( "The information for the 72sp has not yet been added, we will\n",
         "add now!\n" )
  }else{
    cat( "The information for the 72sp has already been added, we will\n",
         "not add anything new!\n" )
  }
  
  # Check nodes in data 2
  for ( i in 1:length( tmp_nodes ) ){
    
    cat( "4.7K sp:Evaluating node ", tmp_nodes[i], "...\n" ) 
    tmp_nums <- gsub( pattern = "t_n", replacement = "", x = colnames( sum_list$mean) )
    tmp_ind  <- which( tmp_nums %in% tmp_nodes[i] )
    colnames( sum_list$mean )[tmp_ind] 
    inp_mat[tmp_pos[i],col_pos[1]] <- mean( sum_list$mean[,tmp_ind])
    inp_mat[tmp_pos[i],col_pos[2]] <- mean( sum_list$qdown[,tmp_ind])
    inp_mat[tmp_pos[i],col_pos[3]] <- mean( sum_list$qup[,tmp_ind])
    cat( "Mean: ", mean( sum_list$mean[,tmp_ind]), "\n",
         "Q-2.5%:", mean( sum_list$qdown[,tmp_ind]), "\n",
         "Q-97.5%", mean( sum_list$qup[,tmp_ind]), "\n\n" )
    
    # Data 1
    if ( ready72 == FALSE ){
      cat( "72 sp:Evaluating node ", tmp_nodes72[i], "...\n" ) 
      tmp_nums72 <- gsub( pattern = "t_n", replacement = "", x = colnames( sum_list72$mean) )
      tmp_ind72  <- which( tmp_nums72 %in% tmp_nodes72[i] )
      colnames( sum_list72$mean )[tmp_ind72] 
      inp_mat[tmp_pos[i],1] <- mean( sum_list72$mean[,tmp_ind72])
      inp_mat[tmp_pos[i],2] <- mean( sum_list72$qdown[,tmp_ind72])
      inp_mat[tmp_pos[i],3] <- mean( sum_list72$qup[,tmp_ind72])
      cat( "Mean: ", mean( sum_list72$mean[,tmp_ind72]), "\n",
           "Q-2.5%:", mean( sum_list72$qdown[,tmp_ind72]), "\n",
           "Q-97.5%", mean( sum_list72$qup[,tmp_ind72]), "\n\n" )
    }

  }
  
  # Return object
  return( inp_mat )
  
}

# 2. Data sets 
## 1. Nuc12CP
nuc12cp_f_sum        <- find_prob_MCMC( num_dirs = 16, num_chains = 1, delcol = 2,
                                        name_dir_subt = paste( wd5, "1", sep = "" ),
                                        num_divt = 71,
                                        node_calib = "Calibs_nodes_72sp.csv",
                                        perc = 0.975, clean = TRUE, new = FALSE, new_rev = TRUE )
nuc12cpF_half1 <- apply( X = nuc12cp_f_sum$mean[1:8,], MARGIN = 2, FUN = mean )
nuc12cpF_half2 <- apply( X = nuc12cp_f_sum$mean[9:16,], MARGIN = 2, FUN = mean )

if( ! dir.exists( "plots" ) ){
  dir.create( "plots" )
}

pdf( "plots/01_convergence_plot_nuc12CP.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "72sp - nuc12CP", mean_divt1 = nuc12cpF_half1,
                  mean_divt2 = nuc12cpF_half2, num_runs = 16 )
dev.off() 

# Reuse function to generate new tsv
tmp_nodes72 <- 73:143
tmp_pos     <- which( colnames( softbounds_72sp_sum$mean ) %in% paste( "t_n", tmp_nodes72, sep = "" ) )
matching_nodes_D1_nuc12mit12 <- sum_tableD1NM( inp_mat = matching_nodes_D1_nuc12mit12,
                                               tmp_nodes = 73:143,
                                               tmp_pos = tmp_pos,
                                               sum_list = nuc12cp_f_sum,
                                               tmp_nodes72 = tmp_nodes72,
                                               sum_list72 = softbounds_72sp_sum,
                                               subtree = "NUC-12CP", ready72 = FALSE, col_pos = c(4, 5, 6) )

## 2. Mit12CP
mit12cp_f_sum        <- find_prob_MCMC( num_dirs = 16, num_chains = 1, delcol = 2,
                                        name_dir_subt = paste( wd5, "2", sep = "" ),
                                        num_divt = 71,
                                        node_calib = "Calibs_nodes_72sp.csv",
                                        perc = 0.975, clean = TRUE, new = FALSE, new_rev = TRUE )
mit12cpF_half1 <- apply( X = mit12cp_f_sum$mean[1:8,], MARGIN = 2, FUN = mean )
mit12cpF_half2 <- apply( X = mit12cp_f_sum$mean[9:16,], MARGIN = 2, FUN = mean )

pdf( "plots/02_convergence_plot_mit12CP.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "72sp - mit12CP", mean_divt1 = mit12cpF_half1,
                  mean_divt2 = mit12cpF_half2, num_runs = 16 )
dev.off() 

# Reuse function to generate new tsv
tmp_nodes72 <- 73:143
tmp_pos     <- which( colnames( softbounds_72sp_sum$mean ) %in% paste( "t_n", tmp_nodes72, sep = "" ) )
matching_nodes_D1_nuc12mit12 <- sum_tableD1NM( inp_mat = matching_nodes_D1_nuc12mit12,
                                               tmp_nodes = 73:143,
                                               tmp_pos = tmp_pos,
                                               sum_list = mit12cp_f_sum,
                                               tmp_nodes72 = tmp_nodes72,
                                               sum_list72 = softbounds_72sp_sum,
                                               subtree = "MIT-12CP",
                                               ready72 = TRUE, col_pos = c(7, 8, 9) )

# ## 3. MitRNA -- NOT USED IN PLOTS OR SUMMARY
# mitrna_f_sum        <- find_prob_MCMC( num_dirs = 16, num_chains = 1, delcol = 2,
#                                        name_dir_subt = paste( wd5, "3", sep = "" ),
#                                        num_divt = 71,
#                                        node_calib = "Calibs_nodes_72sp.csv",
#                                        perc = 0.975, clean = TRUE, new = FALSE, new_rev = TRUE )
# mitrnaF_half1 <- apply( X = mitrna_f_sum$mean[1:8,], MARGIN = 2, FUN = mean )
# mitrnaF_half2 <- apply( X = mitrna_f_sum$mean[9:16,], MARGIN = 2, FUN = mean )
# 
# pdf( "plots/03_convergence_plot_mitRNA.pdf", paper = "a4" )
# plot_convergence( name_dir_subt = "72sp - mitRNA", mean_divt1 = mitrnaF_half1,
#                   mean_divt2 = mitrnaF_half2, num_runs = 16 )
# dev.off() 

## 4. All 3 parts
all3parts_f_sum        <- find_prob_MCMC( num_dirs = 16, num_chains = 1, delcol = 6,
                                          name_dir_subt = paste( wd5, "4", sep = "" ),
                                          num_divt = 71,
                                          node_calib = "Calibs_nodes_72sp.csv",
                                          perc = 0.975, clean = TRUE, new = FALSE, new_rev = TRUE )
all3partsF_half1 <- apply( X = all3parts_f_sum$mean[1:8,], MARGIN = 2, FUN = mean )
all3partsF_half2 <- apply( X = all3parts_f_sum$mean[9:16,], MARGIN = 2, FUN = mean )

pdf( "plots/04_convergence_plot_all3parts.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "72sp - all 3 parts", mean_divt1 = all3partsF_half1,
                  mean_divt2 = all3partsF_half2, num_runs = 16 )
dev.off() 

# Reuse function to generate new tsv
matching_nodes_D1_nuc12mit12 <- sum_tableD1NM( inp_mat = matching_nodes_D1_nuc12mit12,
                                               tmp_nodes = 73:143,
                                               tmp_pos = tmp_pos,
                                               sum_list = all3parts_f_sum,
                                               tmp_nodes72 = tmp_nodes72,
                                               sum_list72 = softbounds_72sp_sum,
                                               subtree = "ALL 3 PARTS",
                                               ready72 = TRUE, col_pos = c(10, 11, 12) )

## FINAL STEP: WRITE A TABLE WITH THE SUMMARY !
write.table( x = matching_nodes_D1_nuc12mit12,
             file = "../../03_Extra_analyses/03_Analyses/Matching_table_72sp_NUC12CP_MIT12CP.tsv",
             sep = "\t", quote = FALSE  )

#----------------------------------#
# PLOTS -- COMPARE ESTIMATED TIMES #  
#----------------------------------#
# Get mean times
mean.divt.72sp           <- apply( X = softbounds_72sp_sum$mean*100, MARGIN = 2, FUN = mean )
mean.qup.72sp            <- apply( X = softbounds_72sp_sum$qup*100, MARGIN = 2, FUN = mean )
mean.qdown.72sp          <- apply( X = softbounds_72sp_sum$qdown*100, MARGIN = 2, FUN = mean )

mean.divt.72sp.NUC12CP   <- apply( X = nuc12cp_f_sum$mean*100, MARGIN = 2, FUN = mean )
mean.qup.72sp.NUC12CP    <- apply( X = nuc12cp_f_sum$qup*100, MARGIN = 2, FUN = mean )
mean.qdown.72sp.NUC12CP  <- apply( X = nuc12cp_f_sum$qdown*100, MARGIN = 2, FUN = mean )

mean.divt.72sp.MIT12CP   <- apply( X = mit12cp_f_sum$mean*100, MARGIN = 2, FUN = mean )
mean.qup.72sp.MIT12CP    <- apply( X = mit12cp_f_sum$qup*100, MARGIN = 2, FUN = mean )
mean.qdown.72sp.MIT12CP  <- apply( X = mit12cp_f_sum$qdown*100, MARGIN = 2, FUN = mean )

# mean.divt.72sp.MITRNA    <- apply( X = mitrna_f_sum$mean*100, MARGIN = 2, FUN = mean )
# mean.qup.72sp.MITRNA     <- apply( X = mitrna_f_sum$qup*100, MARGIN = 2, FUN = mean )
# mean.qdown.72sp.MITRNA   <- apply( X = mitrna_f_sum$qdown*100, MARGIN = 2, FUN = mean )

mean.divt.72sp.3PARTS    <- apply( X = all3parts_f_sum$mean*100, MARGIN = 2, FUN = mean )
mean.qup.72sp.3PARTS     <- apply( X = all3parts_f_sum$qup*100, MARGIN = 2, FUN = mean )
mean.qdown.72sp.3PARTS   <- apply( X = all3parts_f_sum$qdown*100, MARGIN = 2, FUN = mean )

# Do the same for those obtained with matching nodes
#[1] "mean.data1"      "qdown.data1"     "qup.data1"       "mean.data2"      "qdown.data1"     "qup.data2"       "data2-nodenames"
# pos_ind <- which( as.numeric( matching_nodes_d1VSd2[,1] ) != 0 )
# mean.72sp.matched    <-  c( as.numeric( matching_nodes_d1VSd2[,1] )[pos_ind]*100 )
# qdown.72sp.matched   <- c( as.numeric( matching_nodes_d1VSd2[,2] )[pos_ind]*100 )
# qup.72sp.matched     <- c( as.numeric( matching_nodes_d1VSd2[,3] )[pos_ind]*100 )
# mean.4705sp.matched  <- c( as.numeric( matching_nodes_d1VSd2[,4])[pos_ind]*100 )
# qdown.4705sp.matched <- c( as.numeric( matching_nodes_d1VSd2[,5] )[pos_ind]*100 )
# qup.4705sp.matched   <- c( as.numeric( matching_nodes_d1VSd2[,6] )[pos_ind]*100 )
## PLOTS WITH "NA"
# mean.72sp.matched    <-  c( as.numeric( matching_nodes_d1VSd2[,1] )*100 )
# qdown.72sp.matched   <- c( as.numeric( matching_nodes_d1VSd2[,2] )*100 )
# qup.72sp.matched     <- c( as.numeric( matching_nodes_d1VSd2[,3] )*100 )
# mean.4705sp.matched  <- c( as.numeric( matching_nodes_d1VSd2[,4])*100 )
# qdown.4705sp.matched <- c( as.numeric( matching_nodes_d1VSd2[,5] )*100 )
# qup.4705sp.matched   <- c( as.numeric( matching_nodes_d1VSd2[,6] )*100 )

mean.72sp.matched    <-  c( as.numeric( matching_nodes_d1VSd2_2[,1] )*100 )
qdown.72sp.matched   <- c( as.numeric( matching_nodes_d1VSd2_2[,2] )*100 )
qup.72sp.matched     <- c( as.numeric( matching_nodes_d1VSd2_2[,3] )*100 )
mean.4705sp.matched  <- c( as.numeric( matching_nodes_d1VSd2_2[,4])*100 )
qdown.4705sp.matched <- c( as.numeric( matching_nodes_d1VSd2_2[,5] )*100 )
qup.4705sp.matched   <- c( as.numeric( matching_nodes_d1VSd2_2[,6] )*100 )

# Plot the following:
# 1. 72-sp VS stitched tree
# 2. 72-sp VS nuclear data
# 3. 72-sp VS mit12CP data
# 4. mit12CP data vs nuc12CP data
setwd( wd )

pdf( "plots/Comparison_divtimes_72spVS4705sp.pdf" )
par( mfrow = c( 2, 2 ) ) 
par( omi = c( .1,.1,.1,.1), mai = c(.8,.8,.4,0.01) )

## PANEL 1 (72sp VS 4705)
pos.red           <- which( matching_nodes_d1VSd2[,8] == 1 )
col.pan1          <- rep( "black", length( mean.72sp.matched ) )
col.pan1[pos.red] <- "red"
plot( 
  x = mean.72sp.matched,
  y = mean.4705sp.matched, 
  # xlab = expression( '72-sp (data 1) - Posterior mean time '*italic(t)*' (Ma)' ),
  # ylab = expression( '4,705-sp (data 2) - Posterior mean time '*italic(t)*' (Ma)' ),
  xlab = expression( 'Posterior on 72s tree - '*italic(t)*' (Ma)' ),
  ylab = expression( 'Posterior on 4.7K tree: all data - '*italic(t)*' (Ma)' ),
  xaxt = 'n', yaxt = 'n',
  ylim = c(0, 2.50)*100,
  xlim = c(0, 2.50)*100,
  pch = 19, cex = 0.8, cex.lab = 1.3, col = col.pan1 )
axis( 1, at = seq(0, 230, by = 20), las = 2, labels = FALSE )
axis( 2, at = seq(0, 230, by = 20), las = 2, cex = 2 )
# Draw VERTICAL arrows
arrows( 
  mean.72sp.matched, qup.4705sp.matched,
  mean.72sp.matched, qdown.4705sp.matched,
  length = 0, lwd = 1.5, col = "darkgrey")
# Draw HORIZONTAL arrows
arrows( 
  qup.72sp.matched, mean.4705sp.matched,
  qdown.72sp.matched, mean.4705sp.matched,
  length = 0, col = "darkgrey", lwd = 1.5 )
# Draw points again to be in front of arrows
points( 
  x = mean.72sp.matched,
  y = mean.4705sp.matched,
  pch = 19, cex = 0.8, col = col.pan1 )
# Add KPg lines
abline( v = 66, col = "lightgrey", lwd = 1, lty = 5 )
abline( h = 66, col = "lightgrey", lwd = 1, lty = 5 )
# Fit lm going through 0
lm_pan1 <- lm( mean.4705sp.matched ~ 0 + mean.72sp.matched )  # 0.9971
cor( mean.4705sp.matched , mean.72sp.matched ) # 0.996605 ~ 0.997
summary( lm_pan1 )
abline( lm( mean.4705sp.matched ~ 0 + mean.72sp.matched ),
        col="black", lty = 1, lwd = 1 )
text( x = 35, y = 140,
      labels = expression( atop( ''*italic( R^2 )*' = 0.9971' ) ),
      cex = 2 )


## PANEL 2 (72sp-tree VS 4.7K NUC12 ONLY (having extracted 72sp present in bb))
plot( 
  x = c( mean.divt.72sp ), 
  y = c( mean.divt.72sp.NUC12CP ), 
  xlab = expression( 'Posterior on 72s tree - '*italic(t)*' (Ma)' ),
  ylab = expression( 'Posterior on 4.7K tree: nuc 12CP only (72s) - '*italic(t)*' (Ma)' ),
  xaxt = 'n', yaxt = 'n',
  ylim = c(0, 2.50)*100,
  xlim = c(0, 2.50)*100,
  pch = 19, cex = 0.8, cex.lab = 1.3 )
axis( 1, at = seq(0, 230, by = 20), las = 2, labels = FALSE )
axis( 2, at = seq(0, 230, by = 20), las = 2, labels = FALSE )
# Draw VERTICAL arrows
arrows( 
  mean.divt.72sp, mean.qup.72sp.NUC12CP,
  mean.divt.72sp, mean.qdown.72sp.NUC12CP,
  length = 0, lwd = 1.5, col = "darkgrey")
# Draw HORIZONTAL arrows
arrows( 
  mean.qup.72sp, mean.divt.72sp.NUC12CP,
  mean.qdown.72sp, mean.divt.72sp.NUC12CP,
  length = 0, col = "darkgrey", lwd = 1.5 )
# Draw points again to be in front of arrows
points( 
  x = c( mean.divt.72sp ),
  y = c( mean.divt.72sp.NUC12CP ),
  pch = 19, cex = 0.8)
# Add KPg lines
abline( v = 66, col = "lightgrey", lwd = 1, lty = 5 )
abline( h = 66, col = "lightgrey", lwd = 1, lty = 5 )
# Fit lm going through 0
lm_pan2 <- lm( mean.divt.72sp.NUC12CP ~ 0 + mean.divt.72sp ) # 0.9937
cor( mean.divt.72sp.NUC12CP , mean.divt.72sp ) # 0.9924919 ~ 0.992
summary( lm_pan2 )
abline( lm( mean.divt.72sp.NUC12CP ~ 0 + c( mean.divt.72sp ) ),
        col="black", lty = 1, lwd = 1 )
text( x = 35, y = 140,
      labels = expression( atop( ''*italic( R^2 )*' = 0.9937' ) ),
      cex = 2 )

## PANEL 3 (72sp-tree VS 4.7K MIT12 ONLY (having extracted 72sp present in bb))
plot( 
  x = c( mean.divt.72sp ), 
  y = c( mean.divt.72sp.MIT12CP ), 
  # xlab = expression( '72-sp 15,268 genes (data 1) - Posterior mean time '*italic(t)*' (Ma)' ),
  # ylab = expression( '4,705-sp 12 mit-12CP genes (data 2) - Posterior mean time '*italic(t)*' (Ma)' ),
  xlab = expression( 'Posterior on 72s tree - '*italic(t)*' (Ma)' ),
  ylab = expression( 'Posterior on 4.7K tree: mit 12CP only (72s) - '*italic(t)*' (Ma)' ),
  xaxt = 'n', yaxt = 'n',
  ylim = c(0, 2.50)*100,
  xlim = c(0, 2.50)*100,
  pch = 19, cex = 0.8, cex.lab = 1.3 )
axis( 1, at = seq(0, 230, by = 20), las = 2, cex = 2 )
axis( 2, at = seq(0, 230, by = 20), las = 2, cex = 2 )
# Draw VERTICAL arrows
arrows( 
  mean.divt.72sp, mean.qup.72sp.MIT12CP,
  mean.divt.72sp, mean.qdown.72sp.MIT12CP,
  length = 0, lwd = 1.5, col = "darkgrey")
# Draw HORIZONTAL arrows
arrows( 
  mean.qup.72sp, mean.divt.72sp.MIT12CP,
  mean.qdown.72sp, mean.divt.72sp.MIT12CP,
  length = 0, col = "darkgrey", lwd = 1.5 )
# Draw points again to be in front of arrows
points( 
  x = c( mean.divt.72sp ),
  y = c( mean.divt.72sp.MIT12CP ),
  pch = 19, cex = 0.8)
# Add KPg lines
abline( v = 66, col = "lightgrey", lwd = 1, lty = 5 )
abline( h = 66, col = "lightgrey", lwd = 1, lty = 5 )
# Fit lm going through 0
lm_pan3 <- lm( mean.divt.72sp.MIT12CP ~ 0 + mean.divt.72sp ) # 0.9959
cor( mean.divt.72sp.MIT12CP , mean.divt.72sp ) # 0.9953834 ~ 0.995
summary( lm_pan3 )
abline( lm( mean.divt.72sp.MIT12CP ~ 0 + c( mean.divt.72sp ) ),
        col="black", lty = 1, lwd = 1 )
text( x = 35, y = 140,
      labels = expression( atop( ''*italic( R^2 )*' = 0.9959 ' ) ),
      cex = 2 )

## PANEL 4 (4.7K NUC12 ONLY VS 4.7K MIT12 ONLY (having extracted 72sp present in bb))
plot( 
  x = c( mean.divt.72sp.NUC12CP ), 
  y = c( mean.divt.72sp.MIT12CP ), 
  xlab = expression( 'Posterior on 4.7K tree: nuc 12CP only (72s) - '*italic(t)*' (Ma)' ),
  ylab = expression( 'Posterior on 4.7K tree: mit 12CP only (72s) - '*italic(t)*' (Ma)' ),
  xaxt = 'n', yaxt = 'n',
  ylim = c(0, 2.50)*100,
  xlim = c(0, 2.50)*100,
  pch = 19, cex = 0.8, cex.lab = 1.3 )
axis( 1, at = seq(0, 230, by = 20), las = 2, labels = FALSE )
axis( 2, at = seq(0, 230, by = 20), las = 2, cex = 2 )
# Draw VERTICAL arrows
arrows( 
  mean.divt.72sp.NUC12CP, mean.qup.72sp.MIT12CP,
  mean.divt.72sp.NUC12CP, mean.qdown.72sp.MIT12CP,
  length = 0, lwd = 1.5, col = "darkgrey")
# Draw HORIZONTAL arrows
arrows( 
  mean.qup.72sp.NUC12CP, mean.divt.72sp.MIT12CP,
  mean.qdown.72sp.NUC12CP, mean.divt.72sp.MIT12CP,
  length = 0, col = "darkgrey", lwd = 1.5 )
# Draw points again to be in front of arrows
points( 
  x = c( mean.divt.72sp.NUC12CP ),
  y = c( mean.divt.72sp.MIT12CP ),
  pch = 19, cex = 0.8)
# Add KPg lines
abline( v = 66, col = "lightgrey", lwd = 1, lty = 5 )
abline( h = 66, col = "lightgrey", lwd = 1, lty = 5 )
# Fit lm going through 0
lm_pan4 <- lm( mean.divt.72sp.MIT12CP ~ 0 + mean.divt.72sp.NUC12CP ) # 0.9938
cor( mean.divt.72sp.MIT12CP , mean.divt.72sp.NUC12CP ) # 0.9921839 ~ 0.992
summary( lm_pan4 )
abline( lm( mean.divt.72sp.MIT12CP ~ 0 + c( mean.divt.72sp.NUC12CP ) ),
        col="black", lty = 1, lwd = 1 )
text( x = 35, y = 140,
      labels = expression( atop( ''*italic( R^2 )*' = 0.9938 ' ) ),
      cex = 2 )

dev.off()


#---------------#
# CALCULATE ESS #
#---------------#
#> ESS with RStan
# Each column is assumed to be an MCMC. Rows are iterations for parameter X
# Source explaining why it is preferable than the function in coda:
# https://nature.berkeley.edu/~pdevalpine/MCMC_comparisons/nimble_MCMC_comparisons.html

# load( "R_ESS_dat.RData" )
ESS_Afro     <- sum_MCMC_ESS( x = Afro_sum$all_mcmc, coda_fun = TRUE )
ESS_Afro$tab
# Tail-ESS Bulk-ESS coda-ESS
# Med.    67274    35355  68760.48
# Min.    13357     7677  15619.30
# Max.   138681   133312 261369.99
min( ESS_Afro$stats$Rhat )
#[1] 0.9999956
dim( Afro_sum$all_mcmc )
#[1] 320016  59
ESS_Xen      <- sum_MCMC_ESS( x = Xen_sum$all_mcmc, coda_fun = TRUE )
ESS_Xen$tab
# Tail-ESS Bulk-ESS  coda-ESS
# Med. 114443.5 104506.5 209990.43
# Min.  74456.0  46625.0  95625.55
# Max. 118293.0 118979.0 238095.00
min( ESS_Xen$stats$Rhat )
#[1] 0.9999918
dim( Xen_sum$all_mcmc )
#[1] 240012 32
ESS_Mar      <- sum_MCMC_ESS( x = Mar_sum$all_mcmc, coda_fun = TRUE )
ESS_Mar$tab
# Tail-ESS Bulk-ESS  coda-ESS
# Med.   3597.5     2674  5598.162
# Min.   1097.0      375  1502.159
# Max.   4515.0     4579 10337.529
min( ESS_Mar$stats$Rhat )
#[1] 0.9997819
dim( Mar_sum$all_mcmc )
#[1] 8852 306
ESS_Eua      <- sum_MCMC_ESS( x = Eua_sum$all_mcmc, coda_fun = TRUE )
ESS_Eua$tab
# Tail-ESS Bulk-ESS  coda-ESS
# Med.      864      777 1621.9549
# Min.      234       92  253.6148
# Max.     1128     1231 2728.8969
min( ESS_Eua$stats$Rhat )
#[1] 0.9990502
dim( Eua_sum$all_mcmc )
#[1] 2097 485
ESS_Lcet     <- sum_MCMC_ESS( x = Lcet_sum$all_mcmc, coda_fun = TRUE )
ESS_Lcet$tab
# Tail-ESS Bulk-ESS  coda-ESS
# Med.     1656     1452 2874.9432
# Min.      679      363  666.9737
# Max.     1968     2109 4669.3114
min( ESS_Lcet$stats$Rhat )
#[1] 0.9994723
dim( Lcet_sum$all_mcmc )
#[1] 3772 430
ESS_Rsq      <- sum_MCMC_ESS( x = Rsq_sum$all_mcmc, coda_fun = TRUE )
ESS_Rsq$tab
# Tail-ESS Bulk-ESS  coda-ESS
# Med.     8233     5751 12054.381
# Min.     2349     2113  3848.437
# Max.    10284    10177 20546.000
min( ESS_Rsq$stats$Rhat )
#[1] 0.9999032
dim( Rsq_sum$all_mcmc )
#[1] 20546 266
ESS_Rcte     <- sum_MCMC_ESS( x = Rcte_sum$all_mcmc, coda_fun = TRUE )
ESS_Rcte$tab
# Tail-ESS Bulk-ESS  coda-ESS
# Med.    14210    11078 22451.953
# Min.     4644     3156  6256.274
# Max.    16316    16444 32813.128
min( ESS_Rcte$stats$Rhat )
#[1] 0.9999387
dim( Rcte_sum$all_mcmc )
#[1] 32571 209
ESS_Lchiro1  <- sum_MCMC_ESS( x = Lchir1_sum$all_mcmc, coda_fun = TRUE )
ESS_Lchiro1$tab
# Tail-ESS Bulk-ESS    coda-ESS
# Med.    93469    35118  75597.1027
# Min.      600      300    313.6009
# Max.   273655   234873 456114.3438
min( ESS_Lchiro1$stats$Rhat )
#[1] 0.9999969
dim( Lchir1_sum$all_mcmc )
#[1] 633815 255
ESS_Lchiro2  <- sum_MCMC_ESS( x = Lchir2_sum$all_mcmc, coda_fun = TRUE )
ESS_Lchiro2$tab
# Tail-ESS Bulk-ESS   coda-ESS
# Med.     2786      729  1907.5397
# Min.      219       67   114.6598
# Max.    13020    10651 20371.2606
min( ESS_Lchiro2$stats$Rhat )
#[1] 0.9999262
dim( Lchir2_sum$all_mcmc )
#[1] 27022 633
ESS_Rsubt1   <- sum_MCMC_ESS( x = Rod1_sum$all_mcmc, coda_fun = TRUE )
ESS_Rsubt1$tab
# Tail-ESS Bulk-ESS  coda-ESS
# Med.     2932      843  2113.362
# Min.      247      125    53.201
# Max.    14387    13662 26973.711
min( ESS_Rsubt1$stats$Rhat )
#[1] 0.9999344
dim( Rod1_sum$all_mcmc )
#[1] 30169 629
ESS_Rsubt2   <- sum_MCMC_ESS( x = Rod2_sum$all_mcmc, coda_fun = TRUE )
ESS_Rsubt2$tab
# Tail-ESS Bulk-ESS   coda-ESS
# Med.     1329      307   866.4240
# Min.       76       43   122.4566
# Max.     8188     7008 15675.5597
min( ESS_Rsubt2$stats$Rhat )
#[1] 0.9998953
dim( Rod2_sum$all_mcmc )
#[1] 17681 690
ESS_Ltherest <- sum_MCMC_ESS( x = Ltherest_sum$all_mcmc, coda_fun = TRUE )
ESS_Ltherest$tab
# Tail-ESS Bulk-ESS    coda-ESS
# Med.   2600.5    848.5  1886.05527
# Min.     66.0     41.0    53.89446
# Max.  10146.0  10237.0 20961.00000
min( ESS_Ltherest$stats$Rhat )
#[1] 0.9999046
dim( Ltherest_sum$all_mcmc )
#[1] 20961 658
ESS_Lag      <- sum_MCMC_ESS( x = Lag_sum$all_mcmc, coda_fun = TRUE )
ESS_Lag$tab
# Tail-ESS Bulk-ESS    coda-ESS
# Med.     3606      881  2806.94167
# Min.       73       23    34.64206
# Max.    18496    12895 26828.21169
min( ESS_Lag$stats$Rhat )
#[1] 0.9999505
dim( Lag_sum$all_mcmc )
#[1] 40002 87

ESS_nuc12CPf      <- sum_MCMC_ESS( x = nuc12cp_f_sum$all_mcmc, coda_fun = TRUE )
ESS_nuc12CPf$tab
# Tail-ESS Bulk-ESS  coda-ESS
# Med.   104512    80607 159081.42
# Min.    28780    13406  33919.13
# Max.   155872   154352 305462.48
min( ESS_nuc12CPf$stats$Rhat )
#[1] 0.9999939
dim( nuc12cp_f_sum$all_mcmc )
#[1] 320016     71

ESS_mit12CPf      <- sum_MCMC_ESS( x = mit12cp_f_sum$all_mcmc, coda_fun = TRUE )
ESS_mit12CPf$tab
# Tail-ESS Bulk-ESS coda-ESS
# Med.   155395   151143 303420.4
# Min.   133775   123654 243254.0
# Max.   159260   159988 320016.0
min( ESS_mit12CPf$stats$Rhat )
#[1] 0.999994
dim( mit12cp_f_sum$all_mcmc )
#[1] 320016     71

ESS_mitRNAf      <- sum_MCMC_ESS( x = mitrna_f_sum$all_mcmc, coda_fun = TRUE )
ESS_mitRNAf
# Tail-ESS Bulk-ESS coda-ESS
# Med.   157586   156993 317029.7
# Min.   147271   127928 253669.0
# Max.   159595   160951 323374.0
min( ESS_mitRNAf$stats$Rhat )
#[1] 0.999994
dim( mitrna_f_sum$all_mcmc )
#[1] 320016     71

ESS_all3partsf      <- sum_MCMC_ESS( x = all3parts_f_sum$all_mcmc, coda_fun = TRUE )
ESS_all3partsf$tab
# Tail-ESS Bulk-ESS  coda-ESS
# Med.   125335   100993 202613.61
# Min.    47743    17870  26419.67
# Max.   156509   144865 290685.92
min( ESS_all3partsf$stats$Rhat )
#[1] 0.9999938
dim( all3parts_f_sum$all_mcmc )
#[1] 320016     71
