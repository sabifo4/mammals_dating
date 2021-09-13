#-------------------#
# Clean environment #
#-------------------#
rm( list = ls( ) )

#-----------------------#
# Set working directory #
#-----------------------#
library( rstudioapi ) 
# Get the path to current open R script and find main dir
path_to_file <- getActiveDocumentContext()$path
wd <- paste( dirname( path_to_file ), "/", sep = "" )
setwd( wd )

#-----------#
# FUNCTIONS #
#-----------#
# FUNCTION 1
# This function outputs a list with two entries: one with the 
# CI width for each node and another one with the estimated 
# divergence times.
#
# Arguments
#    X   data.frame, object with divtimes
sum_divt <- function( X ){
  
  # Get mean and quantiles and format it
  mean_est_divt  <- apply( X = X, MARGIN = 2, FUN = mean )
  quant_est_divt <- t( apply( X = X, MARGIN = 2, FUN = quantile, probs = c(0.025,0.975) ) )
  est_mat        <- cbind( mean_est_divt, quant_est_divt )
  colnames( est_mat )   <- c( "Posterior.mean", "95%CI-low", "95%CI-up" )
  CIwidth               <- as.matrix( est_mat[,3] - est_mat[,2],
                                      ncol = 1, nrow = dim(est_mat)[1] )
  colnames( CIwidth ) <- "CI-width"
  
  # Return formatted object
  return( list( CI = CIwidth, divt = est_mat ) )
}

#-----------#
# LOAD DATA #
#-----------#
# 1. Create global vars
dat <- c( "sample1sp", "sample10sp", "sample30sp", 
          "sample50sp", "sample100sp",
          "sample500sp", "sample1000sp", 
          
          "sample10sp_all", "sample30sp_all", "sample100sp_all",
          "sample500sp_all", "sample1000sp_all",
          
          "ALL" )
dat.mcmc.1.GBM <- dat.mcmc.2.GBM <- dat.mcmc.3.GBM <- 
  dat.mcmc.1.ILN <- dat.mcmc.2.ILN <- dat.mcmc.3.ILN <- 
  divtimes.1.GBM <- divtimes.2.GBM <- divtimes.3.GBM <- 
  divtimes.1.ILN <- divtimes.2.ILN <- divtimes.3.ILN <- vector( mode = "list", length( dat ) )
names( dat.mcmc.1.GBM ) <- names( dat.mcmc.2.GBM ) <- names( dat.mcmc.3.GBM ) <-
  names( dat.mcmc.1.ILN ) <- names( dat.mcmc.2.ILN ) <- names( dat.mcmc.3.ILN ) <- 
  names( divtimes.1.GBM ) <- names( divtimes.2.GBM ) <- names( divtimes.3.GBM ) <- 
  names( divtimes.1.ILN ) <- names( divtimes.2.ILN ) <- names( divtimes.3.ILN ) <- dat 

# 2. Get data
for( i in 1:length( dat ) ){

  if( dat[i] == "sample1sp" ){
    cat( "Parsing 1 sample! Appending the 1g-conc results to all lists\n\n" )
    dat.mcmc.1.GBM[[ i ]] <- dat.mcmc.2.GBM[[ i ]] <- dat.mcmc.3.GBM[[ i ]] <- 
      read.table( file = paste( "../", dat[i], "/1/mcmc_files_GBM/mcmc_tracer.txt",
                                sep = "" ), stringsAsFactors = FALSE, header = TRUE, sep = "\t" )
    dat.mcmc.1.ILN[[ i ]] <- dat.mcmc.2.ILN[[ i ]] <- dat.mcmc.3.ILN[[ i ]] <- 
      read.table( file = paste( "../", dat[i], "/1/mcmc_files_ILN/mcmc_tracer.txt",
                                sep = "" ),
                  stringsAsFactors = FALSE, header = TRUE, sep = "\t" )
  }else{
    if( i == length( dat ) ){
      cat( "Parsing all data! Appending the 16K-4parts results to all lists\n\n" )
      dat.mcmc.3.GBM[[ i ]] <- dat.mcmc.2.GBM[[ i ]] <- dat.mcmc.1.GBM[[ i ]] <- 
        read.table( file = "../../02_atlantogenata_tarver2016/mcmc_files_GBM/mcmc_tracer.txt",
                    stringsAsFactors = FALSE, header = TRUE, sep = "\t" )
      dat.mcmc.3.ILN[[ i ]] <- dat.mcmc.2.ILN[[ i ]] <- dat.mcmc.1.ILN[[ i ]] <- 
        read.table( file = "../../02_atlantogenata_tarver2016/mcmc_files_ILN/mcmc_tracer.txt",
                    stringsAsFactors = FALSE, header = TRUE, sep = "\t" )
    }else{
      cat( "Parsing GBM: ", dat[i], "- concatenated ... \n" )
      dat.mcmc.1.GBM[[ i ]] <- read.table( file = paste( "../", dat[i], "/1/mcmc_files_GBM/mcmc_tracer.txt",
                                                         sep = "" ),
                                           stringsAsFactors = FALSE, header = TRUE, sep = "\t" )
      cat( "Parsing GBM: ", dat[i], "- 2 partitions ... \n" )
      dat.mcmc.2.GBM[[ i ]] <- read.table( file = paste( "../", dat[i], "/2/mcmc_files_GBM/mcmc_tracer.txt",
                                                         sep = "" ),
                                           stringsAsFactors = FALSE, header = TRUE, sep = "\t" )
      cat( "Parsing GBM: ", dat[i], "- 4 partitions ... \n" )
      dat.mcmc.3.GBM[[ i ]] <- read.table( file = paste( "../", dat[i], "/3/mcmc_files_GBM/mcmc_tracer.txt",
                                                         sep = "" ),
                                           stringsAsFactors = FALSE, header = TRUE, sep = "\t" )
      
      cat( "Parsing ILN: ", dat[i], "- concatenated ... \n" )
      dat.mcmc.1.ILN[[ i ]] <- read.table( file = paste( "../", dat[i], "/1/mcmc_files_ILN/mcmc_tracer.txt",
                                                         sep = "" ),
                                           stringsAsFactors = FALSE, header = TRUE, sep = "\t" )
      cat( "Parsing ILN: ", dat[i], "- 2 partitions ... \n" )
      dat.mcmc.2.ILN[[ i ]] <- read.table( file = paste( "../", dat[i], "/2/mcmc_files_ILN/mcmc_tracer.txt",
                                                         sep = "" ),
                                           stringsAsFactors = FALSE, header = TRUE, sep = "\t" )
      cat( "Parsing ILN: ", dat[i], "- 4 partitions ... \n\n" )
      dat.mcmc.3.ILN[[ i ]] <- read.table( file = paste( "../", dat[i], "/3/mcmc_files_ILN/mcmc_tracer.txt",
                                                         sep = "" ),
                                           stringsAsFactors = FALSE, header = TRUE, sep = "\t" )
    }
  }
  
  
}

# 3. Parse data using `sum_divt` so we can then use the code to plot 
#    the infinite-sites plots below
## NOTE: We added 16K-parts to all "3" entries and the same with 1gene-conc
## so there are no issues with the lists.
for ( i in 1:length( dat ) ){
  divtimes.1.GBM[[ i ]] <- sum_divt( X = dat.mcmc.1.GBM[[ i ]][,-c(1,(dim(dat.mcmc.1.GBM[[ i ]])[2]-2):(dim(dat.mcmc.1.GBM[[ i ]])[2]))] )
  divtimes.2.GBM[[ i ]] <- sum_divt( X = dat.mcmc.2.GBM[[ i ]][,-c(1,(dim(dat.mcmc.2.GBM[[ i ]])[2]-4):(dim(dat.mcmc.2.GBM[[ i ]])[2]))] )
  divtimes.3.GBM[[ i ]] <- sum_divt( X = dat.mcmc.3.GBM[[ i ]][,-c(1,(dim(dat.mcmc.3.GBM[[ i ]])[2]-8):(dim(dat.mcmc.3.GBM[[ i ]])[2]))] )
  
  divtimes.1.ILN[[ i ]] <- sum_divt( X = dat.mcmc.1.ILN[[ i ]][,-c(1,(dim(dat.mcmc.1.ILN[[ i ]])[2]-2):(dim(dat.mcmc.1.ILN[[ i ]])[2]))] )
  divtimes.2.ILN[[ i ]] <- sum_divt( X = dat.mcmc.2.ILN[[ i ]][,-c(1,(dim(dat.mcmc.2.ILN[[ i ]])[2]-4):(dim(dat.mcmc.2.ILN[[ i ]])[2]))] )
  divtimes.3.ILN[[ i ]] <- sum_divt( X = dat.mcmc.3.ILN[[ i ]][,-c(1,(dim(dat.mcmc.3.ILN[[ i ]])[2]-8):(dim(dat.mcmc.3.ILN[[ i ]])[2]))] )
}

# Generate global vars
subsample.divt.1.GBM <- subsample.divt.2.GBM <- subsample.divt.3.GBM <- 
  subsample.divt.1.ILN <- subsample.divt.2.ILN <- subsample.divt.3.ILN <- vector( mode = "list", 7 )
names( subsample.divt.1.GBM ) <- names( subsample.divt.2.GBM ) <- names( subsample.divt.3.GBM ) <-
  names( subsample.divt.1.ILN ) <- names( subsample.divt.2.ILN ) <- names( subsample.divt.3.ILN ) <- dat[c(1,8:13)] 
count <- 0

## NOTE: We are only using `sample1sp`, `sample10sp_all`, `sample30sp_all`,
## `sample100sp_all`, `sample500sp_all`, `sample1000sp_all`, `ALL`. Those are in 
## entries 1 and from 8 to 13 in `divtimes` objects.
for( i in c(1,8:13) ){
  count <- count + 1
  subsample.divt.1.GBM[[ count ]] <- divtimes.1.GBM[[ i ]]
  subsample.divt.2.GBM[[ count ]] <- divtimes.2.GBM[[ i ]]
  subsample.divt.3.GBM[[ count ]] <- divtimes.3.GBM[[ i ]]
  subsample.divt.1.ILN[[ count ]] <- divtimes.1.ILN[[ i ]]
  subsample.divt.2.ILN[[ count ]] <- divtimes.2.ILN[[ i ]]
  subsample.divt.3.ILN[[ count ]] <- divtimes.3.ILN[[ i ]]
} 
  
# Generate global vars
subsample2.divt.1.GBM <- subsample2.divt.2.GBM <- subsample2.divt.3.GBM <- 
  subsample2.divt.1.ILN <- subsample2.divt.2.ILN <- subsample2.divt.3.ILN <- vector( mode = "list", 7 )
names( subsample2.divt.1.GBM ) <- names( subsample2.divt.2.GBM ) <- names( subsample2.divt.3.GBM ) <-
  names( subsample2.divt.1.ILN ) <- names( subsample2.divt.2.ILN ) <- names( subsample2.divt.3.ILN ) <- dat[c(1:3,5:7,13)] 
count <- 0

## NOTE: We are now just getting the data from `sample1sp`, `sample10sp`, `sample30sp`,
## `sample100sp`, `sample500sp`, `sample1000sp`, `ALL`. Those are in 
## entries 1 to 3, 5 to 7, and 13 in `divtimes` objects. This is because we want to 
## compare the effect of sampling from the whole pool of filtered genes (step above)
## and from only those genes present in all 72 mammal taxa (for 500sp and 1000sp objects,
## genes from the big pool were randomly sampled and added to the 645 genes, see R script 
## `01_Analysis_filtered_genes.R` for details).
for( i in c(1:3,5:7,13) ){
  count <- count + 1
  subsample2.divt.1.GBM[[ count ]] <- divtimes.1.GBM[[ i ]]
  subsample2.divt.2.GBM[[ count ]] <- divtimes.2.GBM[[ i ]]
  subsample2.divt.3.GBM[[ count ]] <- divtimes.3.GBM[[ i ]]
  subsample2.divt.1.ILN[[ count ]] <- divtimes.1.ILN[[ i ]]
  subsample2.divt.2.ILN[[ count ]] <- divtimes.2.ILN[[ i ]]
  subsample2.divt.3.ILN[[ count ]] <- divtimes.3.ILN[[ i ]]
} 

#----------------------------#
# No. loci (x) VS post.times #
#----------------------------#
# 73: Mammalia
# 77: Placentalia
# 79: Xenarthra
# 80: Afrotheria
# 86: Chiroptera
# 88: Carnivora 
# 101: Lagomorpha
# 121: Primates 
# 102: Rodentia 
nums     <- c( 73, 77, 79, 80, 86, 88, 101, 121, 102 )
ind.rows <- which( rownames( divtimes.1.GBM$sample30sp$divt ) %in% paste( "t_n", nums, sep = "" ) )

# Set vars for plots 
min_lims <- c( 1.6, 0.7, 0.4, 0.5, 0.4, 0.4, 0.6, 0.4, 0.5 )*100
max_lims <- c( 2.6, 0.9, 0.8, 0.8, 0.6, 0.6, 0.7, 0.6, 0.7 )*100
main_titles <- c( "Mammalia", "Placentalia", "Xenarthra", "Afrotheria",
                  "Chiroptera", "Carnivora", "Lagomorpha", "Primates",
                  "Rodentia" )
chosen <- c(1, 2, 3, 4, 8, 9 )


#-------#
# PLOTS #
#-------#
# FUNCTION
# Function to plot the number of loci randomly sampled VS the posterior mean time 
#
# Arguments:
#
# x_coords   Numeric, vector with the values in the x axis where the 
#            the data for each data set will be plotted.
#            Suggested: c( 10, 20, 30, 40, 50, 60, 70 ).
# dat        List. There are as many entries as data sets to plot.
#            Now, we have 7 data sets, so 7 entries.
# len        Numeric, number of data sets in the list passed to `dat`.
# min_ylim   Numeric, vector with the minimum values of the x axis 
#            for each data set to be plotted.
# max_ylim   Numeric, vector with the minimum values of the y axis 
#            for each data set to be plotted.
# ind.row    Numeric, vector with the position of the nodes for which 
#            estimated times we want to plot.
# cols       Character, vector with the colours to use to plot 
#            the results for each selected node.
# labs       Character, vector with the labels to be used in the 
#            x axis.
plot_lenVSpost2 <- function( x_coords, dat, len, min_ylim, max_ylim,
                             ind.row, cols, labs ){
  
  # Set global vars for divtimes
  # List: one entry per node and each entry has a matrix with the node values
  # per sample (mean, CIdown, CIup)
  tmp.dat    <- vector( mode = "list", length = length( ind.row ) )
  tmp        <- matrix( 0, nrow = len, ncol = 3 )
  node.names <- rownames( dat[[ 1 ]]$divt )[ind.row]
  names( tmp.dat ) <- node.names
  for( j in 1:length( tmp.dat ) ){
    for( i in 1:len ){
      tmp[i,1] <- dat[[ i ]]$divt[ind.row[j],1]*100
      tmp[i,2] <- dat[[ i ]]$divt[ind.row[j],2]*100
      tmp[i,3] <- dat[[ i ]]$divt[ind.row[j],3]*100
      tmp.dat[[ j ]] <- tmp
    }
  }
  
  # Iterate over the amount of samples
  for( i in 1:length( tmp.dat ) ){
    if( i == 1 ){
      plot( x = x_coords,
            y = tmp.dat[[ i ]][,1],
            xlab = c( "Number of loci" ),
            ylab = expression( 'Posterior mean time '*italic(t)*' (Ma)' ),
            cex.lab = 1.5, xaxt = 'n',
            # xlab = '', 
            #ylab = '', #yaxt = 'n', 
            ylim = c(min_ylim,max_ylim),
            pch = 19, cex = 0.8,
            col = cols[i] )
    }else{
      points( x = x_coords, y = tmp.dat[[ i ]][,1], pch = 19, col = cols[i] )
    }
    
    lines(  x = x_coords,
            y = tmp.dat[[ i ]][,2],
            lty = 2, lwd = 0.2, col = adjustcolor( cols[i], alpha.f = 0.2 )
    )
    lines(  x = x_coords,
            y = tmp.dat[[ i ]][,3],
            lty = 2, lwd = 0.2, col = adjustcolor( cols[i], alpha.f = 0.2 )
    )
  }
  
  title( main = "All together" ) 
  axis( side = 1, labels = labs, at = x_coords )
  
  
  # Return object used to plot
  return( tmp.dat )
}

# PLOT 1 (genes present in all 645genes + added from big pool if needed)
dat_plotting.conc.GBM <- plot_lenVSpost2( x_coords = c( 10, 20, 30, 40, 50, 60, 70 ),
                                          dat = subsample2.divt.1.GBM,
                                          len = 7,
                                          min_ylim = min( min_lims[chosen] ), 
                                          max_ylim = max( max_lims[chosen] ),
                                          ind.row = ind.rows[chosen],
                                          cols = c( "black", "red", "green",
                                                    "blue", "orange", "purple"),
                                          labs = c( "1", "10", "30", "100", 
                                                    "500", "1,000", "~16,000" )
)

dat_plotting.2parts.GBM <- plot_lenVSpost2( x_coords = c( 10, 20, 30, 40, 50, 60, 70 ),
                                            dat = subsample2.divt.2.GBM,
                                            len = 7,
                                            min_ylim = min( min_lims[chosen] ), 
                                            max_ylim = max( max_lims[chosen] ),
                                            ind.row = ind.rows[chosen],
                                            cols = c( "black", "red", "green",
                                                      "blue", "orange", "purple"),
                                            labs = c( "1", "10", "30", "100", 
                                                      "500", "1,000", "~16,000" )
)

dat_plotting.4parts.GBM <- plot_lenVSpost2( x_coords = c( 10, 20, 30, 40, 50, 60, 70 ),
                                            dat = subsample2.divt.3.GBM,
                                            len = 7,
                                            min_ylim = min( min_lims[chosen] ), 
                                            max_ylim = max( max_lims[chosen] ),
                                            ind.row = ind.rows[chosen],
                                            cols = c( "black", "red", "green",
                                                      "blue", "orange", "purple"),
                                            labs = c( "1", "10", "30", "100", 
                                                      "500", "1,000", "~16,000" )
)


# PLOT 2 (all genes randomly sampled from the big pool of filtered genes)
dat_plotting2.conc.GBM <- plot_lenVSpost2( x_coords = c( 10, 20, 30, 40, 50, 60, 70 ),
                                          dat = subsample.divt.1.GBM,
                                          len = 7,
                                          min_ylim = min( min_lims[chosen] ), 
                                          max_ylim = max( max_lims[chosen] ),
                                          ind.row = ind.rows[chosen],
                                          cols = c( "black", "red", "green",
                                                    "blue", "orange", "purple"),
                                          labs = c( "1", "10", "30", "100", 
                                                    "500", "1,000", "~16,000" )
)

dat_plotting2.2parts.GBM <- plot_lenVSpost2( x_coords = c( 10, 20, 30, 40, 50, 60, 70),
                                            dat = subsample.divt.2.GBM,
                                            len = 7,
                                            min_ylim = min( min_lims[chosen] ), 
                                            max_ylim = max( max_lims[chosen] ),
                                            ind.row = ind.rows[chosen],
                                            cols = c( "black", "red", "green",
                                                      "blue", "orange", "purple"),
                                            labs = c( "1", "10", "30", "100", 
                                                      "500", "1,000", "~16,000" )
)

dat_plotting2.4parts.GBM <- plot_lenVSpost2( x_coords = c( 10, 20, 30, 40, 50, 60, 70 ),
                                            dat = subsample.divt.3.GBM,
                                            len = 7,
                                            min_ylim = min( min_lims[chosen] ), 
                                            max_ylim = max( max_lims[chosen] ),
                                            ind.row = ind.rows[chosen],
                                            cols = c( "black", "red", "green",
                                                      "blue", "orange", "purple"),
                                            labs = c( "1", "10", "30", "100", 
                                                      "500", "1,000", "~16,000" )
)


