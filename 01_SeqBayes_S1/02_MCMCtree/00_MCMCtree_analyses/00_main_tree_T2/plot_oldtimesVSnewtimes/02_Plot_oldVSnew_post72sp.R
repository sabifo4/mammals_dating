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

#------#
# FUNS #
#------#
load_mcmc <- function( mcmc_file ){
  
  # 1. Load files and get parameters
  paste( "Load combined mcmc.txt ... ...\n", sep = "" )
  run    <- read.table( mcmc_file, header = T, sep = "\t" )
  
  # 2. Summarise parameters
  cat( "Generate objects with summarised estimated divergence times... ...\n")
  divtimes <- run[,-c(1, (dim(run)[2]-8):dim(run)[2])]
  
  mean_est_divt  <- apply( X = divtimes, MARGIN = 2, FUN = mean )
  quant_est_divt <- apply( X = divtimes, MARGIN = 2, FUN = quantile, probs = c( 0.025,0.975) )
  quant_est_divt <- t( quant_est_divt )
  all.equal( quant_est_divt[1,], quantile( divtimes[,1], probs = c( 0.025, 0.975 ) ) )
  
  # 3. Return object 
  cat( "\nTasks done! Return objects\n")
  return( list( divt = divtimes, mean_divt = mean_est_divt, quant_divt = quant_est_divt ) )
  
}

load_mcmc_v2 <- function( mcmc1, subt, delcol = 10, perc = 0.975 ){
  
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

#-----------#
# LOAD DATA #
#-----------#
# Load new divtimes 72sp 
post.new.72sp     <- load_mcmc( mcmc_file = "../01_MCMCtree_posterior/02_atlantogenata_tarver2016/mcmc_files/mcmc_tracer.txt" )
divtimes.new72sp  <- post.new.72sp$divt

# Load old divtimes 72sp 
post.old.72sp    <- load_mcmc( mcmc_file = "../01_MCMCtree_posterior_old/mcmc_files_filt/mcmc_tracer.txt" )
divtimes.old72sp <- post.old.72sp$divt

#------------#
# PARSE DATA #
#------------#
# Calib nodes 72sp
calib_nodes <- paste( "t_n", c( 73, 74, 75, 76, 77, 79, 80, 81, 84, 86, 88, 91, 93,
                                94, 96, 97, 99, 101, 102, 103, 104, 110, 117,
                                121, 122, 125, 129, 131, 132, 140, 141, 142 ), sep = "" )
calib_names <- c( "Mammalia", "Theria", "~Marsupialia~", "Eometatheria",
                  "Placentalia", "Xenarthra", "Afrotheria", "Paenungulata",
                  "Lipotyphla", "Chiroptera", "Carnivora", "Caniformia",
                  "Euungulata", "Artiodactyla", "Cetruminantia", "Bovidae",
                  "Euarchontoglires", "Glires", "Lagomorpha", "Rodentia",
                  "Rod-nosquirrel", "Geomyoidea", "Murinae", "Primates",
                  "Strepsirhini", "Anthropoidea", "Catarrhini", "Cercopithecinae",
                  "Papionini", "Hominidae", "Homininae", "Hominini" )
names( calib_nodes ) <- calib_names

#-------#
# PLOTS #
#-------#
plot_72sp_priorVSpost <- function( new_post, old_post, cal_name, cal_node,
                                   info_legend = c( "New-post72sp", "Old-post72sp" ),
                                   col_legend = c( "black", "red" ) ) {
  
  # 1. Find limit axis
  max_x_newpost  <- round( max( density( new_post )$x ) + 0.5 )
  max_x_oldpost  <- round( max( density( old_post )$x) + 0.5 )
  max_y_newpost  <- round( max( density( new_post )$y ) + 0.5 )
  max_y_oldpost  <- round( max( density( old_post )$y) + 0.5 )
  x_max <- max( max_x_newpost, max_x_oldpost )
  y_lim <- max( max_y_newpost, max_y_oldpost )
  
  # 2. Plot
  ## New posterior
  d   <- density( new_post )#, adj = 1 )
  #d$y <- d$y/max(d$y)
  plot( d,
        xlim = c( 0, x_max ), ylim = c( 0, y_lim ), 
        main=NULL, xlab = '', ylab = '', cex.axis = 0.7, mgp = c(2.5,0,0),
        #xaxt = "n", yaxt = "n", 
        tck = 0.01, col.main = "white" )
  title( main = cal_name, sub = cal_node,
         line = -4, adj = 0.01,
         cex.main = 0.8,
         cex.sub = 0.7 )
  info.legend <- info_legend
  col.legend  <- col_legend
  legend( x = x_max-0.2, y = y_lim/2, legend = info.legend, col = col.legend,
          lty = 1, bty = 'n', cex = 0.7 )
  ## PRIOR-72sp       
  d2   <- density( old_post )#, adj = 1 )
  #d2$y <- d2$y/max(d2$y)
  lines( d2, col = "red" )
  
}

## PLOT ##
# Plot old posterior VS new posterior
all.equal( colnames( divtimes.new72sp ), colnames( divtimes.old72sp ))
# TRUE
pdf( file = "00_Check_oldpostVSnewpost-I.pdf" )
par( mfrow = c(8,2), mai=c(0.1,0.1,0.1,0.1) )
for( i in 1:16 ){
  ind  <- which( colnames( divtimes.new72sp ) == calib_nodes[i] )
  plot_72sp_priorVSpost( new_post = divtimes.new72sp[,ind],
                         old_post = divtimes.old72sp[,ind],
                         cal_name = names( calib_nodes )[i],
                         cal_node = calib_nodes[i],
                         info_legend = c( "New-post72sp", "Old-post72sp" ),
                         col_legend = c( "black", "red" ) )
}
dev.off()

pdf( file = "00_Check_oldpostVSnewpost-II.pdf" )
par( mfrow = c(8,2), mai=c(0.1,0.1,0.1,0.1) )
for( i in 17:32 ){
  ind  <- which( colnames( divtimes.new72sp ) == calib_nodes[i] )
  plot_72sp_priorVSpost( new_post = divtimes.new72sp[,ind],
                         old_post = divtimes.old72sp[,ind],
                         cal_name = names( calib_nodes )[i],
                         cal_node = calib_nodes[i],
                         info_legend = c( "New-post72sp", "Old-post72sp" ),
                         col_legend = c( "black", "red" ))
}
dev.off()

# The distributions are almost identical. Removing the 11 nuclear genes that 
# are common in the second data set is needed for the sequential Bayesian approach 
# but there is no impact in the posterior divergence times obtained with the 
# data set that included them.