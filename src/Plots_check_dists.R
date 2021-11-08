# -------------------------------------------------- # 
#   FUNCTIONS USED WHEN CHECKING POST DISTS AFTER    #
#                  MCMC SAMPLING                     #
# -------------------------------------------------- # 
# The two functions below are used to check if there #
# is not any conflict between node calibrations when #
# sampling from the prior in MCMCtree                #
#                                                    #
# -------------------------------------------------- # 

#=======================================================================#
# FUNCTION TO DEFINE ST CALIBRATIONS -- VSHARPMIN
#-----------------------------------------------------------------------#
# Arguments
#     log_ST   character, path to the log output file printed out 
#              by MCMCtree during the run 
#     ST.fitted.dists.RData  character, path to the RData file output 
#                            when fitting ST distributions to posterior 
#                            time estimates when using 72sp 
#     ST.calib.node.names    character, vector of length as many ST 
#                            calibrations, previously created in the 
#                            script, with the names of each calibrated node
#
define_ST_vsharpmin <- function( log_ST, ST.fitted.dists.RData, ST.calib.nodes.names ){
  
  # 1. Read out.txt of one of the runs and extract info pasted above
  cat( "Loading log screen file saved when running MCMCtree... ...\n")
  outMCMCtree <- readLines( con = log_ST )
  # 2. Extract lines with the ST calibrations
  cat( "Generating objects with ST calibrations... ...\n")
  ST.text <- grep( pattern = "Node ..* ST ..*", x = outMCMCtree, value = TRUE )
  # 3. Extract the position of the calibrated nodes
  ST.calib.nodes <- as.numeric( unlist( stringr::str_split( string = gsub( pattern = 
                                                                             "Node  |Node |: ..*", 
                                                                           replacement = "",
                                                                           x = ST.text ),
                                                            pattern = " " ) ) )
  # 4. Extract the parameters for the ST distributions
  ST.calibs <- lapply( X = stringr::str_split( gsub( pattern = "Node.*ST.*\\( | \\)|,",
                                                     replacement = "",
                                                     x = ST.text ),
                                               pattern = " " ),
                       FUN = as.numeric )
  names( ST.calibs ) <- paste( "t_n", ST.calib.nodes, sep = "" )
  # 5. Name each node accordingly                          
  names( ST.calib.nodes ) <- ST.calib.nodes.names
  # 6. Add ST calibrations used for these nodes 
  for( i in 1:length( ST.calib.nodes ) ){
    
    names( ST.calibs[[ i ]] )[1] <- "xi"    #position ~ approx. to mean age
    names( ST.calibs[[ i ]] )[2] <- "omega" #scale ~ variance (small number narrow, large wider)
    names( ST.calibs[[ i ]] )[3] <- "alpha" #shape (slanding, going left or right)
    names( ST.calibs[[ i ]] )[4] <- "nu"    #degrees of freedom (tail)
    
  }
  # 7. Load RData object with ST calibrations fitted to the post calibrations of 
  #    72sp mammals to find the equivalent node position
  cat( "Load RData object with ST calibrations fitted to posterior time estimates
       when using the 72sp mammals phylogeny... ...\n")
  # 210317: Now the object is "ST.fitted.vsharpmin.dists"
  load( ST.fitted.dists.RData )
  tmp <- FALSE
  nodes.72sp <- rep( "", length( ST.calibs ) )
  cat( "Finding matching node positions with nodes calibrated with ST distributions... ...\n")
  for( i in 1:length( ST.calibs ) ){
    
    for( j in 1:length( ST.fitted.vsharpmin.dists ) ){
      tmp.ndec <- nchar( gsub( pattern = "..*\\.", x = ST.calibs[[i]], replacement = "" ) )
      if( any( tmp.ndec == 4 ) ){
        tmp <- all.equal( round( ST.fitted.vsharpmin.dists[[ j ]], 4 ), round( ST.calibs[[ i ]], 4 ) )
      }else{
        tmp <- all.equal( round( ST.fitted.vsharpmin.dists[[ j ]], 3 ), round( ST.calibs[[ i ]], 3 ) )
      }
      
      if( tmp == TRUE ){
        break
      }
    }
    
    # Add index label and set tmp to FALSE
    nodes.72sp[ i ] <- names( ST.fitted.vsharpmin.dists )[j]
    tmp <- FALSE
    
  }
  # 8. Delete ST.fitted.dists obj
  remove( ST.fitted.vsharpmin.dists )
  
  # 9. Return list with objects
  cat( "\nTasks done! Return objects with ST info\n")
  return( list( ST.calibs = ST.calibs, nodes.72sp = nodes.72sp, ST.calib.nodes = ST.calib.nodes ))
  
}



#=======================================================================#
# FUNCTION TO DEFINE ST CALIBRATIONS 
#-----------------------------------------------------------------------#
# Arguments
#     log_ST   character, path to the log output file printed out 
#              by MCMCtree during the run 
#     ST.fitted.dists.RData  character, path to the RData file output 
#                            when fitting ST distributions to posterior 
#                            time estimates when using 72sp 
#     ST.calib.node.names    character, vector of length as many ST 
#                            calibrations, previously created in the 
#                            script, with the names of each calibrated node
#
define_ST <- function( log_ST, ST.fitted.dists.RData, ST.calib.nodes.names,
                       SN = FALSE, SN.calibs.nodes.names = FALSE ){
  
  # 1. Read out.txt of one of the runs and extract info pasted above
  cat( "Loading log screen file saved when running MCMCtree... ...\n")
  outMCMCtree <- readLines( con = log_ST )
  # 2. Extract lines with the ST calibrations
  cat( "Generating objects with ST calibrations... ...\n")
  ST.text <- grep( pattern = "Node ..* ST ..*", x = outMCMCtree, value = TRUE )
  # 3. Extract the position of the calirbated nodes
  ST.calib.nodes <- as.numeric( unlist( stringr::str_split( string = gsub( pattern = 
                                                                             "Node  |Node |: ..*", 
                                                                           replacement = "",
                                                                           x = ST.text ),
                                                            pattern = " " ) ) )
  # 4. Extract the parameters for the ST distributions
  ST.calibs <- lapply( X = stringr::str_split( gsub( pattern = "Node.*ST.*\\( | \\)|,",
                                                     replacement = "",
                                                     x = ST.text ),
                                               pattern = " " ),
                       FUN = as.numeric )
  names( ST.calibs ) <- paste( "t_n", ST.calib.nodes, sep = "" )
  # 5. Name each node accordingly                          
  names( ST.calib.nodes ) <- ST.calib.nodes.names
  # 6. Add ST calibrations used for these nodes 
  for( i in 1:length( ST.calib.nodes ) ){
    
    names( ST.calibs[[ i ]] )[1] <- "xi"    #position ~ approx. to mean age
    names( ST.calibs[[ i ]] )[2] <- "omega" #scale ~ variance (small number narrow, large wider)
    names( ST.calibs[[ i ]] )[3] <- "alpha" #shape (slanding, going left or right)
    names( ST.calibs[[ i ]] )[4] <- "nu"    #degrees of freedom (tail)
    
  }
  # 7. Load RData object with ST calibrations fitted to the post calibrations of 
  #    72sp mammals to find the equivalent node position
  cat( "Load RData object with ST calibrations fitted to posterior time estimates
       when using the 72sp mammals phylogeny... ...\n")
  load( ST.fitted.dists.RData )
  tmp <- FALSE
  nodes.72sp <- rep( "", length( ST.calibs ) )
  cat( "Finding matching node positions with nodes calibrated with ST distributions... ...\n")
  for( i in 1:length( ST.calibs ) ){
    
    for( j in 1:length( ST.fitted.dists ) ){
      tmp.ndec <- nchar( gsub( pattern = "..*\\.", x = ST.calibs[[i]], replacement = "" ) )
      if( any( tmp.ndec == 4 ) ){
        tmp <- all.equal( round( ST.fitted.dists[[ j ]], 4 ), round( ST.calibs[[ i ]], 4 ) )
      }else{
        tmp <- all.equal( round( ST.fitted.dists[[ j ]], 3 ), round( ST.calibs[[ i ]], 3 ) )
      }
      
      if( tmp == TRUE ){
        break
      }
    }
    
    # Add index label and set tmp to FALSE
    nodes.72sp[ i ] <- names( ST.fitted.dists )[j]
    tmp <- FALSE
    
  }
  # 8. Delete ST.fitted.dists obj
  remove( ST.fitted.dists )
  
  # Check if SN available !
  if ( SN == TRUE ){
    # 1. Extract lines with the SN calibrations
    cat( "Generating objects with SN calibrations... ...\n")
    SN.text <- grep( pattern = "Node ..* SN ..*", x = outMCMCtree, value = TRUE )
    # 3. Extract the position of the calibrated nodes
    SN.calib.nodes <- as.numeric( unlist( stringr::str_split( string = gsub( pattern = 
                                                                               "Node  |Node |: ..*", 
                                                                             replacement = "",
                                                                             x = SN.text ),
                                                              pattern = " " ) ) )
    # 4. Extract the parameters for the SN distributions
    SN.calibs <- lapply( X = stringr::str_split( gsub( pattern = "Node.*SN.*\\( | \\)|,",
                                                       replacement = "",
                                                       x = SN.text ),
                                                 pattern = " " ),
                         FUN = as.numeric )
    names( SN.calibs ) <- paste( "t_n", SN.calib.nodes, sep = "" )
    # 5. Name each node accordingly                          
    names( SN.calib.nodes ) <- SN.calibs.nodes.names
    # 6. Add SN calibrations used for these nodes 
    for( i in 1:length( SN.calib.nodes ) ){
      
      names( SN.calibs[[ i ]] )[1] <- "xi"    #position ~ approx. to mean age
      names( SN.calibs[[ i ]] )[2] <- "omega" #scale ~ variance (small number narrow, large wider)
      names( SN.calibs[[ i ]] )[3] <- "alpha" #shape (slanding, going left or right)
      
    }
    # 7. Load RData object with ST calibrations fitted to the post calibrations of 
    #    72sp mammals to find the equivalent node position. The SN calibration has
    #    the same 3 first values
    cat( "Load RData object with ST calibrations fitted to posterior time estimates
       when using the 72sp mammals phylogeny... ...\n")
    load( ST.fitted.dists.RData )
    tmp <- FALSE
    nodes.SN.72sp <- rep( "", length( SN.calibs ) )
    cat( "Finding matching node positions with nodes calibrated with SN distributions... ...\n")
    for( i in 1:length( SN.calibs ) ){
      
      for( j in 1:length( ST.fitted.dists ) ){
        tmp.ndec <- nchar( gsub( pattern = "..*\\.", x = SN.calibs[[i]], replacement = "" ) )
        if( any( tmp.ndec == 4 ) ){
          tmp <- all.equal( round( ST.fitted.dists[[ j ]][1:3], 4 ), round( SN.calibs[[ i ]], 4 ) )
        }else{
          tmp <- all.equal( round( ST.fitted.dists[[ j ]][1:3], 3 ), round( SN.calibs[[ i ]], 3 ) )
        }
        
        if( tmp == TRUE ){
          break
        }
      }
      
      # Add index label and set tmp to FALSE
      nodes.SN.72sp[ i ] <- names( ST.fitted.dists )[j]
      tmp <- FALSE
      
    }
  }
  # 9. Return list with objects
  if( SN == TRUE ){
    cat( "\nTasks done! Return objects with ST and SN info\n")
    tmp.len.ST   <- length( ST.calibs )
    tmp.len.SN   <- length( ST.calibs )
    total.calibs <- tmp.len.ST + tmp.len.SN 
    ST.SN.calibs     <- ST.SN.calib.nodes <- vector( mode = "list", length = total.calibs )
    nodes.all.72sp   <- vector( mode = "character", length = total.calibs )
    for( i in 1:tmp.len.ST ){
      ST.SN.calibs[[ i ]]      <- ST.calibs[[ i ]]
      nodes.all.72sp[i]        <- nodes.72sp[i]
      ST.SN.calib.nodes[[ i ]] <- ST.calib.nodes[[ i ]]
    }
    names( ST.SN.calibs )[1:tmp.len.ST]      <- names( ST.calibs )
    names( ST.SN.calib.nodes )[1:tmp.len.ST] <- names( ST.calib.nodes )
    count <- 0
    for( i in (tmp.len.ST+1):total.calibs ){
      count <- count + 1
      ST.SN.calibs[[ i ]]      <- SN.calibs[[ count ]]
      nodes.all.72sp[i]        <- nodes.SN.72sp[count]
      ST.SN.calib.nodes[[ i ]] <- SN.calib.nodes[[ count ]]
    }
    names( ST.SN.calibs )[(tmp.len.ST+1):total.calibs]      <- names( SN.calibs )
    names( ST.SN.calib.nodes )[(tmp.len.ST+1):total.calibs] <- names( SN.calib.nodes )
    return( list( ST.SN.calibs = ST.SN.calibs, nodes.72sp = nodes.all.72sp, ST.SN.calib.nodes = ST.SN.calib.nodes
                  ) )
  }else{
    cat( "\nTasks done! Return objects with ST info\n")
    return( list( ST.calibs = ST.calibs, nodes.72sp = nodes.72sp, ST.calib.nodes = ST.calib.nodes ))
  }
  
  
}



#=======================================================================#
# FUNCTION TO DEFINE SN CALIBRATIONS 
#-----------------------------------------------------------------------#
# Arguments
#     log_SN   character, path to the log output file printed out 
#              by MCMCtree during the run 
#     SN.fitted.dists.RData  character, path to the RData file output 
#                            when fitting SN distributions to posterior 
#                            time estimates when using 72sp 
#     SN.calib.node.names    character, vector of length as many SN 
#                            calibrations, previously created in the 
#                            script, with the names of each calibrated node
#
define_SN <- function( log_SN, SN.fitted.dists.RData, SN.calib.nodes.names ){
  
  # 1. Read out.txt of one of the runs and extract info pasted above
  cat( "Loading log screen file saved when running MCMCtree... ...\n")
  outMCMCtree <- readLines( con = log_SN )
  # 2. Extract lines with the SN calibrations
  cat( "Generating objects with SN calibrations... ...\n")
  SN.text <- grep( pattern = "Node ..* SN ..*", x = outMCMCtree, value = TRUE )
  # 3. Extract the position of the calibrated nodes
  SN.calib.nodes <- as.numeric( unlist( stringr::str_split( string = gsub( pattern = 
                                                                             "Node  |Node |: ..*", 
                                                                           replacement = "",
                                                                           x = SN.text ),
                                                            pattern = " " ) ) )
  # 4. Extract the parameters for the ST distributions
  SN.calibs <- lapply( X = stringr::str_split( gsub( pattern = "Node.*SN.*\\( | \\)|,",
                                                     replacement = "",
                                                     x = SN.text ),
                                               pattern = " " ),
                       FUN = as.numeric )
  names( SN.calibs ) <- paste( "t_n", SN.calib.nodes, sep = "" )
  # 5. Name each node accordingly                          
  names( SN.calib.nodes ) <- SN.calib.nodes.names
  # 6. Add SN calibrations used for these nodes 
  for( i in 1:length( SN.calib.nodes ) ){
    
    names( SN.calibs[[ i ]] )[1] <- "xi"    # mean ~ approx. to mean age
    names( SN.calibs[[ i ]] )[2] <- "omega" # sd ~ small number narrow, large wider
    names( SN.calibs[[ i ]] )[3] <- "alpha" # gamma1

  }
  # 7. Load RData object with SN calibrations fitted to the post calibrations of 
  #    72sp mammals to find the equivalent node position
  cat( "Load RData object with SN calibrations fitted to posterior time estimates
       when using the 72sp mammals phylogeny... ...\n")
  ## 210719 - SAC
  ## Just in case different names are given to the "SN.fitted.dists.RData" 
  ## as I tried different SN values, I use the `assign` function as
  ## suggested here: https://stackoverflow.com/questions/5577221/how-can-i-load-an-object-into-a-variable-name-that-i-specify-from-an-r-data-file
  #load( SN.fitted.dists.RData )
  assign( 'SN.fitted.dists', get( load( SN.fitted.dists.RData ) ) )
  tmp <- FALSE
  nodes.72sp <- rep( "", length( SN.calibs ) )
  cat( "Finding matching node positions with nodes calibrated with SN distributions... ...\n")
  for( i in 1:length( SN.calibs ) ){
    
    for( j in 1:length( SN.fitted.dists ) ){
      tmp.ndec <- nchar( gsub( pattern = "..*\\.", x = SN.calibs[[i]], replacement = "" ) )
      ## 210719 - SAC
      ## Need to rename "mean", "s.d." and "gamma1" to "xi", "omega", and "alpha"; respectively
      names( SN.fitted.dists[[ j ]] ) <- names( SN.calibs[[ i ]] )
      if( any( tmp.ndec == 4 ) ){
        tmp <- all.equal( round( SN.fitted.dists[[ j ]], 4 ), round( SN.calibs[[ i ]], 4 ) )
      }else{
        tmp <- all.equal( round( SN.fitted.dists[[ j ]], 3 ), round( SN.calibs[[ i ]], 3 ) )
      }
      
      if( tmp == TRUE ){
        break
      }
    }
    
    # Add index label and set tmp to FALSE
    nodes.72sp[ i ] <- names( SN.fitted.dists )[j]
    tmp <- FALSE
    
  }
  # 8. Delete SN.fitted.dists obj
  remove( SN.fitted.dists )
  
  # 9. Return list with objects
  cat( "\nTasks done! Return objects with SN info\n")
  return( list( SN.calibs = SN.calibs, nodes.72sp = nodes.72sp, SN.calib.nodes = SN.calib.nodes ))
  
}



#======================================================#
# FUNCTION: Load data and calculate mean and quantiles # 
#------------------------------------------------------#
# Arguments: 
#      mcmc1  character, path to mcmc.txt from run1 
#      mcmc2  character, path to mcmc.txt from run2
load_data <- function( mcmc1, mcmc2 ){
  
  # 1. Load files and get parameters
  cat( "Load mcmc.txt from run 1... ...\n")
  run1    <- read.table( mcmc1, header = T, sep = "\t" )
  cat( "Load mcmc.txt from run 2... ...\n")
  run2    <- read.table( mcmc2, header = T, sep = "\t" )
  
  # 2. Summarise parameters for both runs 
  cat( "Generate objects with summarised estimated divergence times... ...\n")
  divtimes1 <- run1[,-c(1, dim(run1)[2])]
  divtimes2 <- run2[,-c(1, dim(run2)[2])]
  divtimes  <- divtimes.prior <- rbind( divtimes1, divtimes2 )

  mean_est_divt  <- apply( X = divtimes, MARGIN = 2, FUN = mean )
  quant_est_divt <- apply( X = divtimes, MARGIN = 2, FUN = quantile, probs = c( 0.025,0.975) )
  quant_est_divt <- t( quant_est_divt )
  all.equal( quant_est_divt[1,], quantile( divtimes[,1], probs = c( 0.025, 0.975 ) ) )
  
  # 3. Return object 
  cat( "\nTasks done! Return objects\n")
  return( list( divt1 = divtimes1, divt2 = divtimes2, divt = divtimes, 
                mean_divt = mean_est_divt, quant_divt = quant_est_divt ) )
}



#=======================================================================#
# FUNCTION: Plotting the posterior dists, ST analytical, and post. 72sp #
#-----------------------------------------------------------------------#
## 211019 -- Now it checks if objects passed with calibrations have 3 or 
##           4 elements (ST or SN)
# Arguments: 
#    post.cal       data.frame, posterior divergence times with the data of Seq.Bayes.S2.
#    ind.post.cal   numeric, vector of length number of node calibrations used in the study.
#    ST.cal         list, vector of length number of node calibrations used in the study, and 
#                   each entry has the corresponding ST distribution for the node position.
#                   E.g., object `ST.calibs.??sp`.
#    ST.cal.nodes   character, vector of length number of node calibrations used in the study. 
#                   Each element contains the node position in the format "t_nXXX".
#                   E.g., object `names( ST.calibs.??sp )`.
#    post.72sp      data.frame, posterior divergence times with the data of Seq.Bayes.S1, 72sp.
#    ind.post.72sp  numeric, vector of length number of node calibrations used in this study. 
#                   Each element is the column in `post.72sp` that corresponds to the 
#                   node for which we want to retrieve the samples.
#    names.calibs   character, vector of length number of node calibrations used in this study.
#                   Each element is the name of the node that has been calibrated in the phylogeny.
#    p.nrow         numeric, vector of length 1 with the number of rows to use when plotting 
#                   the graphs.
#    p.ncol         numeric, vector of length 1 with the number of columns to use when plotting 
#                   the graphs.
#    out            boolean, TRUE if you want to save the plot. If you just want to display it, 
#                   then FALSE.
#    outname        character, name for the output file if you wanted to have one with the plots.
plot_func <- function( post.cal = divtimes, ind.post.cal = ind.post.cal,
                       ST.cal = ST.calibs, ST.cal.nodes = names( ST.calibs ),
                       post.72sp = divtimes.post72sp, ind.post.72sp = ind.post.72sp,
                       names.calibs = ST.calib.nodes.names,
                       p.nrow = p.nrow, p.ncol = p.ncol, legend = TRUE, 
                       out = TRUE, outname ){
  
  if( length( ind.post.cal) != length( ind.post.72sp ) ){
    stop( "Same calibrations required for post with 72sp and current" )
  }
  
  if( out == TRUE ){
    # The following code enables to save the same plot 
    # in both pdf and png format
    pdf( file = paste( outname, ".pdf", sep = "" ),
         paper = "a4" )
    curr_plot <- dev.cur()
    png( filename = paste( outname, ".png", sep = "" ),
         width = 1024, height = 768 )
    dev.control( "enable" )
    par( mfrow = c( p.nrow,p.ncol ), mai = c(0,0,0,0) )
  }else{
    par( mfrow = c( p.nrow,p.ncol ), mai = c(0,0,0,0) )
  }
  flag <- rep( 0 , length( ind.post.cal ) )
  for( i in 1:length( ind.post.cal ) ){
    
    plot( density(  post.cal[,ind.post.cal[i]], adj = 1 ), 
          main = paste( ST.cal.nodes[i], " = ", names.calibs[i], sep = "" ),
          xlab = '', ylab = '' )
    if( length( ST.cal[[ i ]] ) == 3 ){
      curve( dsn( x, xi =  ST.cal[[ i ]][1], omega =  ST.cal[[ i ]][2],
                  alpha =  ST.cal[[ i ]][3] ),
             from = 0, to = 5,
             n = 1e4, add = TRUE, col = "red" ) 
      flag[i] <- 1
    }else{
      curve( dst( x, xi =  ST.cal[[ i ]][1], omega =  ST.cal[[ i ]][2],
                  alpha =  ST.cal[[ i ]][3], nu =  ST.cal[[ i ]][4] ),
             from = 0, to = 5,
             n = 1e4, add = TRUE, col = "red" ) 
      flag[i] <- 0
    }
    lines ( density(  post.72sp[,ind.post.72sp[i]], adj = 1 ),
            col = "darkgrey" )
    if( i == 1 ){
      if( legend == TRUE ){
        # Add legend
        if( 1 %in% flag ){
          info.legend <- c( "Post.", "Analyt. ST/SN, 72sp",
                            "Post., 72sp")
        }else{
          info.legend <- c( "Post.", "Analyt. ST, 72sp",
                            "Post., 72sp")
        }
        
        col.legend  <- c( "black", "red", "grey" )
        legend( x=2.2, y = 1.5, legend = info.legend, col = col.legend,
                lty = 1, bty = 'n' )
      }
      
    }
  }
  if( out == TRUE ){
    # Copy the plot in pdf to png and then close 
    # both graphics
    dev.copy( which = curr_plot )
    dev.off()
    dev.off()
  }
  
}



#=======================================================================#
# FUNCTION: Plotting the posterior dists, SN analytical, and post. 72sp #
#-----------------------------------------------------------------------#
# Arguments: 
#    post.cal       data.frame, posterior divergence times with the data of Seq.Bayes.S2.
#    ind.post.cal   numeric, vector of length number of node calibrations used in the study.
#    SN.cal         list, vector of length number of node calibrations used in the study, and 
#                   each entry has the corresponding SN distribution for the node position.
#                   E.g., object `SN.calibs.??sp`.
#    SN.cal.nodes   character, vector of length number of node calibrations used in the study. 
#                   Each element contains the node position in the format "t_nXXX".
#                   E.g., object `names( SN.calibs.??sp )`.
#    ST.cal         list, vector of length number of node calibrations used in the study, and 
#                   each entry has the corresponding ST distribution for the node position.
#                   E.g., object `ST.calibs.??sp`.
#    ST.cal.nodes   character, vector of length number of node calibrations used in the study. 
#                   Each element contains the node position in the format "t_nXXX".
#                   E.g., object `names( ST.calibs.??sp )`.
#    post.72sp      data.frame, posterior divergence times with the data of Seq.Bayes.S1, 72sp.
#    ind.post.72sp  numeric, vector of length number of node calibrations used in this study. 
#                   Each element is the column in `post.72sp` that corresponds to the 
#                   node for which we want to retrieve the samples.
#    names.calibs   character, vector of length number of node calibrations used in this study.
#                   Each element is the name of the node that has been calibrated in the phylogeny.
#    p.nrow         numeric, vector of length 1 with the number of rows to use when plotting 
#                   the graphs.
#    p.ncol         numeric, vector of length 1 with the number of columns to use when plotting 
#                   the graphs.
#    out            boolean, TRUE if you want to save the plot. If you just want to display it, 
#                   then FALSE.
#    outname        character, name for the output file if you wanted to have one with the plots.
plot_func_SN <- function( post.cal = divtimes, ind.post.cal = ind.post.cal,
                          SN.cal = SN.calibs, SN.cal.nodes = names( SN.calibs ),
                          ST.cal = ST.calibs, ST.cal.nodes = names( ST.calibs ),
                          post.72sp = divtimes.post72sp, ind.post.72sp = ind.post.72sp,
                          names.calibs = SN.calib.nodes.names,
                          p.nrow = p.nrow, p.ncol = p.ncol, legend = TRUE, 
                          out = TRUE, outname ){
  
  if( length( ind.post.cal) != length( ind.post.72sp ) ){
    stop( "Same calibrations required for post with 72sp and current" )
  }
  
  if( out == TRUE ){
    # The following code enables to save the same plot 
    # in both pdf and png format
    pdf( file = paste( outname, ".pdf", sep = "" ),
         paper = "a4" )
    curr_plot <- dev.cur()
    png( filename = paste( outname, ".png", sep = "" ),
         width = 1024, height = 768 )
    dev.control( "enable" )
    par( mfrow = c( p.nrow,p.ncol ), mai = c(0,0,0,0) )
  }else{
    par( mfrow = c( p.nrow,p.ncol ), mai = c(0,0,0,0) )
  }
  for( i in 1:length( ind.post.cal ) ){
    
    plot( density(  post.cal[,ind.post.cal[i]], adj = 1 ), 
          main = paste( SN.cal.nodes[i], " = ", names.calibs[i], sep = "" ),
          xlab = '', ylab = '' )
    curve( dsn( x, xi =  SN.cal[[ i ]][1], omega =  SN.cal[[ i ]][2],
                alpha =  SN.cal[[ i ]][3] ),
           from = 0, to = 5,
           n = 1e4, add = TRUE, col = "red" ) 
    curve( dst( x, xi =  ST.cal[[ i ]][1], omega =  ST.cal[[ i ]][2],
                alpha =  ST.cal[[ i ]][3], nu = ST.cal[[ i ]][4] ),
           from = 0, to = 5,
           n = 1e4, add = TRUE, col = "purple" ) 
    lines ( density(  post.72sp[,ind.post.72sp[i]], adj = 1 ),
            col = "darkgrey" )
    if( i == 1 ){
      if( legend == TRUE ){
        # Add legend
        info.legend <- c( "Post.", "Analyt. SN, 72sp",
                          "Analyt. ST, 72sp", "Post., 72sp")
        col.legend  <- c( "black", "red", "purple", "grey" )
        legend( x=2.2, y = 1.5, legend = info.legend, col = col.legend,
                lty = 1, bty = 'n' )
      }
      
    }
  }
  if( out == TRUE ){
    # Copy the plot in pdf to png and then close 
    # both graphics
    dev.copy( which = curr_plot )
    dev.off()
    dev.off()
  }
  
}



#============================================================================#
# FUNCTION: Plotting the posterior dists, SN & ST analytical, and post. 72sp #
#----------------------------------------------------------------------------#
# Arguments: 
#    post.cal         data.frame, posterior divergence times with the data of Seq.Bayes.S2.
#    ind.post.cal     numeric, vector of length number of node calibrations used in the study.
#    SN.ST.cal        list, vector of length number of node calibrations used in the study, and 
#                     each entry has the corresponding SN & ST distribution for the node position.
#                     E.g., object `SN.ST.calibs.??sp`.
#    SN.ST.cal.nodes  character, vector of length number of node calibrations used in the study. 
#                     Each element contains the node position in the format "t_nXXX".
#                     E.g., object `names( SN.ST.calibs.??sp )`.
#    ST.cal           list, vector of length number of node calibrations used in the study, and 
#                     each entry has the corresponding SN & ST distributions for the node position.
#                     E.g., object `SN.ST.calibs.??sp`.
#    ST.cal.nodes     character, vector of length number of node calibrations used in the study. 
#                     Each element contains the node position in the format "t_nXXX".
#                     E.g., object `names( SN.ST.calibs.??sp )`.
#    post.72sp        data.frame, posterior divergence times with the data of Seq.Bayes.S1, 72sp.
#    ind.post.72sp    numeric, vector of length number of node calibrations used in this study. 
#                     Each element is the column in `post.72sp` that corresponds to the 
#                     node for which we want to retrieve the samples.
#    names.calibs     character, vector of length number of node calibrations used in this study.
#                     Each element is the name of the node that has been calibrated in the phylogeny.
#    p.nrow           numeric, vector of length 1 with the number of rows to use when plotting 
#                     the graphs.
#    p.ncol           numeric, vector of length 1 with the number of columns to use when plotting 
#                     the graphs.
#    out              boolean, TRUE if you want to save the plot. If you just want to display it, 
#                     then FALSE.
#    outname          character, name for the output file if you wanted to have one with the plots.
plot_func_SN_ST <- function( post.cal = divtimes, ind.post.cal = ind.post.cal,
                             SN.ST.cal = SN.ST.calibs, SN.ST.cal.nodes = names( SN.ST.calibs ),
                             post.72sp = divtimes.post72sp, ind.post.72sp = ind.post.72sp,
                             names.calibs = SN.ST.calib.nodes.names,
                             p.nrow = p.nrow, p.ncol = p.ncol, legend = TRUE, 
                             out = TRUE, outname ){
  
  if( length( ind.post.cal) != length( ind.post.72sp ) ){
    stop( "Same calibrations required for post with 72sp and current" )
  }
  
  if( out == TRUE ){
    # The following code enables to save the same plot 
    # in both pdf and png format
    pdf( file = paste( outname, ".pdf", sep = "" ),
         paper = "a4" )
    curr_plot <- dev.cur()
    png( filename = paste( outname, ".png", sep = "" ),
         width = 1024, height = 768 )
    dev.control( "enable" )
    par( mfrow = c( p.nrow,p.ncol ), mai = c(0,0,0,0) )
  }else{
    par( mfrow = c( p.nrow,p.ncol ), mai = c(0,0,0,0) )
  }
  for( i in 1:length( ind.post.cal ) ){
    
    plot( density(  post.cal[,ind.post.cal[i]], adj = 1 ), 
          main = paste( SN.ST.cal.nodes[i], " = ", names.calibs[i], sep = "" ),
          xlab = '', ylab = '' )
    if( length( SN.ST.cal[[ i ]] ) == 3 ){
      curve( dsn( x, xi =  SN.ST.cal[[ i ]][1], omega =  SN.ST.cal[[ i ]][2],
                  alpha =  SN.ST.cal[[ i ]][3] ),
             from = 0, to = 5,
             n = 1e4, add = TRUE, col = "purple" ) 
    }
    else if( length( SN.ST.cal[[ i ]] ) == 4 ){
      curve( dst( x, xi =  SN.ST.cal[[ i ]][1], omega =  SN.ST.cal[[ i ]][2],
                  alpha =  SN.ST.cal[[ i ]][3], nu = SN.ST.cal[[ i ]][4] ),
             from = 0, to = 5,
             n = 1e4, add = TRUE, col = "red" ) 
    }
    lines ( density(  post.72sp[,ind.post.72sp[i]], adj = 1 ),
            col = "darkgrey" )
    if( legend == TRUE ){
        # Add legend
      if( length( SN.ST.cal[[ i ]] ) == 3 ){
        info.legend <- c( "Post.", "Analyt. SN, 72sp",
                          "Post., 72sp")
        col.legend  <- c( "black", "purple", "grey" )
      } else if ( length( SN.ST.cal[[ i ]] ) == 4 ){ 
        info.legend <- c( "Post.", "Analyt. ST, 72sp",
                          "Post., 72sp")
        col.legend  <- c( "black", "red", "grey" )
      }
      legend( "bottom", legend = info.legend, col = col.legend,
              lty = 1, bty = 'n' )
    }
  }
  if( out == TRUE ){
    # Copy the plot in pdf to png and then close 
    # both graphics
    dev.copy( which = curr_plot )
    dev.off()
    dev.off()
  }
  
}

#=======================================================================#
# FUNCTION: Plotting the posterior dists, ST analytical, and post. 72sp #
#-----------------------------------------------------------------------#
## 211019 -- Now it checks if objects passed with calibrations have 3 or 
##           4 elements (ST or SN)
# Arguments: 
#    post.cal.1     data.frame, posterior divergence times with the data of Seq.Bayes.S2, r1.
#    post.cal.2     data.frame, posterior divergence times with the data of Seq.Bayes.S2, r2.
#    ind.post.cal   numeric, vector of length number of node calibrations used in the study.
#    ST.cal         list, vector of length number of node calibrations used in the study, and 
#                   each entry has the corresponding ST distribution for the node position.
#                   E.g., object `ST.calibs.??sp`.
#    ST.cal.nodes   character, vector of length number of node calibrations used in the study. 
#                   Each element contains the node position in the format "t_nXXX".
#                   E.g., object `names( ST.calibs.??sp )`.
#    post.72sp      data.frame, posterior divergence times with the data of Seq.Bayes.S1, 72sp.
#    ind.post.72sp  numeric, vector of length number of node calibrations used in this study. 
#                   Each element is the column in `post.72sp` that corresponds to the 
#                   node for which we want to retrieve the samples.
#    names.calibs   character, vector of length number of node calibrations used in this study.
#                   Each element is the name of the node that has been calibrated in the phylogeny.
#    p.nrow         numeric, vector of length 1 with the number of rows to use when plotting 
#                   the graphs.
#    p.ncol         numeric, vector of length 1 with the number of columns to use when plotting 
#                   the graphs.
#    out            boolean, TRUE if you want to save the plot. If you just want to display it, 
#                   then FALSE.
#    outname        character, name for the output file if you wanted to have one with the plots.
plot_func_mcmc <- function( post.cal.1, post.cal.2, 
                            ind.post.cal = ind.post.cal,
                            ST.cal, ST.cal.nodes,
                            post.72sp = divtimes.post72sp, ind.post.72sp = ind.post.72sp,
                            names.calibs = ST.calib.nodes.names,
                            p.nrow = p.nrow, p.ncol = p.ncol, legend = TRUE, out = TRUE, outname ){
  
  if( length( ind.post.cal) != length( ind.post.72sp ) ){
    stop( "Same calibrations required for post with 72sp and current" )
  }
  
  if( out == TRUE ){
    # The following code enables to save the same plot 
    # in both pdf and png format
    pdf( file = paste( outname, ".pdf", sep = "" ),
         paper = "a4" )
    curr_plot <- dev.cur()
    png( filename = paste( outname, ".png", sep = "" ),
         width = 1024, height = 768 )
    dev.control( "enable" )
    par( mfrow = c( p.nrow,p.ncol ), mai = c(0,0,0,0) )
  }else{
    par( mfrow = c( p.nrow,p.ncol ), mai = c(0,0,0,0) )
  }
  flag <- rep( 0 , length( ind.post.cal ) )
  for( i in 1:length( ind.post.cal ) ){
    
    plot( density( post.cal.1[,ind.post.cal[i]], adj = 1 ), 
          main = paste( ST.cal.nodes[i], " = ", names.calibs[i], sep = "" ),
          lty = 3, xlab = '', ylab = '' )
    lines ( density( post.cal.2[,ind.post.cal[i]], adj = 1 ),
            lty = 5, col = "purple" )
    if( length( ST.cal[[ i ]] == 3 ) ){
      curve( dsn( x, xi =  ST.cal[[ i ]][1], omega =  ST.cal[[ i ]][2],
                  alpha =  ST.cal[[ i ]][3] ),
             from = 0, to = 5,
             n = 1e4, add = TRUE, col = "darkgray" ) 
      flag[i] <- 1
    }else{
      curve( dst( x, xi =  ST.cal[[ i ]][1], omega =  ST.cal[[ i ]][2],
                  alpha =  ST.cal[[ i ]][3], nu =  ST.cal[[ i ]][4] ),
             from = 0, to = 5,
             n = 1e4, add = TRUE, col = "darkgray" ) 
      flag[i] <- 0
    }
    
    
    if( i == 1 ){
      if( legend == TRUE ){
        # Add legend
        if( 1 %in% flag ){
          info.legend <- c( "post. r1", "post. r2",
                            "Analy. ST/SN, 72sp")
        }else{
          info.legend <- c( "post. r1", "post. r2",
                            "Analy. ST, 72sp")
        }
        col.legend  <- c( "black", "purple", "darkgrey" )
        legend( x=2.2, y = 1.5, legend = info.legend, col = col.legend,
                lty = 1, bty = 'n' )
      }
      
    }
  }
  if( out == TRUE ){
    # Copy the plot in pdf to png and then close 
    # both graphics
    dev.copy( which = curr_plot )
    dev.off()
    dev.off()
  }
  
}



#=======================================================================#
# FUNCTION: Plotting the posterior dists, SN analytical, and post. 72sp #
#-----------------------------------------------------------------------#
# Arguments: 
#    post.cal.1     data.frame, posterior divergence times with the data of Seq.Bayes.S2, r1.
#    post.cal.2     data.frame, posterior divergence times with the data of Seq.Bayes.S2, r2.
#    ind.post.cal   numeric, vector of length number of node calibrations used in the study.
#    SN.cal         list, vector of length number of node calibrations used in the study, and 
#                   each entry has the corresponding SN distribution for the node position.
#                   E.g., object `SN.calibs.??sp`.
#    SN.cal.nodes   character, vector of length number of node calibrations used in the study. 
#                   Each element contains the node position in the format "t_nXXX".
#                   E.g., object `names( SN.calibs.??sp )`.
#    post.72sp      data.frame, posterior divergence times with the data of Seq.Bayes.S1, 72sp.
#    ind.post.72sp  numeric, vector of length number of node calibrations used in this study. 
#                   Each element is the column in `post.72sp` that corresponds to the 
#                   node for which we want to retrieve the samples.
#    names.calibs   character, vector of length number of node calibrations used in this study.
#                   Each element is the name of the node that has been calibrated in the phylogeny.
#    p.nrow         numeric, vector of length 1 with the number of rows to use when plotting 
#                   the graphs.
#    p.ncol         numeric, vector of length 1 with the number of columns to use when plotting 
#                   the graphs.
#    out            boolean, TRUE if you want to save the plot. If you just want to display it, 
#                   then FALSE.
#    outname        character, name for the output file if you wanted to have one with the plots.
plot_func_mcmc_SN <- function( post.cal.1, post.cal.2, 
                               ind.post.cal = ind.post.cal,
                               SN.cal, SN.cal.nodes,
                               post.72sp = divtimes.post72sp, ind.post.72sp = ind.post.72sp,
                               names.calibs,
                               p.nrow = p.nrow, p.ncol = p.ncol, legend = TRUE, out = TRUE, outname ){
  
  if( length( ind.post.cal) != length( ind.post.72sp ) ){
    stop( "Same calibrations required for post with 72sp and current" )
  }
  
  if( out == TRUE ){
    # The following code enables to save the same plot 
    # in both pdf and png format
    pdf( file = paste( outname, ".pdf", sep = "" ),
         paper = "a4" )
    curr_plot <- dev.cur()
    png( filename = paste( outname, ".png", sep = "" ),
         width = 1024, height = 768 )
    dev.control( "enable" )
    par( mfrow = c( p.nrow,p.ncol ), mai = c(0,0,0,0) )
  }else{
    par( mfrow = c( p.nrow,p.ncol ), mai = c(0,0,0,0) )
  }
  for( i in 1:length( ind.post.cal ) ){
    
    plot( density( post.cal.1[,ind.post.cal[i]], adj = 1 ), 
          main = paste( SN.cal.nodes[i], " = ", names.calibs[i], sep = "" ),
          lty = 3, xlab = '', ylab = '' )
    lines ( density( post.cal.2[,ind.post.cal[i]], adj = 1 ),
            lty = 5, col = "purple" )
    curve( dsn( x, xi =  SN.cal[[ i ]][1], omega =  SN.cal[[ i ]][2],
                alpha =  SN.cal[[ i ]][3] ),
           from = 0, to = 5,
           n = 1e4, add = TRUE, col = "darkgray" ) 
    
    if( i == 1 ){
      if( legend == TRUE ){
        # Add legend
        info.legend <- c( "post. r1", "post. r2",
                          "Analy. SN, 72sp")
        col.legend  <- c( "black", "purple", "darkgrey" )
        legend( x=2.2, y = 1.5, legend = info.legend, col = col.legend,
                lty = 1, bty = 'n' )
      }
      
    }
  }
  if( out == TRUE ){
    # Copy the plot in pdf to png and then close 
    # both graphics
    dev.copy( which = curr_plot )
    dev.off()
    dev.off()
  }
  
}



#============================================================================#
# FUNCTION: Plotting the posterior dists, SN & ST analytical, and post. 72sp #
#----------------------------------------------------------------------------#
# Arguments: 
#    post.cal.1        data.frame, posterior divergence times with the data of Seq.Bayes.S2, r1.
#    post.cal.2        data.frame, posterior divergence times with the data of Seq.Bayes.S2, r2.
#    ind.post.cal      numeric, vector of length number of node calibrations used in the study.
#    SN.ST.cal         list, vector of length number of node calibrations used in the study, and 
#                      each entry has the corresponding SN & ST distributions for the node position.
#                      E.g., object `SN.ST.calibs.??sp`.
#    SN.ST.cal.nodes   character, vector of length number of node calibrations used in the study. 
#                      Each element contains the node position in the format "t_nXXX".
#                      E.g., object `names( SN.ST.calibs.??sp )`.
#    post.72sp         data.frame, posterior divergence times with the data of Seq.Bayes.S1, 72sp.
#    ind.post.72sp     numeric, vector of length number of node calibrations used in this study. 
#                      Each element is the column in `post.72sp` that corresponds to the 
#                      node for which we want to retrieve the samples.
#    names.calibs      character, vector of length number of node calibrations used in this study.
#                      Each element is the name of the node that has been calibrated in the phylogeny.
#    p.nrow            numeric, vector of length 1 with the number of rows to use when plotting 
#                      the graphs.
#    p.ncol            numeric, vector of length 1 with the number of columns to use when plotting 
#                      the graphs.
#    out               boolean, TRUE if you want to save the plot. If you just want to display it, 
#                      then FALSE.
#    outname           character, name for the output file if you wanted to have one with the plots.
plot_func_mcmc_SN_ST <- function( post.cal.1, post.cal.2, 
                                  ind.post.cal = ind.post.cal,
                                  SN.ST.cal = SN.ST.calibs, SN.ST.cal.nodes = names( SN.calibs ),
                                  post.72sp = divtimes.post72sp, ind.post.72sp = ind.post.72sp,
                                  names.calibs = SN.ST.calib.nodes.names,
                                  p.nrow = p.nrow, p.ncol = p.ncol, legend = TRUE, out = TRUE, outname ){
  
  if( length( ind.post.cal) != length( ind.post.72sp ) ){
    stop( "Same calibrations required for post with 72sp and current" )
  }
  
  if( out == TRUE ){
    # The following code enables to save the same plot 
    # in both pdf and png format
    pdf( file = paste( outname, ".pdf", sep = "" ),
         paper = "a4" )
    curr_plot <- dev.cur()
    png( filename = paste( outname, ".png", sep = "" ),
         width = 1024, height = 768 )
    dev.control( "enable" )
    par( mfrow = c( p.nrow,p.ncol ), mai = c(0,0,0,0) )
  }else{
    par( mfrow = c( p.nrow,p.ncol ), mai = c(0,0,0,0) )
  }
  for( i in 1:length( ind.post.cal ) ){
    
    plot( density( post.cal.1[,ind.post.cal[i]], adj = 1 ), 
          main = paste( SN.ST.cal.nodes[i], " = ", names.calibs[i], sep = "" ),
          lty = 3, xlab = '', ylab = '' )
    if( length( SN.ST.cal[[ i ]] ) == 3 ){
      lines ( density( post.cal.2[,ind.post.cal[i]], adj = 1 ),
              lty = 5, col = "purple" )
      curve( dsn( x, xi =  SN.ST.cal[[ i ]][1], omega = SN.ST.cal[[ i ]][2],
                  alpha =  SN.ST.cal[[ i ]][3] ),
             from = 0, to = 5,
             n = 1e4, add = TRUE, col = "darkgray" ) 
    }else if( length( SN.ST.cal[[ i ]] ) ==  4){
      lines ( density( post.cal.2[,ind.post.cal[i]], adj = 1 ),
              lty = 5, col = "red" )
      curve( dst( x, xi =  SN.ST.cal[[ i ]][1], omega = SN.ST.cal[[ i ]][2],
                  alpha =  SN.ST.cal[[ i ]][3], nu = SN.ST.cal[[ i ]][4] ),
             from = 0, to = 5,
             n = 1e4, add = TRUE, col = "darkgray" ) 
    }
    
    if( legend == TRUE ){
        # Add legend
        if( length( SN.ST.cal[[ i ]] ) == 3 ){
          info.legend <- c( "post. r1", "post. r2",
                            "Analy. SN, 72sp")
          col.legend  <- c( "black", "purple", "darkgrey" )
        } else if( length( SN.ST.cal[[ i ]] ) == 4 ){
          info.legend <- c( "post. r1", "post. r2",
                            "Analy. ST, 72sp")
          col.legend  <- c( "black", "red", "darkgrey" )
        }
        # legend( x=2.2, y = 1.5, legend = info.legend, col = col.legend,
        #         lty = 1, bty = 'n' )
        legend( "bottom", legend = info.legend, col = col.legend,
                lty = 1, bty = 'n' )
    }
  }
  if( out == TRUE ){
    # Copy the plot in pdf to png and then close 
    # both graphics
    dev.copy( which = curr_plot )
    dev.off()
    dev.off()
  }
  
}

