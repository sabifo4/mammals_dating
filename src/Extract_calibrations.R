# Extract calibrations from the log file saved during MCMCtree
#
# Arguments: 
#      log_ST                 character, path to the location of the log file saved
#                             during and MCMC run when running MCMCtree. It expects 
#                             this specific format.
#      ST.calib.nodes.names   Rdata, object previously loaded in the code that contains
#                             the names of the ST calibrations. This location can be found 
#                             inside `scripts/inp_dat` and it is the same generated when 
#                             checking there is no conflict when running MCMCtree sampling
#                             from the prior without data
which_STnames <- function( log_ST, ST.calib.nodes.names ){
  
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
  
  
  # 9. Return list with objects
  cat( "Tasks done! Return objects with ST info\n")
  return( list( ST.calibs = ST.calibs, ST.calib.nodes = ST.calib.nodes, names = names( ST.calib.nodes ) ) )
  
}