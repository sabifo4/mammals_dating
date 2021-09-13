#======================#
#  0. CLEAN ENVIRONMENT
#======================# 
rm( list = ls( ) )

#===========================#
#  1. SET WORKING DIRECTORY
#===========================#
# Load package to help find the path to this source file 
library(rstudioapi) 
# Get the path of current open file
path_to_file <- getActiveDocumentContext()$path 
# Get working directory path
wd      <- paste( dirname( path_to_file ), "/", sep = "" )
wd.name <- dirname( path_to_file )
# Set wd
setwd( wd )

#==========================================#
# 2. NAVIGATE THROUGH DIRS AND SCAN TREES 
#==========================================#
# Get a vector with the directory names 
dirs.names <- list.dirs(path = ".", full.names = F)

# If an empty dir has been caught, erase it from
# the character vector
empty.dirs <- rep( 0, length( dirs.names ) )
for( i in seq( 1:length( dirs.names ) ) ){
  
  # Find empty dirs 
  if( dirs.names[ i ] == "" | dirs.names[ i ] == "." ){
    empty.dirs[ i ] <- i
  }
  # Last iteration remove empty dirs
  if( i == length( dirs.names ) ){
    rm.ind     <- which( empty.dirs == 0 )
    empty.dirs <- empty.dirs[-rm.ind]
    if ( length( empty.dirs ) != 0 ){
      dirs.names <- dirs.names[ -empty.dirs ]
    }
  }
  
}

# Get trees 
trees.list <- vector( mode = "list", 7 )
names( trees.list ) <- dirs.names

for( i in 1:length( dirs.names ) ){
  
  # Read path to tree file
  tree.file <- list.files( paste( wd, dirs.names[i], "/", sep = "" ), pattern = "calibrated.tree$" )
  # Read tree and put it in corresponding entry in the list 
  trees.list[[ i ]] <- readLines( paste( dirs.names[i], "/", tree.file, sep = "" ) )
  # Get phylip header for later, only with first it 
  if( i == 1 ){
    phylip.header <- trees.list[[ i ]][1]
  }
  # Remove first line with PHYLIP header and then replace square brackets with nothing
  trees.list[[ i ]] <- trees.list[[ i ]][-1]
  # Make sure the name corresponds to the loaded tree !
  names( trees.list )[i] <- dirs.names[i]
  
}

#======================================#
#    3. GET FILE WITH CALIBRATIONS 
#======================================#
# Read tree and rename cols
calibrations <- read.table( file = paste( wd, "Calibrations_mammalia.txt", sep = "" ),
                            stringsAsFactors = F, sep = "|", blank.lines.skip = T )
colnames( calibrations ) <- c( "name", "MCMCtree calib" )

#======================================#
#  4. REPLACE NAMES WITH CALIBRATIONS
#======================================#
# Replace calibration names with corresponding calibration
trees.calibs <- trees.list
for( i in 1:length( trees.list ) ){
  for( j in 1:dim( calibrations )[1] ){
    #cat( "Calibrations for tree ", names( trees.calibs )[i], "is: ", calibrations[j,1], calibrations[j,2], "\n" )
    trees.calibs[[ i ]] <- gsub( pattern = paste0("\\[",calibrations[j,1],"\\]"),
                                 x = trees.calibs[[ i ]],
                                 replacement = paste( "'", calibrations[j,2], "'", sep = "" ) )
    
  
  }
}


#=======================#
#  5. WRITE NEW FILES !
#=======================#
for( i in 1:length( trees.calibs ) ){
  
  tmp.tree.name <- gsub( pattern = paste( "0", i, "_", sep = "" ), replacement = "", 
                         x = names( trees.calibs )[i] )
  write( x = phylip.header,
         file = paste( dirs.names[i], "/", "72sp_", tmp.tree.name,
                       "_MCMCtree_calib_v2.tree", sep = "" ) )
  write( x = trees.calibs[[ i ]],
         file = paste( dirs.names[i], "/", "72sp_", tmp.tree.name,
                       "_MCMCtree_calib_v2.tree", sep = "" ),
         append = TRUE )
  
}

