#---------------------#
#  CLEAN ENVIRONMENT  #
#---------------------#
rm( list = ls( ) )

#-------------------------#
#  SET WORKING DIRECTORY  #
#-------------------------#
# Load package to help find the path to this source file 
library(rstudioapi) 
# Get the path of current open file
path_to_file <- getActiveDocumentContext()$path 
# Get working directory path
wd      <- paste( dirname( path_to_file ), "/", sep = "" )
wd.name <- dirname( path_to_file )
# Set wd
setwd( wd )
# Set name for output tree
tname <- "MAMMALS_atlantogenata_tarver2016"

#------------#
# SCAN TREES #
#------------#
# Get trees 
tree_original     <- readLines( con = "01_trees_with_ST_calibs/Rinp/FigTree_72sp_nodelabels.tree" )
tree              <- tree_2              <- tree_original
tree_71ST         <- tree_71ST_2         <- tree_original
tree_71ST_rounded <- tree_71ST_rounded_2 <- tree_original

#------------------------#
# GET FILE WITH node.cal #
#------------------------#
# Read tree and rename cols
node.cal <- read.table( file = "01_trees_with_ST_calibs/Rinp/Node_calibrations_mammalia_72sp.txt",
                        stringsAsFactors = F, sep = "|", blank.lines.skip = T, na.strings = "NAN" )
colnames( node.cal ) <- c( "name", "Uniform.cal", "Node.calib" )

ST.fit <- read.table( file = "00_fitST/Rout/ST.fitted.dists.G2.40.tsv",
                      header = TRUE, stringsAsFactors = F, sep = "\t" )
rownames( ST.fit ) <- as.numeric( gsub( pattern =  "t_n", replacement = "",
                                        x = row.names( ST.fit ) ) )

#-------------------------------#
#  REPLACE NAMES WITH node.cal  #
#-------------------------------#

## ---- NOT USED ANYMORE ---- ##
# USE 32 CALIBRATIONS ONLY 
# Replace calibration names with corresponding calibration
# for( i in sort( node.cal$Node.calib, decreasing = TRUE ) ){
#   
#   # Find calibration and replace in tree with ST-fit
#   index.ST <- which( i == as.numeric( rownames( ST.fit ) ) )
#   tree     <- gsub( pattern = paste0(" ",i," "),
#                     x = tree, replacement = paste( " '", ST.fit[index.ST,5], "' ", sep = "" ) )
#   
#   # Last step -- update sp labels and remove nodes without calibrations
#   if( i == min( node.cal$Node.calib ) ){
#     num_sp <- 72:1
#     for( j in num_sp ){
#       tree <- gsub( pattern = paste0( j,"_" ), x = tree, replacement = "" )
#     }
#     uncalibrated.nodes <- setdiff( as.numeric( rownames( ST.fit ) ), sort( node.cal$Node.calib, decreasing = TRUE ) )
#     for( k in uncalibrated.nodes ){
#       tree <- gsub( pattern = paste0(" ",k," "),
#                     x = tree, replacement = "" )
#     }
#   }
# 
# }
## ---- NOT USED ANYMORE ---- ##

# USE 71 CALIBRATIONS ONLY 
# Replace calibration names with corresponding calibration
node.labels.ST <- sort( as.numeric( rownames( ST.fit ) ), decreasing = TRUE )
for( i in node.labels.ST ){
  
  # Find calibration and replace in tree with ST-fit
  index.ST <- which( i == as.numeric( rownames( ST.fit ) ) )
  tree_71ST     <- gsub( pattern = paste0(" ",i," "),
                    x = tree_71ST, replacement = paste( " '", ST.fit[index.ST,5], "' ", sep = "" ) )
  
  # Last step -- update sp labels
  if( i == min( node.cal$Node.calib ) ){
    num_sp <- 72:1
    for( j in num_sp ){
      tree_71ST <- gsub( pattern = paste0( j,"_" ), x = tree_71ST, replacement = "" )
    }
  }
  
}


# USE 71 CALIBRATIONS ONLY - rounded
# Replace calibration names with corresponding calibration
node.labels.ST <- sort( as.numeric( rownames( ST.fit ) ), decreasing = TRUE )
for( i in node.labels.ST ){
  
  # Find calibration and replace in tree with ST-fit
  index.ST <- which( i == as.numeric( rownames( ST.fit ) ) )
  tree_71ST_rounded <- gsub( pattern = paste0(" ",i," "),
                             x = tree_71ST_rounded, replacement = paste( " '", ST.fit[index.ST,6], "' ", sep = "" ) )
  
  # Last step -- update sp labels
  if( i == min( node.cal$Node.calib ) ){
    num_sp <- 72:1
    for( j in num_sp ){
      tree_71ST_rounded <- gsub( pattern = paste0( j,"_" ), x = tree_71ST_rounded, replacement = "" )
    }
  }
  
}

#-----------------#
# WRITE NEW FILES #
#-----------------#
num_sp <- c( stringr::str_count( string = tree_original, pattern = "," ) + 1 )
phylip_header <- paste( num_sp, "  1", sep = "" )

## ---- NOT USED ANYMORE ---- ##
# write( x = phylip_header, file = paste( num_sp, "sp_", tname, "_32STcalib.tree", sep = "" ) )
# write( x = tree, file = paste( num_sp, "sp_", tname, "_32STcalib.tree", sep = "" ),
#        append = TRUE )
## ---- NOT USED ANYMORE ---- ##

write( x = phylip_header, file = paste( "01_trees_with_ST_calibs/Rout/",
                                        num_sp, "sp_", tname, "_71STcalib.tree", sep = "" ) )
write( x = tree_71ST, file = paste( "01_trees_with_ST_calibs/Rout/",
                                    num_sp, "sp_", tname, "_71STcalib.tree", sep = "" ),
       append = TRUE )

write( x = phylip_header, file = paste( "01_trees_with_ST_calibs/Rout/",
                                        num_sp, "sp_", tname, "_71STcalib_rounded.tree", sep = "" ) )
write( x = tree_71ST_rounded, file = paste( "01_trees_with_ST_calibs/Rout/",
                                            num_sp, "sp_", tname, "_71STcalib_rounded.tree", sep = "" ),
       append = TRUE )



## --- EXTRA 210812 | Only get tree with specific calibs --- ##

# USE 71 CALIBRATIONS ONLY 
# Replace calibration names with corresponding calibration
node.labels.ST.2 <- sort( as.numeric( rownames( ST.fit ) ), decreasing = TRUE )
check.indexes    <- sort( as.numeric( node.cal$Node.calib ), decreasing = TRUE )
for( i in node.labels.ST.2 ){
  
  # Find calibration and replace in tree with ST-fit
  if( i %in% check.indexes ){
    index.ST    <- which( i == as.numeric( rownames( ST.fit ) ) )
    tree_71ST_2 <- gsub( pattern = paste0(" ",i," "),
                         x = tree_71ST_2, replacement = paste( " '", ST.fit[index.ST,5], "' ", sep = "" ) )
  }else{ 
    tree_71ST_2 <- gsub( pattern = paste0(" ",i," "), x = tree_71ST_2, replacement = "" )
  }
  
  
  # Last step -- update sp labels
  if( i == min( node.cal$Node.calib ) ){
    num_sp <- 72:1
    for( j in num_sp ){
      tree_71ST_2 <- gsub( pattern = paste0( j,"_" ), x = tree_71ST_2, replacement = "" )
    }
  }
  
}


# USE 71 CALIBRATIONS ONLY - rounded
# Replace calibration names with corresponding calibration
node.labels.ST.2 <- sort( as.numeric( rownames( ST.fit ) ), decreasing = TRUE )
check.indexes    <- sort( as.numeric( node.cal$Node.calib ), decreasing = TRUE )
for( i in node.labels.ST.2 ){
  
  # Find calibration and replace in tree with ST-fit
  if( i %in% check.indexes ){
    index.ST <- which( i == as.numeric( rownames( ST.fit ) ) )
    tree_71ST_rounded_2 <- gsub( pattern = paste0(" ",i," "),
                               x = tree_71ST_rounded_2, replacement = paste( " '", ST.fit[index.ST,6], "' ", sep = "" ) )
  }else{
    tree_71ST_rounded_2 <- gsub( pattern = paste0(" ",i," "), x = tree_71ST_rounded_2, replacement = "" )
  }
  
  # Last step -- update sp labels
  if( i == min( node.cal$Node.calib ) ){
    num_sp <- 72:1
    for( j in num_sp ){
      tree_71ST_rounded_2 <- gsub( pattern = paste0( j,"_" ), x = tree_71ST_rounded_2, replacement = "" )
    }
  }
  
}

#-----------------#
# WRITE NEW FILES #
#-----------------#
num_sp <- c( stringr::str_count( string = tree_original, pattern = "," ) + 1 )
phylip_header <- paste( num_sp, "  1", sep = "" )

## ---- NOT USED ANYMORE ---- ##
# write( x = phylip_header, file = paste( num_sp, "sp_", tname, "_32STcalib.tree", sep = "" ) )
# write( x = tree, file = paste( num_sp, "sp_", tname, "_32STcalib.tree", sep = "" ),
#        append = TRUE )
## ---- NOT USED ANYMORE ---- ##

write( x = phylip_header, file = paste( "01_trees_with_ST_calibs/Rout/",
                                        num_sp, "sp_", tname, "_71STcalib_v2.tree", sep = "" ) )
write( x = tree_71ST_2, file = paste( "01_trees_with_ST_calibs/Rout/",
                                    num_sp, "sp_", tname, "_71STcalib_v2.tree", sep = "" ),
       append = TRUE )

write( x = phylip_header, file = paste( "01_trees_with_ST_calibs/Rout/",
                                        num_sp, "sp_", tname, "_71STcalib_rounded_v2.tree", sep = "" ) )
write( x = tree_71ST_rounded_2, file = paste( "01_trees_with_ST_calibs/Rout/",
                                            num_sp, "sp_", tname, "_71STcalib_rounded_v2.tree", sep = "" ),
       append = TRUE )
