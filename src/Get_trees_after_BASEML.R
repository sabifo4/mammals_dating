# Function to parse all the BASEML trees with calibrations 
# manually added. Expected to run in a for loop so in 
# each iteration one tree is parsed.
#
# Arguments:
#    tt               character, vector with the path to the 
#                     tree inside "rooted_calibs".
#    calibs           character, vector with the path to the text file  
#                     with the calibrations to be included in output trees.
#    tname            character, vector with the name of the subtree that
#                     will be used to find the path where dummy alns are 
#                     to be saved.
#    name.tt          character, vector with the tree file name.
#    path.dummy.alns  character, path to the directory where dummy alns 
#                     are to be saved. Needs "/" at the end of the vector.
#    out              character, path to the directory where the output
#                     calibrated trees will be saved. Needs to include "/"
#                     at the end of the vector.
calib_tree <- function( tt, calibs, tname, name.tt, path.dummy.alns, out, num.part ){
  
  #=================#
  # 0. SCAN TREES   #
  #=================#
  cat( "Scannig tree ", tt, "and calibrations file ... ...\n" )
  # Get trees 
  tree_original <- readLines( con = tt )
  tree <- tree_original
  
  # Get sp names 
  spnames <- ape::read.tree( tt )$tip.label
  
  #======================================#
  #    1. GET FILE WITH CALIBRATIONS 
  #======================================#
  # Read tree and rename cols
  calibrations <- read.table( file = calibs,
                              stringsAsFactors = F, sep = "|", blank.lines.skip = T )
  colnames( calibrations ) <- c( "name", "calib" )
  
  #======================================#
  #  2. REPLACE NAMES WITH CALIBRATIONS
  #======================================#
  cat( "Generating calibrated tree ... ...\n" )
  # Replace calibration names with corresponding calibration
  for( i in 1:dim( calibrations )[1] ){
    
    # If there is an original calibration, use the original
    tree <- gsub( pattern = paste0("\\'",calibrations[i,1],"\\'"),
                  x = tree,
                  replacement = paste( "'", calibrations[i,2], "'", sep = "" ) )
  }
  
  #=======================#
  #  3. WRITE NEW FILES ! #
  #=======================#
  num_sp_count <- c( stringr::str_count( string = tree_original, pattern = "," ) + 1 )
  num_sp <- length( ape::read.tree( tt )$tip.label )
  all.equal( num_sp_count, num_sp )
  phylip_header <- paste( num_sp, "  1", sep = "" )
  
  phylip_header_aln <- paste( num_sp, "  2\n", sep = "" )
  spnames.2nuc      <- paste( spnames, "     AT", sep = "" )
  
  write( x = phylip_header, file = paste( out, num_sp, "sp_", tname,
                                          "_MCMCtree_calib_", num.part, ".tree", sep = "" ) )
  write( x = tree, file = paste( out, num_sp, "sp_", tname, "_MCMCtree_calib_", num.part,
                                 ".tree", sep = "" ),
         append = TRUE )
  write( x = spnames, file = paste( out, num_sp, "sp_", tname, "_spnameslist_", num.part,
                                    ".txt", sep = "" ) )
  
  # Generate dummy aln
  cat( "Generatng dummy alignment ... ...\n" )
  parts     <- "parts"
  if( ! dir.exists( paste( path.dummy.alns, parts, sep = "" ) ) ){
    dir.create( paste( path.dummy.alns, parts, sep = "" ) )
  }
  write( x = phylip_header_aln, file = paste( path.dummy.alns, parts,
                                              "/dummy_aln_", num.part, ".aln", sep = "" ) )
  write( x = spnames.2nuc, file = paste( path.dummy.alns, parts,
                                         "/dummy_aln_", num.part, ".aln", sep = "" ),
         append = TRUE )
  
  cat( "All tasks finished for tree ", tt, "!\n\n" )
}
