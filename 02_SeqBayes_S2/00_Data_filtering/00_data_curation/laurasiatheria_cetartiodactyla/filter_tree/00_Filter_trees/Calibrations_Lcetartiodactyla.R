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
wd2 <- gsub( pattern = "00_data_curation..*", replacement = "",
             x = wd )
# Set name 
tname <- "laurasiatheria_cetartiodactyla"

#=================#
# 2. SCAN TREES   #
#=================#
# Get trees 
tree_original <- readLines( con = "laurasiatheria_cetartiodactyla_calibnames.tree" )
tree <- tree_original

# Get sp names 
spnames <- ape::read.tree( "laurasiatheria_cetartiodactyla_calibnames.tree" )$tip.label

#======================================#
#    3. GET FILE WITH CALIBRATIONS 
#======================================#
# Read tree and rename cols
calibrations <- read.table( file = "Calibrations_LaurCetart.txt",
                            stringsAsFactors = F, sep = "|", blank.lines.skip = T )
colnames( calibrations ) <- c( "name", "calib" )

#======================================#
#  4. REPLACE NAMES WITH CALIBRATIONS
#======================================#
# Replace calibration names with corresponding calibration
for( i in 1:dim( calibrations )[1] ){
  
  # If there is an original calibration, use the original
  tree <- gsub( pattern = paste0("\\'",calibrations[i,1],"\\'"),
                x = tree,
                replacement = paste( "'", calibrations[i,2], "'", sep = "" ) )
}

#=======================#
#  5. WRITE NEW FILES ! #
#=======================#
num_sp_count <- c( stringr::str_count( string = tree_original, pattern = "," ) + 1 )
num_sp <- length( ape::read.tree( "laurasiatheria_cetartiodactyla_calibnames.tree" )$tip.label )
all.equal( num_sp_count, num_sp )
phylip_header <- paste( num_sp, "  1", sep = "" )

phylip_header_aln <- paste( num_sp, "  2\n", sep = "" )
spnames.2nuc      <- paste( spnames, "     AT", sep = "" )

write( x = phylip_header, file = paste( num_sp, "sp_", tname, "_MCMCtree_calib.tree", sep = "" ) )
write( x = tree, file = paste( num_sp, "sp_", tname, "_MCMCtree_calib.tree", sep = "" ),
       append = TRUE )
write( x = spnames, file = paste( num_sp, "sp_", tname, "_spnameslist.txt", sep = "" ) )


# Generate dummy aln
if( ! dir.exists( paste( wd2, "01_alignments/01_mammal_dummy_alns/", sep = "" ) ) ){
  dir.create( paste( wd2, "01_alignments/01_mammal_dummy_alns/", sep = "" ) )
}
if( ! dir.exists( paste( wd2, "01_alignments/01_mammal_dummy_alns/",
                         tolower( tname ), sep = "" ) ) ){
  dir.create( paste( wd2, "01_alignments/01_mammal_dummy_alns/",
                     tolower( tname ), sep = "" ) )
}
write( x = phylip_header_aln, file = paste( wd2, "01_alignments/01_mammal_dummy_alns/", 
                                            tolower( tname ), "/dummy_aln.aln", sep = "" ) )
write( x = spnames.2nuc, file = paste( wd2, "01_alignments/01_mammal_dummy_alns/", 
                                       tolower( tname ), "/dummy_aln.aln", sep = "" ), append = TRUE )

