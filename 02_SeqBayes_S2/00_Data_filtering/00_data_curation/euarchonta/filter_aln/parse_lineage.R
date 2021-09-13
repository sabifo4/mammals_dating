#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-------------------------------------#
# AUTOMATICALLY SET WORKING DIRECTORY #
#-------------------------------------#
# Load package to help find the path to this source file 
library(rstudioapi) 
# Get the path to current open R script and find main dir
path_to_file <- getActiveDocumentContext()$path
wd <- gsub( pattern = "/R/scripts", replacement = "/", x = dirname( path_to_file ) )
# Set wd
setwd( wd )
# Load function
source( "../../../../../src/Filter_lineages.R" )

#-----------------#
# LOAD GENES FILE #
#-----------------#
genes <- readLines( "../../genes.txt" )

#-------------------#
# LOAD LINEAGE FILE #
#-------------------#
lineage.raw <- readLines( "lineage.txt" )

# Order
lineage.ord.filt1 <- grep( pattern = "order", x = lineage.raw, value = T )
lineage.ord.filt2 <- gsub( pattern = "..*order/", replacement = "order/", x = lineage.ord.filt1 )
# Family
lineage.fam.filt1 <- grep( pattern = "/family/", x = lineage.raw, value = T )
lineage.fam.filt2 <- gsub( pattern = "..*/family/", replacement = "family/", x = lineage.fam.filt1 )
# Subfamily
lineage.subf.filt1 <- grep( pattern = "subfamily/", x = lineage.raw, value = T )
lineage.subf.filt2 <- gsub( pattern = "..*subfamily/", replacement = "subfamily/", x = lineage.subf.filt1 )
# Genus
lineage.gen.filt1 <- grep( pattern = "/genus/", x = lineage.raw, value = T )
lineage.gen.filt2 <- gsub( pattern = "..*genus/", replacement = "genus/", x = lineage.gen.filt1 )

# Create vector with all lineages to check
check.lineages <- vector( mode = "list", length = 4 )
check.lineages[[ 1 ]] <- lineage.ord.filt2 
check.lineages[[ 2 ]] <- lineage.fam.filt2 
check.lineages[[ 3 ]] <- lineage.subf.filt2
check.lineages[[ 4 ]] <- lineage.gen.filt2
names( check.lineages ) <- c( "order", "family", "subfamily", "genus" )

# Now run the function for each of the lineages
levels.checked <- vector( mode = "list", 4 )
names( levels.checked ) <- names( check.lineages )
for( i in 1:length( check.lineages ) ){
  
  tmp.pattern <- paste( "~..*|", tolower( names( check.lineages )[i] ), "/", sep = "" ) # E.g. "~..*|order/"
  lins        <- gsub( pattern = tmp.pattern, replacement = "", x = check.lineages[[ i ]] )
  unique.lin  <- unique( lins )
  name.lin    <- names( check.lineages )[i]
  sp          <- gsub( pattern = "..*species/|[A-Z][a-z].*\\|", replacement = "", x = check.lineages[[ i ]] )
  
  levels.checked[[ i ]] <- filter.lineages( unique.lin = unique.lin, sp = sp,
                                            genes = genes, name.lin = name.lin )
  
}

save( levels.checked, file = "levels.checked.RData" )
# load( "levels.checked.RData" )
