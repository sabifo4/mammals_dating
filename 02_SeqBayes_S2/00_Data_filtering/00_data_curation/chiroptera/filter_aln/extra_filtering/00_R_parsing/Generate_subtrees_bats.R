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

#-------#
# TASKS #
#-------#
# 1. Load bats
bats_tree <- ape::read.tree( "main_bats_uncalib.tree" )

# 2. Chop outgroup from subtree
tmp.chop   <- bats_tree
which( bats_tree$tip.label == "nycteris_javanica" ) # 630
##> NOTE: keep the outgroup there
subt1         <- ape::drop.tip( tmp.chop, tip = c(1:630) )
ape::write.tree( subt1, file = "ChiroSubt1_megachiroptera.tree")

# 3. Generate now second subtree
subt2         <- ape::drop.tip( tmp.chop, tip = c(631:885) )
ape::write.tree( subt2, file = "ChiroSubt2_microchiroptera.tree")

# 4. Read the calibrated tree now and count species. Output a text 
#    file with the list of taxa included here
length(subt1$tip.label) # 255
writeLines( text = subt1$tip.label, con = "255sp_Lchiro_megachiroptera_list.txt" )

# 5. Repeat the same with the second subset
length(subt2$tip.label) # 630
labels_with_outg <- c(subt2$tip.label, "zaglossus_bruijni", "tachyglossus_aculeatus", "ornithorhynchus_anatinus" ) #633
writeLines( text = labels_with_outg, con = "633sp_Lchiro_microchiroptera_list.txt" )

##>> NOTE: This was an exploratory analysis. 
##   After knowing at which part the tree had to be divided into two 
##   and which taxa were to be kept in each data subset. We did not 
##   generate the calibrated trees until the extra taxa were 
##   added.