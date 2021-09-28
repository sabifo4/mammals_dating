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
# 1. Load rodents tree with 1,314 sp
rodents_tree <- ape::read.tree( "main_rodents_uncalib.tree" )

# 2. Chop outgroup from subtree
tmp.chop   <- rodents_tree
which( rodents_tree$tip.label == "neusticomys_monticolus" ) # 685
##> NOTE: keep the outgroup there
subt1         <- ape::drop.tip( tmp.chop, tip = c(1:686) )
ape::write.tree( subt1, file = "Muridae.tree")

# 3. Generate now second subtree
subt2         <- ape::drop.tip( tmp.chop, tip = c(687:1311) )
ape::write.tree( subt2, file = "Eumoroida_rest.tree")

# 4. Read the calibrated tree now and count species. Output a text 
#    file with the list of taxa included here
length(subt1$tip.label) # 628
writeLines( text = subt1$tip.label, con = "628sp_muridae_list.txt" )

# 5. Repeat the same with the second subset
length(subt2$tip.label) # 689
writeLines( text = subt2$tip.label, con = "689sp_eumoroidarest_list.txt" )

##>> NOTE: This was an exploratory analysis. 
##   After knowing at which part the tree had to be divided into two 
##   and which taxa were to be kept in each data subset, we 
##   then manually generated the calibrated trees.
##   The trees were saved as `628sp_Muridae_calib_MCMCtree.tree`
##   and `689sp_Eumoroida_calib_MCMCtree.tree`.
##>> END