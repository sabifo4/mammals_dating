#----------------------#
# 0. CLEAN ENVIRONMENT #
#----------------------#
rm( list = ls( ) )

#-------------------#
# 1. LOAD LIBRARIES #
#-------------------#
library( ape )

#--------------------------#
# 2. SET WORKING DIRECTORY #
#--------------------------#
library( rstudioapi ) 
# Get the path to current open R script and find main dir "00_Gene_filtering"
path_to_file <- getActiveDocumentContext()$path
script_wd <- paste( dirname( path_to_file ), "/", sep = "" )
wd <- gsub( pattern = "/scripts", replacement = "", x = script_wd )
setwd( wd )

#-----------#
# 3. TASKS  #
#-----------#
## 3.1. Read trees from "filtered_genes" and save them in a list
trees.dir <- paste( wd, "filtered_genes/", sep = "" )
num.genes <- length( list.files( trees.dir ) )

mammal.trees <- vector( mode = "list", num.genes )
names( mammal.trees ) <- list.files( trees.dir )

count <- 0
for ( i in 1:num.genes ){
  
  files <- list.files( paste( trees.dir,
                              list.files( trees.dir )[i],
                              sep = "" ) )
  check.files <-  grep( pattern = "bestTree", files )
  
  if( length( check.files ) != 0 ){
    
    count <- count + 1
    mammal.trees[[ i ]] <- ape::read.tree( file = paste( trees.dir,
                                                         list.files( trees.dir )[i],
                                                         "/RAxML_bestTree.",
                                                         list.files( trees.dir )[i],
                                                         sep = "" ) )
    # Ensure the name of the entry in the list is the same gene that has been added
    # tho this entry
    names( mammal.trees )[i] <- list.files( trees.dir )[i]
  }
  
  cat( count, "\n" )
  
}

# Save RData objects with filtered genes
save( mammal.trees, file = "out_RData/mammal.trees.RData" )

# Load the RData objects with filtered genes if you 
# are not generating the data using the command above
load( "out_RData/mammal.trees.RData" )

## 3.2. Sum branch lengths for each tree to get the tree length 
blengths  <- as.data.frame( matrix( 0, nrow = length( mammal.trees ), ncol = 5 ) )
rownames( blengths ) <- names( mammal.trees )
colnames( blengths ) <- c( "Sum.blengths", "Largbr.rel2tree",
                           "Tip.name", "H-M.dist",
                           "Delete?" )

for( i in 1:length( mammal.trees ) ){

  blengths[i,1] <- sum( mammal.trees[[i]]$edge.length )

}

# Manually checking that the tree length saved in column 1 in object
# `blengths` is the same than output by `RAxML`
#
# >> CHECK 1: gene ENSG00000000003 in RAxML = RAxML: 2.992765
ind.check1 <- which( rownames( blengths ) == "ENSG00000000003" )
blengths[ind.check1,1] # R = 2.992765
# >> SUCCESS!
#
# >> CHECK 2: gene ENSG00000006432 in RAxML = 5.895389
ind.check2 <- which( rownames( blengths ) == "ENSG00000006432" )
blengths[ind.check2,1] # R = 5.895389
# >> SUCCESS!


#-------------------------------------#
# 4. CALCULATE RELATIVE BRANCH LENGTH #
#-------------------------------------#
# Take each branch length (i.e., mammal.trees[[i]]$edge.length)
# and divide it into the corresponding tree length that was stored
# in column 1 for each gene (i.e., blenghts[i,1]). These will result 
# into a vector of relative branch lengths which, if multiplied
# by 100, you can use a useful statistic to detect possible
# alignment errors. If something went very
# very wrong (threshold set to relative branch lengths >= 60%),
# then the alignment makes no sense and the tree has
# a very incredible long branch for that species.
counter <- 0
for( i in 1:length( mammal.trees ) ){
  
  # Divide all branches into the tree length so we can obtain the relative 
  # branch lengths, i.e., r_{ij} = b{ij} / SUM_j (b_{ij}); as detailed in
  # dos Reis et al. 2012.
  # Then, take the largest relative branch length (in %) for each gene and store 
  # it in second  column in data.frame `blengths`.
  blengths[i,2] <- max( mammal.trees[[i]]$edge.length/blengths[i,1] )*100
  print( max( mammal.trees[[i]]$edge.length/blengths[i,1] )*100 )
  
  
  # If the largest relative branch length is equal to or larger than 60%...
  if ( round( blengths[i,2] ) >= 60 ){
    cat( "Gene ", rownames( blengths)[i], " - row", i,
         "- has a branch length longer than 60% of the total tree length\n" )
    write( paste( "Gene ", rownames( blengths)[i], " - row", i,
                  "- has a branch length longer than 60% of the total tree length: ",
                  blengths[i,2], "%\n", sep = ""),
           file = "out_logs/log_04_R_genes_blenghth_longer_60pc.txt",
           append = T )
    counter <- counter + 1
    blengths[i,5] <- c("Y")
    
    # 1. Find which position has largest branch length
    ind1   <- which( round( mammal.trees[[i]]$edge.length/blengths[i,1] ) >= 0.6 )
    # 2. In case more than one branch length was larger than 0.6 and stored as
    #   `ind1`, find the largest and save it as `ind1.1`
    ind1.1 <- max( c(mammal.trees[[i]]$edge.length/blengths[i,1])[ind1] )
    # 3. Find the position of this largest branch which index is `ind1.1`, just
    #    in case there were more than one saved in step 1
    ind1.2 <- which( mammal.trees[[i]]$edge.length/blengths[i,1] == ind1.1 )
    # 4. Find the node position of the individual which branch length is `ind1.2`
    ind2   <- mammal.trees[[i]]$edge[ind1.2,2]
    blengths[i,3] <- mammal.trees[[i]]$tip.label[ind2]
    
    # 5. Now we want to find the corresponding taxon name. But if
    #    ind2 is larger than the amount of taxa available in
    #    this gene tree, then increase ind1.2 by 1 until
    #    it finds the correct taxon name.
    if( ind2 >= length( mammal.trees[[i]]$tip.label ) ){
      while( ind2 >= length( mammal.trees[[i]]$tip.label ) ){
        ind1.2 <- ind1.2 + 1
        ind2   <- mammal.trees[[i]]$edge[ind1.2,2]
      }
      blengths[i,3] <- paste( mammal.trees[[i]]$tip.label[ind2], ",",
                              mammal.trees[[i]]$tip.label[ind2+1], sep = "" )
    }else{
      blengths[i,3] <- mammal.trees[[i]]$tip.label[ind2]
    }
  }else{
    blengths[i,3] <- "NULL"
    blengths[i,5] <- c("N")
  }
  
  # Append that the analysis has finished in the log file!
  if( i == length( mammal.trees ) ){
    cat("\nA total of ", counter, "genes had a branch length longer than 60% of
        the total tree length\n")
    write( paste( "\nA total of ", counter, "genes had a branch length longer than 60% of
                  the total tree length\n", sep = "" ),
           file = "out_logs/log_04_R_genes_blenghth_longer_60pc.txt",
           append = T )
  }
  
}

# Now, identify the 133 genes identified to have at least one branch length 
# longer than 60% of the total tree length.
# This will be later erased, but not now.
ind.delete <- which( blengths[,5] == "Y" )

# Save RData objects with filtered genes
save( blengths, file = "out_RData/mammals_summary_matrix_filtstep1.RData" )

# Load the RData objects with filtered genes if you 
# are not generating the data using the command above
load( "out_RData/mammals_summary_matrix_filtstep1.RData" )
## NOTE: If you are loading the data, please uncomment the
## next command so you can update `ind.delete` so it has the 
## indexes for the 133 genes identified to have at least 
## one branch length longer than 60% of the total tree length
# ind.delete <- which( blengths[,5] == "Y" )

#------------------------------------------------------#
# 5. Calculating pairwise distance from mouse to human #
#------------------------------------------------------#
## 5.1. Read *_mouse_human.aln files in "baseml" and save them in a list
baseml.dir <- c( "baseml/" )

MH_aln.aln <- vector( mode = "list", num.genes )
names( MH_aln.aln ) <- list.files( trees.dir )

count <- 0
for ( i in 1:num.genes ){
  
  files <- list.files( paste( baseml.dir, i, sep = "" ) )
  check.files <-  grep( pattern = "mouse_human", files )
  
  if( length( check.files ) != 0 ){
    
    count <- count + 1
    # Read the file
    aln <- read.table( file = paste( baseml.dir, i,
                                     "/", files[check.files], sep = "" ),
                       sep = "\t", stringsAsFactors = FALSE )
    # Get species names
    sp_names <- aln[,1]
    
    # Get length of sequences and construct the matrix that will be later
    # used to calculate the H-M distance
    num_seqs <- length( unlist( strsplit(aln[1,2], split = "" ) ) )
    aln_mat  <- matrix( 0, nrow = dim( aln )[1], ncol = num_seqs )
    rownames( aln_mat ) <- sp_names
    
    # Append to matrix
    for( j in seq( 1:dim(aln)[1] ) ){
      
      aln_mat[j,] <- unlist( strsplit(aln[j,2], split = "" ) )
      
    }
    
    # Put in the list
    MH_aln.aln[[ i ]] <- ape::as.DNAbin( aln_mat )
    
  }
  
  cat( count, "\n" )
  
}

## 5.2. Save RData objects with MH alignments as DNAbin format.
##      If you have already generated this file, you can run the 
##      second command and load the corresponding RData file
save( MH_aln.aln, file = "out_RData/MHalns_DNAbin.RData" )
load( "out_RData/MHalns_DNAbin.RData" )

## 5.3. Go through every matrix in the MH_aln.aln list 
##      and compute the MH distance
## We picked the TN93 model because this is the one that BASEML uses:
## Page 26 PAML doc:
## ** The pairwise sequence distances are included in the output as well, and also in
## ** a separate file called 2base.t. This is a lower-diagonal distance matrice, readable
## ** by the NEIGHBOR program in Felesenstein's PHYLIP package (Felsenstein 2005).
## ** For models JC69, K80, F81, F84, the appropriate distance formulas are used,
## ** while for more complex models, the TN93 formula is used. BASEML is mainly a maximum
## ** likelihood program, and the distance matrix is printed out for convenience and really
## ** has nothing to do with the later likelihood calculation.

# Use model TN93
count <- 0
for( i in 1:length( MH_aln.aln ) ){
  
  count <- count + 1
  MH.dist <- ape::dist.dna( MH_aln.aln[[ i ]], model = "TN93" )
  blengths[i,4] <- MH.dist
  
  # Safety check
  files <- list.files( paste( baseml.dir, i, sep = "" ) )
  check.files <- grep( pattern = "mouse_human", files )
  cat( "Gene as name in list: ", names( MH_aln.aln )[i], "\n",
       "Gene as name in order it was appended: ", files[check.files], "\n",
       "Gene as name in matrix blengths: ", rownames( blengths )[i], "\n" )
  cat( count, "\n" )
  
  write( paste( "Gene as name in list: ", names( MH_aln.aln )[i], "\n",
                "Gene as name in order it was appended: ", files[check.files], "\n",
                "Gene as name in matrix blengths: ", rownames( blengths )[i], "\n",
                sep = "" ),
         file = "out_logs/log_05_R_check_genes_included_MH_aln_TN93.txt",
         append = TRUE )
  
}

# Do the same with JC69 and raw models. So extend the data frame two more columns.
# This is why the object `blengths2` is created here.
blengths2       <- as.data.frame( matrix( 0, nrow = length( mammal.trees ), ncol = 6 ) )
blengths2[,1:4] <- blengths[,1:4]
rownames( blengths2 ) <- names( mammal.trees )
colnames( blengths2 ) <- c( "Sum.blengths", "Largbr.rel2tree", "Tip.name",
                            "H-M.dist.TN93", "H-M.dist.JC69", "H-M.dist.raw" )

count <- 0
for( i in 1:length( MH_aln.aln ) ){
  
  count <- count + 1
  MH.dist.JC69   <- ape::dist.dna( MH_aln.aln[[ i ]], model = "JC69" )
  blengths2[i,5] <- MH.dist.JC69
  MH.dist.raw    <- ape::dist.dna( MH_aln.aln[[ i ]], model = "raw" )
  blengths2[i,6] <- MH.dist.raw
  
  # Safety check
  files <- list.files( paste( baseml.dir, i, sep = "" ) )
  check.files <-  grep( pattern = "mouse_human", files )
  cat( "Gene as name in list: ", names( MH_aln.aln )[i], "\n",
       "Gene as name in order it was appended: ", files[check.files], "\n",
       "Gene as name in matrix blengths2: ", rownames( blengths2 )[i], "\n" )
  cat( count, "\n" )
  write( paste( "Gene as name in list: ", names( MH_aln.aln )[i], "\n",
                "Gene as name in order it was appended: ", files[check.files], "\n",
                "Gene as name in matrix blengths2: ", rownames( blengths2 )[i], "\n",
                sep = "" ),
         file = "out_logs/log_05_R_check_genes_included_MH_aln_JC69_raw.txt",
         append = TRUE )
  
}

# Save the object so you can avoid running the loop above
save( blengths2,
      file = "out_RData/mammals_summary_matrix.TN93.JC69.raw_filtstep2.RData" )
# Load the file if you have already generated it
load( "out_RData/mammals_summary_matrix.TN93.JC69.raw_filtstep2.RData" )

#------------------------------------------------------------#
# 6. Find those genes which the largest single branch is     #
#    larger than or equal to 60% relative to the tree length #
#------------------------------------------------------------#
# Find index of largest bl 
ind.larg.b  <- which( round(blengths2$Largbr.rel2tree) >= 60 )
# >> Check that we have 133 and that they are the same than found before
length( ind.larg.b ) # 133
setdiff( ind.larg.b, ind.delete ) # 0
# >> SUCCESS! 
filt.blengths <- blengths2[-c( ind.larg.b ),]
# >> Safety check: Double check that now all rows in col 3 are "NULL"
length( which( filt.blengths[,3] == "NULL" ) ) == dim( filt.blengths )[1]
# Double check and OK!
filt.df.blengths <- filt.blengths[,-3] # Remove 3rd column as double check done!
# >> SUCCESS!

# Which are the tips which bl is too large?
large.blengths <- blengths2[ind.larg.b,]
tips.longb <- factor( blengths2[ind.larg.b,3] ) #24 different
#plot( tips.longb, las = 3 )

#-------------------------------------------#
# 7. Order filtered genes from slow to fast #
#-------------------------------------------#
# Create objects with ordered genes (fast- to slow-evolving genes) according to 
# the three models tested
ord.blengths.f2s.TN93 <- filt.df.blengths[order( round( filt.df.blengths$`H-M.dist.TN93`, 
                                                        15) ),]
ord.blengths.f2s.JC69 <- filt.df.blengths[order( round( filt.df.blengths$`H-M.dist.JC69`, 
                                                        15) ),]
ord.blengths.f2s.raw <- filt.df.blengths[order( round( filt.df.blengths$`H-M.dist.raw`, 
                                                       15) ),]

# Create objects with ordered genes (slow- to fast-evolving genes) according to 
# the three models tested
ord.blengths.s2f.TN93 <- filt.df.blengths[order( round( filt.df.blengths$`H-M.dist.TN93`, 
                                                        15 ),decreasing = TRUE ),]
ord.blengths.s2f.JC69 <- filt.df.blengths[order( round( filt.df.blengths$`H-M.dist.JC69`, 
                                                        15 ),decreasing = TRUE ),]
ord.blengths.s2f.raw <- filt.df.blengths[order( round( filt.df.blengths$`H-M.dist.raw`, 
                                                       15 ),decreasing = TRUE ),]

# >> CHECK | Which are NaN or Inf --> very different sequences.
#    Do this check in the three data frames previously created with ordered genes. 

# 1. Use object `ord.blengths.s2f.TN93`, which genes were ordered from slow- 
#    to fast-evolving according to the model `TN93`.
#> Column H-M.dist.TN93 (index = 3)
which( ord.blengths.s2f.TN93[,3] == "Inf" ) 
# integer(0)
which( ord.blengths.s2f.TN93[,3] == "NaN" ) 
#[1] 15435 15436
which( ord.blengths.s2f.TN93[,3] >= 0.75 )
#[1] 1 2
ord.blengths.s2f.TN93[c( 1, 2, 15435, 15436 ),]
rownames( ord.blengths.s2f.TN93 )[ c( 1, 2, 15435, 15436 )]
#[1] "ENSG00000120937" "ENSG00000204544" "ENSG00000132185" "ENSG00000236699"

#> Column H-M.dist.JC69 (index = 4)
which( ord.blengths.s2f.TN93[,4] == "Inf" )
# integer(0)
which( ord.blengths.s2f.TN93[,4] == "NaN" ) 
#[1] 15436
which( ord.blengths.s2f.TN93[,4] >= 0.75 )
#[1] 15435
ord.blengths.s2f.TN93[c( 15435, 15436 ),]
rownames( ord.blengths.s2f.TN93 )[ c( 15435, 15436 )]
#[1] "ENSG00000132185" "ENSG00000236699"

#> Column H-M.dist.raw (index = 5)
which( ord.blengths.s2f.TN93[,5] == "Inf" ) #column H-M.dist.raw
# integer(0)
which( ord.blengths.s2f.TN93[,5] == "NaN" ) 
#[1] 15436
which( ord.blengths.s2f.TN93[,5] >= 0.75 )
# integer(0)
ord.blengths.s2f.TN93[c( 15436 ),]
rownames( ord.blengths.s2f.TN93 )[ c( 15436 )]
#[1] "ENSG00000236699"


# 2. Use object `ord.blengths.s2f.JC69`, which genes were ordered from slow- 
#    to fast-evolving according to the model `JC69`.
# Column H-M.dist.TN93 (index = 3)
which( ord.blengths.s2f.JC69[,3] == "Inf" )
# integer(0)
which( ord.blengths.s2f.JC69[,3] == "NaN" ) 
#[1] 1 15436
which( ord.blengths.s2f.JC69[,3] >= 0.75 )
#[1] 4 8
ord.blengths.s2f.JC69[c( 1, 4, 8, 15436 ),]
rownames( ord.blengths.s2f.JC69 )[ c( 1, 4, 8, 15436 )]
#[1] "ENSG00000132185" "ENSG00000204544" "ENSG00000120937" "ENSG00000236699"

#> Column H-M.dist.JC69 (index = 4)
which( ord.blengths.s2f.JC69[,4] == "Inf" ) 
# integer(0)
which( ord.blengths.s2f.JC69[,4] == "NaN" ) 
#[1] 15436
which( ord.blengths.s2f.JC69[,4] >= 0.75 )
#[1] 1
ord.blengths.s2f.JC69[c( 1, 15436 ),]
rownames( ord.blengths.s2f.JC69 )[ c( 1, 15436 )]
#[1] "ENSG00000132185" "ENSG00000236699"

#> Column H-M.dist.raw (index = 5)
which( ord.blengths.s2f.JC69[,5] == "NaN" ) 
#[1] 15436
which( ord.blengths.s2f.JC69[,5] == "Inf" )
# integer(0)
which( ord.blengths.s2f.JC69[,5] >= 0.75 )
# integer(0)
ord.blengths.s2f.JC69[c( 15436 ),]
rownames( ord.blengths.s2f.JC69 )[ c( 15436 )]
#[1] "ENSG00000236699"

# 3. Use object `ord.blengths.s2f.raw`, which genes were ordered from slow-
#     to fast-evolving according to the model `raw`.
#> Column H-M.dist.TN93 (index = 3)
which( ord.blengths.s2f.raw[,3] == "Inf" )
# integer(0)
which( ord.blengths.s2f.raw[,3] == "NaN" ) 
#[1] 1 15436
which( ord.blengths.s2f.raw[,3] >= 0.75 )
#[1] 4 8
ord.blengths.s2f.raw[c( 1, 4, 8, 15436 ),]
rownames( ord.blengths.s2f.raw )[ c( 1, 4, 8, 15436 )]
#[1] "ENSG00000132185" "ENSG00000204544" "ENSG00000120937" "ENSG00000236699"

#> Column H-M.dist.JC69 (index = 4)
which( ord.blengths.s2f.raw[,4] == "NaN" ) #column H-M.dist.JC69
#[1] 15436
which( ord.blengths.s2f.raw[,4] == "Inf" ) #column H-M.dist.JC69
# integer(0)
which( ord.blengths.s2f.raw[,4] >= 0.75 )
#[1] 1
ord.blengths.s2f.raw[c( 1, 15436 ),]
rownames( ord.blengths.s2f.raw )[ c( 1, 15436 )]
#[1] "ENSG00000132185" "ENSG00000236699"

#> Column H-M.dist.raw (index = 5)
which( ord.blengths.s2f.raw[,5] == "Inf" ) #column H-M.dist.raw
# integer(0)
which( ord.blengths.s2f.raw[,5] == "NaN" ) #column H-M.dist.raw
#[1] 15436
which( ord.blengths.s2f.raw[,5] >= 0.75 )
# integer(0)
ord.blengths.s2f.raw[c( 15436 ),]
rownames( ord.blengths.s2f.raw )[ c( 15436 )]
#[1] "ENSG00000236699"

# >> END CHECK

# See the genes for which at least one of the models had issues to compute the distance
# between mouse and human sequences or this distance was equal or higher than 0.75.
# Genes: "ENSG00000132185" "ENSG00000204544" "ENSG00000120937" "ENSG00000236699"
ord.blengths.s2f.TN93[c( 1, 2, 15435, 15436 ),]

#                 Sum.blengths Largbr.rel2tree H-M.dist.TN93 H-M.dist.JC69 H-M.dist.raw
# ENSG00000120937     5.334116       12.569959     0.7603284     0.6553421    0.4369748
# ENSG00000204544     5.230688       14.173488     0.7556994     0.7018490    0.4557957
# ENSG00000132185     6.601777        8.205249           NaN     1.6479184    0.6666667
# ENSG00000236699     5.758692       54.273060           NaN           NaN          NaN

# Delete genes
ord.blengths.s2f.TN93.filt <- ord.blengths.s2f.TN93[ -c( 1, 2, 15435, 15436 ),]
dim( ord.blengths.s2f.TN93.filt ) # [1] 15432     5

# Write table with the matrix
write.table( ord.blengths.s2f.TN93.filt,
             file = "out_RData/filtered_ordered_matrix.csv",
             quote = F, sep = "," )

## 7.1. Plot filtered and not filtered matrix with relative and absolute branch lengths
# Output a pdf with the plot
pdf( file = "out_RData/check_relblVStreelength.pdf", paper = "a4" )
plot( ord.blengths.s2f.TN93.filt[,1], ord.blengths.s2f.TN93.filt[,2],
      main = "Tree length VS Relative branch lengths",
      xlab = "Tree length", ylab = "Relative branch lengths" )
dev.off()
pdf( file = "out_RData/check_logtreelengthVSlargestbl.pdf", paper = "a4" )
plot( log(ord.blengths.s2f.TN93.filt[,1]), ord.blengths.s2f.TN93.filt[,2],
      main = "Quality control check of relative branch lengths",
      xlab = "Total tree length (log-scale)", ylab = "Relative branch lengths (%)" )
dev.off()

# There seems to be an outlier and, by manually checking this alignment,
# it has two very long branch lenghts. Therefore, we remove this gene.
ind.outlier <- which( max(ord.blengths.s2f.TN93.filt[,1]) ==
                        ord.blengths.s2f.TN93.filt[,1] )
rownames( ord.blengths.s2f.TN93.filt )[ind.outlier] # [1] "ENSG00000176973"
ord.blengths.s2f.TN93.filt <- ord.blengths.s2f.TN93.filt[-ind.outlier,]

# Output a pdf with the plot
pdf( file = "out_RData/check_relblVStreelength_nooutlier.pdf", paper = "a4" )
plot( ord.blengths.s2f.TN93.filt[,1], ord.blengths.s2f.TN93.filt[,2],
      main = "Tree length VS Relative branch lengths",
      xlab = "Tree length", ylab = "Relative branch lengths",
      xlim = c(0, 60) )
dev.off()
pdf( file = "out_RData/check_logtree_relblVStreelength_nooutlier.pdf", paper = "a4" )
plot( log(ord.blengths.s2f.TN93.filt[,1]), ord.blengths.s2f.TN93.filt[,2],
      main = "Quality control check of relative branch lengths - no outlier",
      xlab = "Total tree length (log-scale)", ylab = "Relative branch lengths (%)" )
dev.off()

# Final size
dim( ord.blengths.s2f.TN93.filt )[1] # 15431

save( ord.blengths.s2f.TN93.filt,
      file = "out_RData/mammals_summary_matrix.TN93_filtstep3.RData" )

# Code used to plot the comparison between the filtered genes before and after
# removing the outlier
par( mfrow = c( 1, 2 ) )
plot( ord.blengths.s2f.TN93[ -c( 1, 2, 15435, 15436 ),1],
      ord.blengths.s2f.TN93[ -c( 1, 2, 15435, 15436 ),2],
      main = "Before removing outlier",
      xlab = "Tree length", ylab = "Relative branch lengths" )
plot( ord.blengths.s2f.TN93.filt[,1], ord.blengths.s2f.TN93.filt[,2],
      main = "After removing outlier",
      xlab = "Tree length", ylab = "Relative branch lengths",
      xlim = c(0, 60) )

# Load object
load( "out_RData/mammals_summary_matrix.TN93_filtstep3.RData" )

#--------------------------------------------#
# 8. Get a summary statistics of the filters #
#--------------------------------------------#
# Start counter with starting number of genes. Then proceed with the sum stats.
start.num.genes <- 15904
##============================================================================##
## UPDATE -- We added an extra filter to account for those genes that had to  ##
## be removed because they are present in the second data set too. This is    ##
## the new extra column. We know that there are a total of 15,268 genes that  ##
## are not overlapping with the second data set (explained in the README.md   ##
## file in the GitHub repository as bash scripts were used for that purpose). ##
## Therefore, the number in the last column is manualy added to `sum.stats`   ##
## as you can see in line 564 to compute the summary statistics.              ##
##============================================================================##
#sum.stats <- as.data.frame( matrix( 0, nrow = 5, ncol = 4 ) )
#colnames(sum.stats)   <- c( "Raw", "First_filter", "Second_filter", "Third filter" )
sum.stats <- as.data.frame( matrix( 0, nrow = 5, ncol = 5 ) )
colnames(sum.stats)   <- c( "Raw", "First_filter", "Second_filter", "Third filter", "Fourth filter" )
rownames( sum.stats ) <- c( "Num_genes", "Lost_genes", "%_Loss",
                            "Cum_lost_genes", "%_Cum_loss" )
sum.stats[1,] <- c( start.num.genes, dim( blengths)[1], dim(filt.df.blengths)[1],
                    dim(ord.blengths.s2f.TN93.filt)[1], 15268 )
sum.stats[2,] <- c( sum.stats[1,1]-sum.stats[1,1],
                    sum.stats[1,1]-sum.stats[1,2],
                    sum.stats[1,2]-sum.stats[1,3],
                    sum.stats[1,3]-sum.stats[1,4],
                    sum.stats[1,4]-sum.stats[1,5])
sum.stats[3,] <- c( round( (100 - (sum.stats[1,1]*100/sum.stats[1,1])),3 ),
                    round( (100 - (sum.stats[1,2]*100/sum.stats[1,1])),3 ),
                    round( (100 - (sum.stats[1,3]*100/sum.stats[1,2])),3 ),
                    round( (100 - (sum.stats[1,4]*100/sum.stats[1,3])),3 ),
                    round( (100 - (sum.stats[1,5]*100/sum.stats[1,4])),3 ))
sum.stats[4,] <- c( start.num.genes-sum.stats[1,1],
                    start.num.genes-sum.stats[1,2],
                    start.num.genes-sum.stats[1,3],
                    start.num.genes-sum.stats[1,4],
                    start.num.genes-sum.stats[1,5])
sum.stats[5,] <- c( round( (100 - (sum.stats[1,1]*100/start.num.genes)),3 ),
                    round( (100 - (sum.stats[1,2]*100/start.num.genes)),3 ),
                    round( (100 - (sum.stats[1,3]*100/start.num.genes)),3 ),
                    round( (100 - (sum.stats[1,4]*100/start.num.genes)),3 ),
                    round( (100 - (sum.stats[1,5]*100/start.num.genes)),3 ))

# Write summary table
#write.table(file = "out_RData/filtering_sum_stats.csv", sum.stats, quote = F, sep = ",")
write.table(file = "out_RData/filtering_sum_stats_nonucgenes.csv", sum.stats, quote = F, sep = ",")

#----------------------------------------------------------#
# 9. Copy filtered genes but saving them in an ordered way #
#    according to fast-slow evolving                       #
#----------------------------------------------------------#
# Create objects 
name_ord_filt_genes <- rownames( ord.blengths.s2f.TN93.filt )
name_notordered     <- names( MH_aln.aln )
list_correspondance <- rep(0, length( name_ord_filt_genes ) )
# Create a list with the indexes of the genes not ordered
# corresponding to the position in ordered list
for ( i in 1:length(name_ord_filt_genes) ){
  list_correspondance[i] <- which( name_notordered == name_ord_filt_genes[i] )
}
# >> CHECK | Double check everything alright !
all.equal( name_ord_filt_genes, name_notordered[list_correspondance] ) # TRUE!
# E.g. If i = 1
k <- 1
ind.k <- which( name_notordered == name_ord_filt_genes[k] )
all.equal( name_notordered[ind.k], name_ord_filt_genes[k] )
#[1] TRUE
# E.g.2 
name_notordered[list_correspondance[1:10]]
name_ord_filt_genes[1:10]
all.equal( name_notordered[list_correspondance[1:10]], name_ord_filt_genes[1:10] )
#[1] TRUE
# >> CHECK FINISHED 

# Generate the data
# Check if directory where to copy files exist. Otherwise, create it.
if ( ! dir.exists( "filtered_genes_step2/" ) ){
  dir.create( "filtered_genes_step2/" )
}
counter = 0 
for( i in seq( 1:length( name_ord_filt_genes ) ) ){
  
  counter = counter + 1
  # Check if directory where to copy files exist. Otherwise, create it.
  if ( ! dir.exists( paste( "filtered_genes_step2/", i, "/", sep = "") ) ){
    dir.create( paste( "filtered_genes_step2/", i, "/", sep = "" ) )
  }
  # Copy files
  file.copy( paste( "baseml/", list_correspondance[i], "/",
                    "partitions12_", name_notordered[ list_correspondance[i] ],
                    ".aln", sep = "" ),
             paste( "filtered_genes_step2/", i, "/", sep = "") )
  file.copy( paste( "baseml/", list_correspondance[i], "/",
                    name_notordered[ list_correspondance[i] ], ".tree", sep = "" ),
             paste( "filtered_genes_step2/", i, "/", sep = "") )
  cat( counter, "\n" )
  write( paste( "Gene ", counter, " is: ",
                name_notordered[ list_correspondance[i] ],
                ". Copy *aln and *tree from folder ",
                list_correspondance[i], " in baseml dir! ", "\n", sep = "" ),
         "out_logs/log_06_R_copy_step2_filtered_genes.txt",
         append = TRUE, sep = "\n" )
}

#-------------------------------------------------------------------------------#
# 10. Copy filtered genes which are present in all 72 species and in an ordered #
#     way according to slow-fast evolving                                       #
#-------------------------------------------------------------------------------#
## 10.1. Prase trees from "filtered_genes_all72sp" 
# Read trees from "filtered_genes_all72sp" and save them in a vector
trees72sp.dir <- c( "filtered_genes_all72sp/" )
num.genes72sp <- length( list.files( trees72sp.dir ) )

genes72sp.ind.ordered    <- which( name_ord_filt_genes %in% list.files( trees72sp.dir ) )
name_ord_filt_genes_72sp <- name_ord_filt_genes[ genes72sp.ind.ordered ]
list_correspondance_72sp <- rep(0, length( genes72sp.ind.ordered ) )
# Create a list with the corresponding indexes of the genes ordered
# but in the not ordered list.
# Now we have 645 genes as 3 might have been trimmed due to the previous
# filterings applied !
for ( i in 1:length( genes72sp.ind.ordered ) ){
  list_correspondance_72sp[i] <- which( name_notordered == name_ord_filt_genes_72sp[i] )
}

# >> CHECK 
k <- 1
which( name_notordered == name_ord_filt_genes_72sp[k] )
#[1] 11609
name_notordered[11609]
#[1] "ENSG00000172817"
genes72sp.ind.ordered[1]
#[1] 1237
name_ord_filt_genes[ genes72sp.ind.ordered[1] ]
#[1] "ENSG00000172817"
# >> END CHECK

# Get blengths ordered object with only genes present in 72 species
ord.blengths.s2f.TN93.72sp <- ord.blengths.s2f.TN93.filt[genes72sp.ind.ordered,]
dim( ord.blengths.s2f.TN93.72sp ) # [1] 645   5
# Write table with the matrix
write.table(file = "out_RData/filtered_ordered_matrix_all72sp.csv",
            ord.blengths.s2f.TN93.72sp, quote = F, sep = ",")

## 10.2. Generate the data
# Start counter
counter <- 0 
# Check if directory where to copy files exist. Otherwise, create it.
if ( ! dir.exists( "filtered_genes_step2_all72sp/" ) ){
  dir.create( "filtered_genes_step2_all72sp/"  )
}
for( i in seq( 1:length( name_ord_filt_genes_72sp ) ) ){
  
  counter = counter + 1
  # Check if directory where to copy files exist. Otherwise, create it.
  if ( ! dir.exists( paste( "filtered_genes_step2_all72sp/", i, "/", sep = "") ) ){
    dir.create( paste( "filtered_genes_step2_all72sp/", i, "/", sep = "" ) )
  }
  # Copy files
  file.copy( paste( "baseml/", list_correspondance_72sp[i], "/",
                    "partitions12_", name_notordered[ list_correspondance_72sp[i] ],
                    ".aln", sep = "" ),
             paste( "filtered_genes_step2_all72sp/", i, "/", sep = "") )
  file.copy( paste( "baseml/", list_correspondance_72sp[i], "/",
                    name_notordered[ list_correspondance_72sp[i] ], ".tree", sep = "" ),
             paste( "filtered_genes_step2_all72sp/", i, "/", sep = "") )
  cat( counter, "\n" )
  write( paste( "Gene ", counter, " is: ",
                name_notordered[ list_correspondance_72sp[i] ],
                ". Copy *aln and *tree from folder ",
                list_correspondance_72sp[i], " in baseml dir! ", "\n", sep = "" ),
         "out_logs/log_06_R_copy_step2_filtered_all72sp_genes.txt",
         append = TRUE, sep = "\n" )
}

#---------------------------------------------------------------------------#
# 11. Create subsets for the Bayesian model selection analysis by randomly  #
#     sampling the final filtered pool of genes.                            #
#     NOTE: This was done before we filtered out the nuclear genes that     #
#     overlap with the second data set (only 11). This does not affect the  # 
#     Bayesian model selection analysis as it has nothing to do with the    #
#     sequential Bayesian dating analysis (it is an independent analysis).  #
#---------------------------------------------------------------------------#
## 11.1 Define functions to generate data subsets
# FUNCTION 1
# Function that samples genes common in all 72 mammal taxa and saves them 
# in the specified directory by the user.
#
# Argumnets:
# range_num   Numeric, vector generated as a first step defining the position 
#             of the genes in `list_correspondance_72sp` that are to be sampled.
# name_dir    Character, name of the directory in which the sampled genes will 
#             be saved.
# out         Character, tag used to append in the log files to identify how 
#             many genes have been sampled. E.g., "30sp".
# extra       Boolean. FALSE when sampling data subsets with less than 645 genes 
#             as sampling will have happened from this gene pool. For data subsets 
#             with more genes, these 645 genes are assumed to be already included, 
#             so they are already copied. The rest of the genes were sampled from 
#             the big pool of filtered genes until reaching the number of genes to 
#             be sampled for that specific data subset (e.g., 1,000 genes,
#             10,000 genes).
# only_report Boolean, default is FALSE. If you do not want to copy the genes
#             (e.g., you have already done that) and you just want a summary output 
#             then set to TRUE. Otherwise, set to FALSE.
subsample <- function( range_num, name_dir, out, extra = FALSE, only_report = FALSE ){
  
  if( extra == FALSE ){
    counter <- 0
    for( i in range_num ){
      counter <- counter + 1
      # Check if directory where to copy files exist. Otherwise, create it
      if ( ! dir.exists( paste( name_dir, "/", counter, "/", sep = "") ) ){
        dir.create( paste( name_dir, "/", counter, "/", sep = "" ) )
      }
      
      # Copy files
      name_gene <- gsub( x = list.files( paste( "filtered_genes_step2_all72sp/", i, "/", sep = "" ),
                                         pattern = "tree" ),
                         pattern = ".tree", replacement = "" )
      if( only_report == TRUE ){
        write( name_gene, paste( name_dir, "/log_R_subsample", out, "_genes.txt", sep = "" ),
               append = TRUE, sep = "\n" )
      }else if( only_report == FALSE ){
        cp_files <- list.files( paste( "filtered_genes_step2_all72sp/", i, "/", sep = "" ), 
                                pattern = "tree|aln", recursive = TRUE )
        file.copy( paste( "filtered_genes_step2_all72sp/", i, "/", cp_files[1], sep = "" ),
                   paste( name_dir, "/", counter, "/", sep = ""), recursive = TRUE )
        file.copy( paste( "filtered_genes_step2_all72sp/", i, "/", cp_files[2], sep = "" ),
                   paste( name_dir, "/", counter, "/", sep = ""), recursive = TRUE )
        write( paste( "Gene ", counter, " -- ", i, " :", name_gene, sep = "" ),
               paste( name_dir, "/log_R_subsample", out, ".txt", sep = "" ),
               append = TRUE, sep = "\n" )
      }
      
    }
  }else if( extra == TRUE ){
    
    count <- 645
    for( i in range_num ){
      count <- count + 1
      # Check if directory where to copy files exist. Otherwise, create it
      if ( ! dir.exists( paste( name_dir, "/", count, "/", sep = "") ) ){
        dir.create( paste( name_dir, "/", count, "/", sep = "" ) )
      }
      
      # Copy files
      name_gene <- gsub( x = list.files( paste( "filtered_genes_step2/", i, "/", sep = "" ),
                                         pattern = "tree" ),
                         pattern = ".tree", replacement = "" )
      if( only_report == TRUE ){
        write( name_gene, paste( name_dir, "/log_R_subsample", out, "_genes.txt", sep = "" ),
               append = TRUE, sep = "\n" )
      }else if( only_report == FALSE ){
        cp_files <- list.files( paste( "filtered_genes_step2/", i, "/", sep = "" ), 
                                pattern = "tree|aln", recursive = TRUE )
        file.copy( paste( "filtered_genes_step2/", i, "/", cp_files[1], sep = "" ),
                   paste( name_dir, "/", count, "/", sep = ""), recursive = TRUE )
        file.copy( paste( "filtered_genes_step2/", i, "/", cp_files[2], sep = "" ),
                   paste( name_dir, "/", count, "/", sep = ""), recursive = TRUE )
        write( paste( "Gene ", count, " -- ", i, " :", name_gene, sep = "" ),
               paste( name_dir, "/log_R_subsample", out, ".txt", sep = "" ),
               append = TRUE, sep = "\n" )
      }
    }
    
  }
  
}

# FUNCTION 2
# Function that samples genes from the big pool of genes that have been filtered
# and saves them in the specified directory by the user.
#
# Argumnets:
# range_num   Numeric, vector generated as a first step defining the position 
#             of the genes in `list_correspondance_72sp` that are to be sampled.
# name_dir    Character, name of the directory in which the sampled genes will 
#             be saved.
# out         Character, tag used to append in the log files to identify how 
#             many genes have been sampled. E.g., "30sp".
# only_report Boolean, default is FALSE. If you do not want to copy the genes
#             (e.g., you have already done that) and you just want a summary output 
#             then set to TRUE. Otherwise, set to FALSE.
subsample_fromall <- function( range_num, name_dir, out, only_report = FALSE ){
  
  counter <- 0
  for( i in range_num ){
    counter <- counter + 1
    # Check if directory where to copy files exist. Otherwise, create it
    if ( ! dir.exists( paste( name_dir, "/", counter, "/", sep = "") ) ){
      dir.create( paste( name_dir, "/", counter, "/", sep = "" ) )
    }
    
    # Copy files
    name_gene <- gsub( x = list.files( paste( "filtered_genes_step2/", i, "/", sep = "" ),
                                       pattern = "tree" ),
                       pattern = ".tree", replacement = "" )
    if( only_report == TRUE ){
      write( name_gene, paste( name_dir, "/log_R_subsample", out, "_genes.txt", sep = "" ),
             append = TRUE, sep = "\n" )
    }else if( only_report == FALSE ){
      cp_files <- list.files( paste( "filtered_genes_step2/", i, "/", sep = "" ), 
                              pattern = "tree|aln", recursive = TRUE )
      file.copy( paste( "filtered_genes_step2/", i, "/", cp_files[1], sep = "" ),
                 paste( name_dir, "/", counter, "/", sep = ""), recursive = TRUE )
      file.copy( paste( "filtered_genes_step2/", i, "/", cp_files[2], sep = "" ),
                 paste( name_dir, "/", counter, "/", sep = ""), recursive = TRUE )
      write( paste( "Gene ", counter, " -- ", i, " :", name_gene, sep = "" ),
             paste( name_dir, "/log_R_subsample", out, ".txt", sep = "" ),
             append = TRUE, sep = "\n" )
    }
    
  }
  
  
}

## 11.1 Find indexes for the genes that will be later sampled
# Randomly create subsets of these genes. Note that sampling is done from the 
# already filtered and ordered genes. Therefore, use `sort` here as the order matters.
set.seed( 12345 )
r30sp  <- sort( sample( x = 1:length( list_correspondance_72sp ), size = 30 ) )
r50sp  <- sort( c( r30sp, sample( x = setdiff( (1:length( list_correspondance_72sp ) ), r30sp ), size = 20 ) ) )
r100sp <- sort( c( r50sp, sample( x = setdiff( (1:length( list_correspondance_72sp ) ), r50sp ), size = 50 ) ) )
r500sp <- sort( c( r100sp, sample( x = setdiff( (1:length( list_correspondance_72sp ) ), r100sp ), size = 400 ) ) )

# Now, add 355 additional genes randomly samples from the big pool -- not shared among all 72sp
# so we can have a subsample of up to 1,000 genes. The same: use `sort` as order matters.
files_all <- list.files( path = "filtered_genes_step2/", pattern = ".tree", recursive = TRUE )
names_all <- gsub( x = stringr::str_sort( files_all, numeric = TRUE ),
                   pattern = "[0-9].*/|.tree", replacement = "" )
names( names_all ) <- 1:length( names_all )
names_645 <- gsub( x = list.files( path = "filtered_genes_step2_all72sp/", pattern = ".tree", recursive = TRUE ),
                   pattern = "[0-9].*/|.tree", replacement = "" )
names( names_645 ) <- 1:length( names_645 )
ind_rmv        <- which( names_all %in% names_645 )
names_all_subs <- names_all[-ind_rmv]
ind_num        <- as.numeric( names( names_all_subs ) )
r355sp_1000    <- sort( sample( x = ind_num, size = 355 ) )
##>> CHECK
which( names_all[ r355sp_1000 ] %in% names_645 )
which( names_all[ r355sp_1000 ] %in% names_645[r500sp] )
##>> END CHECK

# Add extra data:
r1sp   <- sort( sample( x = r30sp, size = 1 ) )
r10sp  <- sort( c( r1sp, sample( x = setdiff( r30sp, r1sp ), size = 9 ) ) )
r9355sp_10000  <- sort( c( r355sp_1000, sample( x = setdiff( ind_num, r355sp_1000 ), size = 9000 ) ) )
##>> CHECK
which( names_all[ r9355sp_10000 ] %in% names_645 )
which( names_all[ r9355sp_10000 ] %in% names_645[r500sp] )
##>> END CHECK

# Sample now from the whole data set!
r10sp_all  <- sort( sample( x = 1:length( list_correspondance ), size = 10 ) )
r30sp_all  <- sort( sample( x = 1:length( list_correspondance ), size = 30 ) )
r100sp_all <- sort( sample( x = 1:length( list_correspondance ), size = 100 ) )
r500sp_all <- sort( sample( x = 1:length( list_correspondance ), size = 500 ) )
r1000sp_all <- sort( sample( x = 1:length( list_correspondance ), size = 1000 ) )

## 11.2 Generate data subsets using indexes based on genes that
#       are present in all 72 mammal taxa
if ( ! dir.exists( "sample1sp" ) ){
  dir.create( "sample1sp" )
}
subsample( range_num = r1sp, name_dir = "sample1sp", out = "1sp",
           extra = FALSE, only_report = TRUE )

if ( ! dir.exists( "sample10sp" ) ){
  dir.create( "sample10sp" )
}
subsample( range_num = r10sp, name_dir = "sample10sp", out = "10sp",
           extra = FALSE, only_report = TRUE )

if ( ! dir.exists( "sample30sp" ) ){
  dir.create( "sample30sp" )
}
subsample( range_num = r30sp, name_dir = "sample30sp", out = "30sp",
           extra = FALSE, only_report = TRUE )

if ( ! dir.exists( "sample50sp" ) ){
  dir.create( "sample50sp" )
}
subsample( range_num = r50sp, name_dir = "sample50sp", out = "50sp",
           extra = FALSE, only_report = TRUE )

if ( ! dir.exists( "sample100sp" ) ){
  dir.create( "sample100sp" )
}
subsample( range_num = r100sp, name_dir = "sample100sp", out = "100sp",
           extra = FALSE, only_report = TRUE )

if ( ! dir.exists( "sample500sp" ) ){
  dir.create( "sample500sp" )
}
subsample( range_num = r500sp, name_dir = "sample500sp", out = "500sp",
           extra = FALSE, only_report = TRUE )

if ( ! dir.exists( "sample1000sp" ) ){
  dir.create( "sample1000sp" )
}
subsample( range_num = r355sp_1000, name_dir = "sample1000sp", out = "1000sp",
           extra = TRUE, only_report = FALSE )

if ( ! dir.exists( "sample10000sp" ) ){
  dir.create( "sample10000sp" )
}
subsample( range_num = r9355sp_10000, name_dir = "sample10000sp", out = "10000sp",
           extra = TRUE, only_report = FALSE )

##>> CHECK
genes_samples <- vector( mode = "list", length = 6 )
dirs_avail <- c( "sample30sp", "sample50sp", "sample100sp", "sample500sp", "sample1000sp",
                 "sample10000sp" )
names( genes_samples ) <- dirs_avail 
for( i in 1:length( dirs_avail ) ){
  genes_samples[i] <- read.table( paste( dirs_avail[i], "/log_R_sub", dirs_avail[i], "_genes.txt", sep = "" ) )
}

intersect( genes_samples$sample500sp, genes_samples$sample1000sp )
intersect( genes_samples$sample500sp, genes_samples$sample10000sp )
##>> END CHECK

## 11.2 Generate data subsets using indexes based on ALL genes
if ( ! dir.exists( "sample10sp_all" ) ){
  dir.create( "sample10sp_all" )
}
subsample_fromall( range_num = r10sp_all, name_dir = "sample10sp_all",
                   out = "10sp_all", only_report = FALSE )

if ( ! dir.exists( "sample30sp_all" ) ){
  dir.create( "sample30sp_all" )
}
subsample_fromall( range_num = r30sp_all, name_dir = "sample30sp_all",
                   out = "30sp_all", only_report = FALSE )

if ( ! dir.exists( "sample100sp_all" ) ){
  dir.create( "sample100sp_all" )
}
subsample_fromall( range_num = r100sp_all, name_dir = "sample100sp_all",
                   out = "100sp_all", only_report = FALSE )

if ( ! dir.exists( "sample500sp_all" ) ){
  dir.create( "sample500sp_all" )
}
subsample_fromall( range_num = r500sp_all, name_dir = "sample500sp_all",
                   out = "500sp_all", only_report = FALSE )

if ( ! dir.exists( "sample1000sp_all" ) ){
  dir.create( "sample1000sp_all" )
}
subsample_fromall( range_num = r1000sp_all, name_dir = "sample1000sp_all",
                   out = "1000sp_all", only_report = FALSE )
 
#-----------------------------------------------------------------------------------------#
# 12. PARTITION THE ALIGNMENT (use bash script) --> ONLY ONE PARTITION WITH 1ST+2ND CPs!! #
#-----------------------------------------------------------------------------------------#