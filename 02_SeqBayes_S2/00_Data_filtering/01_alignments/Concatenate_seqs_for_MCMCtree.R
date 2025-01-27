#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------#
# SET WORKING DIR #
#-----------------#
# Load package to help find the path to this source file 
library(rstudioapi) 
# Get the path of current open file
path_to_file <- getActiveDocumentContext()$path 
# Get working directory path
wd      <- paste( dirname( path_to_file ), "/", sep = "" )
wd.name <- dirname( path_to_file )
# Set wd
setwd( wd )

#------------#
# LOAD FILES #
#------------#
# 0. Set paths and main subtree name 
subt    <- "xenarthra"
aln_dir <- paste( "../00_data_curation/", subt, "/filter_aln", sep = "" )

# 1. Find directories and create and object that contains 
#    the names of the directories of alignments in this 
#    directory and also those with only the 3nt in the 
#    other directory

## All nucs
path.to.files.raw <- list.files( path = aln_dir, all.files = TRUE, recursive = TRUE )
ind.path.to.files <- stringr::str_detect( string = path.to.files.raw, pattern = "t.phylip|checked_aln/alignment_nt3cp.phylip" )
path.to.files     <- path.to.files.raw[ ind.path.to.files ]
name.dirs <- subt

# Add main folder name to path 
path.to.files.all <- paste( aln_dir, path.to.files, sep = "/" )

## Only 4parts and only 3nt part
path.to.files     <- grep( pattern = "t.phylip", x = path.to.files.all, value = TRUE )
# Uncomment and run next line only for Euarchonta and Afrotheria if
# they have `alignment.phylip` files in `checked_aln`
# path.to.files     <- grep( pattern = "filter_aln/alignment.phylip", x = path.to.files.all, value = TRUE )
path.to.files.3nt <- grep( pattern = "nt3cp", x = path.to.files.all, value = TRUE )
## Partition data 
ind.path.to.partitions <- stringr::str_detect( string = path.to.files.raw, pattern = "partitions.txt" )
path.to.partitions     <- sort( path.to.files.raw[ ind.path.to.partitions] )

# Add main folder name to path 
path.to.partitions <- paste( aln_dir, path.to.partitions, sep = "/" )

#------------------------#
# Function to load files #
#------------------------#
# Arguments:
# path_aln = character, vector with the path to the
#            "alignment.phylip" files. 
# path_aln_3nt = character, vector with the path to the 
#                "alignment_3nt.phylip" files. Default: FALSE, as 
#                it assumes the alignment is performed without the 3nt.
# path_partitions = path to the "partitions.txt" files provided in
#                   the "checked" directories, which show when each
#                   partition finishes.
# dir_names = character, vector with the names of the directories
#             in the order they are read in the path
load.files <- function( path_aln, path_aln_3nt = FALSE, path_partitions, dir_names ){
  
  # 0. Break if you do not have the same length of alignments
  #    with 3nt and the rest
  if( path_aln_3nt != FALSE ){
    if( length( path_aln ) != length( path_aln_3nt) ){
      stop( "The name of alignments with 3nt and the rest should be the same" )
    }
  }else{
    cat( "You have not provided a partition with 3nt. Therefore, this will not be added!\n")
  }
  
  # 1. Create global vars
  if( path_aln_3nt != FALSE ){
    aln.list <- aln_3nt.list <- sp.list <- sp_3nt.list <- 
      nuc.list <- nuc_3nt.list <- partitions.list <-
      vector( mode = "list", length = length( path_aln ) )
  }else{
    aln.list <- sp.list <- nuc.list <- partitions.list <-
      vector( mode = "list", length = length( path_aln ) )
  }
  
  
  # 2. Get dummy name lists, later replace with correct names
  names( aln.list ) <- names( sp.list ) <- names( nuc.list ) <-
    names( partitions.list ) <- 1:length( path_aln )
  
  
  # 3. Read files and save them in a list
  for( i in 1:length( path_aln ) ){
    
    cat( "Loading file ", path_aln[i], "...\n" )
    tmp.aln        <- readLines( path_aln[i] )
    tmp.aln.header <- tmp.aln[1]
    sp_nuc <- length( strsplit( tmp.aln.header, " " )[[1]] )
    sp_ind <- sp_nuc - 1 
    tmp.sp  <- as.numeric( strsplit( tmp.aln.header, " " )[[1]][sp_ind] )
    tmp.nuc <- as.numeric( strsplit( tmp.aln.header, " " )[[1]][sp_nuc] )
    ## SAC-200330> OLD CODE. 
    # if( strsplit( tmp.aln.header, " " )[[1]][1] != "" ){
    #   tmp.sp  <- as.numeric( strsplit( tmp.aln.header, " " )[[1]][1] )
    #   tmp.nuc <- as.numeric( strsplit( tmp.aln.header, " " )[[1]][2] )
    # }else{
    #   tmp.sp  <- as.numeric( strsplit( tmp.aln.header, " " )[[1]][2] )
    #   tmp.nuc <- as.numeric( strsplit( tmp.aln.header, " " )[[1]][3] )
    # }
    
    tmp.aln <- tmp.aln[2:length(tmp.aln)]
    aln.list[[ i ]] <- tmp.aln
    sp.list[[ i ]]  <- tmp.sp
    nuc.list[[ i ]] <- tmp.nuc
    
    if( path_aln_3nt != FALSE ){
      cat( "Loading file ", path_aln_3nt[i], "...\n" )
      tmp.aln.3nt        <- readLines( path_aln_3nt[i] )
      tmp.aln.3nt.header <- tmp.aln.3nt[1]
      sp_nuc.3nt <- length( strsplit( tmp.aln.3nt.header, " " )[[1]] )
      sp_ind.3nt <- sp_nuc - 1 
      # tmp.sp.3nt  <- as.numeric( strsplit( tmp.aln.3nt.header, " " )[[1]][1] )
      tmp.sp.3nt  <- as.numeric( strsplit( tmp.aln.3nt.header, " " )[[1]][sp_ind.3nt] )
      all.equal( tmp.sp.3nt, tmp.sp )
      # tmp.nuc.3nt <- as.numeric( strsplit( tmp.aln.3nt.header, " " )[[1]][2] )
      tmp.nuc.3nt <- as.numeric( strsplit( tmp.aln.3nt.header, " " )[[1]][sp_nuc] )
      tmp.aln.3nt <- tmp.aln.3nt[2:length(tmp.aln.3nt)]
      aln_3nt.list[[ i ]] <- tmp.aln.3nt
      sp_3nt.list[[ i ]]  <- tmp.sp.3nt
      nuc_3nt.list[[ i ]] <- tmp.nuc.3nt
    }
    
    # Read partitions 
    tmp.parts.raw <- readLines( path_partitions[i] )
    tmp.parts       <- gsub( pattern = "..*= |\n", replacement = "", x = tmp.parts.raw )
    tmp.names.parts <- gsub( pattern = "DNA, | =..*", replacement = "", x = tmp.parts.raw )
    tmp.matrix <- matrix( 0, ncol = 2, nrow = 4 )
    rownames( tmp.matrix ) <- tmp.names.parts
    for( pos in 1:length( tmp.parts) ){
      tmp.rowcol <- as.numeric( strsplit( x = tmp.parts[pos], split = "-" )[[1]] )
      tmp.matrix[pos,1]      <- tmp.rowcol[1]
      tmp.matrix[pos,2]      <- tmp.rowcol[2]
      partitions.list[[ i ]] <- tmp.matrix
    }
    
    # Adding names to the entry lists 
    cat( "Directory", dir_names[i], "finished!!\n\n")
    if( path_aln_3nt != FALSE ){
      names( aln.list )[i] <- names( sp.list )[i] <- names( nuc.list )[i] <-
        names( aln_3nt.list )[i] <- names( sp_3nt.list )[i] <- names( nuc_3nt.list )[i] <-
        names( partitions.list )[i] <- dir_names[i]
    }else{
      names( aln.list )[i] <- names( sp.list )[i] <- names( nuc.list )[i] <-
        names( partitions.list )[i] <- dir_names[i]
    }
    
    
  }
  
  # 4. Return objects
  if( path_aln_3nt != FALSE ){
    return( list( aln = aln.list, aln_3nt = aln_3nt.list,
                  sp  = sp.list,  sp_3nt  = sp_3nt.list,
                  nuc = nuc.list, nuc_3nt = nuc_3nt.list,
                  parts = partitions.list ) )
  }else{
    return( list( aln = aln.list,
                  sp  = sp.list,
                  nuc = nuc.list,
                  parts = partitions.list ) )
  }
}

# 2. Read files and save them in object
aln.object <- load.files( path_aln = path.to.files, path_aln_3nt = path.to.files.3nt,
                          dir_names = name.dirs, path_partitions = path.to.partitions )

if ( ! dir.exists( "Rout/" ) ){
  dir.create( "Rout/"  )
  dir.create( "Rout/Rdata" ) 
  dir.create( "Rout/log_concatenation" ) 
}

save( aln.object, file = paste( "Rout/Rdata/aln.object.", subt, ".RData", sep = "" ) )
# Uncomment if you just want to load this object!
# load( paste( "Rout/Rdata/aln.object.", subt, ".RData", sep = "" ) )

# 3. Prepare list vectors to save sequences

#-------------------------------------------#
# Function to save the alignments that have #
# been previously read in correct format    #
#-------------------------------------------#
# Argument:
#   x   = list, object generated by function "load.files"
#   nt3 = boolean. Default is FALSE as it assumes the 3nt is not available
store.aln <- function( x, nt3 = FALSE ){
  
  # 0. Get vars out of list passed through argument "x"
  tmp.aln <- x$aln
  tmp.sp  <- x$sp
  tmp.nuc <- x$nuc
  if( nt3 != FALSE ){
    tmp.aln.3nt <- x$aln_3nt
    tmp.sp.3nt  <- x$sp_3nt
    tmp.nuc.3nt <- x$nuc_3nt
    all.equal( names( tmp.sp ), names( tmp.sp.3nt) )
  }
  tmp.parts <- x$parts
  
  # SAC-20/03/30> Create list to save returned matrices to check alignment 3nt has been
  # properly appended
  mat.return <- vector( mode = "list", length( tmp.aln ) )
  
  # 1. Check if directory where to copy files exist. Otherwise, create it
  if ( ! dir.exists( "00_mammal_alns/" ) ){
    dir.create( "00_mammal_alns/"  )
  }
  
  
  # 2. Start for loop that goes over each of the files read
  for( k in 1:length( tmp.aln ) ){
    
    # 0. Get objects in position "k" and store them in tmp vars
    #    and create dir for output aln
    cat( "Parsing alignments for ", names( tmp.aln )[k], "...\n" )
    write( cat( "Parsing alignments for ", names( tmp.aln )[k], "...\n" ),
           paste( "Rout/log_concatenation/log_file_concatenate_seqs_", subt, ".txt", sep = "" ),
           append = TRUE, sep = "\n" )
    tmp.aln.k <- tmp.aln[[ k ]]
    tmp.sp.k  <- tmp.sp[[ k ]]
    tmp.nuc.k <- tmp.nuc[[ k ]]
    if( nt3 != FALSE ){
      tmp.aln.3nt.k <- tmp.aln.3nt[[ k ]]
      tmp.sp.3nt.k  <- tmp.sp.3nt[[ k ]]
      tmp.nuc.3nt.k <- tmp.nuc.3nt[[ k ]]
    }
    
    tmp.parts.k  <- tmp.parts[[ k ]]
    
    # Check if directory where to copy files exist. Otherwise, create it
    if ( ! dir.exists( paste( "00_mammal_alns/", names( tmp.aln )[k], sep = "" ) ) ){
      dir.create( paste( "00_mammal_alns/", names( tmp.aln )[k], sep = "" ) )
    }
    tmp.dirname.k <- paste( "00_mammal_alns/", names( tmp.aln )[k], "/", sep = "" )
    
    # 1. Create temporary elements
    tmp.list           <- vector( mode = "list", length = tmp.sp.k )
    tmp.aln.k.names.sp <- rep( x = "", times = tmp.sp.k )
    if( nt3 != FALSE ){
      tmp.3nt.names.sp   <- rep( x = "", times = tmp.sp.3nt.k )
      cat( "Same species in both 3nt and main alignment?\n" )
      write( paste( "Same species in both 3nt and main alignment?", sep = "" ),
             paste( "Rout/log_concatenation/log_file_concatenate_seqs_", subt, ".txt", sep = "" ),
             append = TRUE, sep = "\n" )
      all.equal( tmp.aln.k.names.sp, tmp.3nt.names.sp )
      tmp.txt <- all.equal( tmp.aln.k.names.sp, tmp.3nt.names.sp )
      write( tmp.txt, 
             paste( "Rout/log_concatenation/log_file_concatenate_seqs_", subt, ".txt", sep = "" ),
             append = TRUE, sep = "\n" )
    }
    
    if( nt3 != FALSE ){
      tmp.nt12cp.k <-tmp.mt12cp.k <- tmp.mt3cp.k <- tmp.mtrna.k <- 
        tmp.nt3cp.k <- vector( mode = "list", length = tmp.sp.k )
    }else{
      tmp.nt12cp.k <-tmp.mt12cp.k <- tmp.mt3cp.k <- tmp.mtrna.k <- 
        vector( mode = "list", length = tmp.sp.k )
    }
    
    # 2. Start with main alignment: read line by line and process
    breaks <- 0 
    cat( "Reading sequences in main alignment...\n" )
    write( paste( "Reading sequences in main alignment...\n", sep = "" ),
           paste( "Rout/log_concatenation/log_file_concatenate_seqs_", subt, ".txt", sep = "" ), 
           append = TRUE, sep = "\n"  )
    for( i in 1:length( tmp.aln.k ) ){
      
      # I have noticed that some of the alignments.phylip have 
      # spaces and others do not. 
      if( i == 1 & length( tmp.aln.k ) == tmp.sp.k ){
        
        cat( "We have been dealing with the other format now where no spaces!\n" )
        new.format <- TRUE
        
      }else{
        new.format <- FALSE
      }
      # While there are seqnames (first 63 lines for Afrotheria, for instance)
      if( i <= tmp.sp.k ){
        
        tmp.line              <- strsplit( x = tmp.aln.k[i], split = "   " )
        tmp.aln.k.names.sp[i] <- tmp.line[[1]][1]
        last.entry            <- length( tmp.line[[1]] )
        # SAC-20/03/23> Needed to update this as otherwise it was getting only spaces
        #               with the new format in some phylip files!
        #tmp.list[[ i ]]  <- tmp.line[[1]][2]
        #tmp.list[[ i ]]  <- tmp.line[[1]][last.entry]
        # Make sure it does not have spaces nor tabs !
        tmp.checked.aln <- gsub( pattern = " ", replacement = "", x = tmp.line[[1]][last.entry] )
        tmp.checked.aln <- gsub( pattern = "\t", replacement = "", x = tmp.checked.aln )
        tmp.list[[ i ]] <- tmp.checked.aln
        
        cat( "Species visited: ", tmp.line[[1]][1], "\n" )
        write( paste( "Species visited: ", tmp.line[[1]][1], sep = "" ),
               paste( "Rout/log_concatenation/log_file_concatenate_seqs_", subt, ".txt", sep = "" ), 
               append = TRUE, sep = "\n" )
        
        if ( i == tmp.sp.k ){
          names( tmp.list ) <- tmp.aln.k.names.sp
        }
        
      }else if( tmp.aln.k[i] == "" ){ # If empty line, create j and add sequences to list
        # This message is just to be printed for every blank line,
        # excluding the last line of the file
        if ( i < length( tmp.aln.k ) ){
          cat( "We are in a break line now! More sequences to be appended now to
               the corresponding species...\n" )
          write( paste( "We are in a break line now! More sequences to be appended now to
               the corresponding species...\n", sep = "" ),
                 paste( "Rout/log_concatenation/log_file_concatenate_seqs_", subt, ".txt", sep = "" ),
                 append = TRUE, sep = "\n"  )
        }
        j <- 0
        breaks <- breaks + 1
        
      }else{ # Now append the sequences to the list in corresponding order
        
        if( new.format != TRUE ){
          j <- j+1
          cat( "Appending sequence to taxon ", j, ":",
               tmp.aln.k.names.sp[j], "\n" )
          write( paste( "Appending sequence to taxon ", j, ":",
                        tmp.aln.k.names.sp[j], sep = "" ),
                 paste( "Rout/log_concatenation/log_file_concatenate_seqs_", subt, ".txt", sep = "" ),
                 append = TRUE, sep = "\n" )
          tmp.list[[ j ]] <- paste( tmp.list[[ j ]],
                                    tmp.aln.k[i],sep = "" )
        }
        
      }
      
      # Print summary stats
      if( i == length( tmp.aln.k )  ){
        cat( "\nTotal num of breaks = ", breaks, "\n" )
        write( paste( "\nTotal num of breaks = ", breaks, "\n", sep = "" ),
               paste( "Rout/log_concatenation/log_file_concatenate_seqs_", subt, ".txt", sep = "" ),
               append = TRUE, sep = "\n" )
      }
      
    }
    
    
    # 2.1 // # SAC-20/03/20> This step is required as of 20/03/23 as we do not have yet 3nt data !
    #        # This will be skipped as soon as we have the 3nt as this bit is included
    #        # in next for loop in step 4!
    if( nt3 == FALSE ){
      order.sp <- rep( 0, tmp.sp.k )
      cat( "\nThis block should only run if 12CP only!!\n" )
      cat( "\nAdding tmp alignments for each partition...\n" )
      write( paste( "\nThis block should only run if 12CP only!!\n", sep = "" ),
             paste( "Rout/log_concatenation/log_file_concatenate_seqs_", subt, ".txt", sep = "" ),
             append = TRUE, sep = "\n" )
      write( paste( "\nAdding tmp alignments for each partition...\n", sep = "" ),
             paste( "Rout/log_concatenation/log_file_concatenate_seqs_", subt, ".txt", sep = "" ),
             append = TRUE, sep = "\n" )
      for( i in 1:tmp.sp.k ){
  
        # Now take advantage of the situation and extract as well individual data!
        tmp.line <- strsplit( x = tmp.aln.k[i], split = "   " )
        tmp.spname <- tmp.line[[1]][1]
        index <- which( names( tmp.list ) == tmp.spname )
        order.sp[i] <- index
        # First, get partitioned data for species in pos "index"
        single.nucs <- strsplit( x = tmp.list[[ index ]], split = "" )[[1]]
  
        tmp.nt12cp.k[[ index ]] <- paste0( single.nucs[ tmp.parts.k[1,1]:tmp.parts.k[1,2] ], collapse = "" )
        nucs.nt12cp.k           <- length( single.nucs[ tmp.parts.k[1,1]:tmp.parts.k[1,2] ] )
        tmp.mt12cp.k[[ index ]] <- paste0( single.nucs[ tmp.parts.k[2,1]:tmp.parts.k[2,2] ], collapse = "" )
        nucs.mt12cp.k           <- length( single.nucs[ tmp.parts.k[2,1]:tmp.parts.k[2,2] ] )
        tmp.mt3cp.k[[ index ]]  <- paste0( single.nucs[ tmp.parts.k[3,1]:tmp.parts.k[3,2] ], collapse = "" )
        nucs.mt3cp.k            <- length( single.nucs[ tmp.parts.k[3,1]:tmp.parts.k[3,2] ] )
        tmp.mtrna.k[[ index ]]  <- paste0( single.nucs[ tmp.parts.k[4,1]:tmp.parts.k[4,2] ], collapse = "" )
        nucs.mtrna.k            <- length( single.nucs[ tmp.parts.k[4,1]:tmp.parts.k[4,2] ] )
  
      }
    }
    
    # 3. Add names to the lists for the partitions
    if( nt3 != FALSE ){
      names( tmp.nt12cp.k ) <- names( tmp.mt12cp.k ) <- names( tmp.mt3cp.k ) <- 
        names( tmp.mtrna.k ) <- names( tmp.nt3cp.k ) <- names( tmp.list )
    }else{
      names( tmp.nt12cp.k ) <- names( tmp.mt12cp.k ) <- names( tmp.mt3cp.k ) <- 
        names( tmp.mtrna.k ) <- names( tmp.list )
    }
    
    
    # 4. Add 3nt to alignment
    if( nt3 != FALSE ){
      breaks <- 0
      order.sp <- bool.check <- sp.check <- rep(0, tmp.sp.3nt.k )
      cat( "Reading sequences with 3nt to append them to previous alignment...\n" )
      write( paste( "Reading sequences with 3nt to append them to previous alignment...\n", sep = "" ),
             paste( "Rout/log_concatenation/log_file_concatenate_seqs_", subt, ".txt", sep = "" ),
             append = TRUE, sep = "\n" )
      for( i in 1:length( tmp.aln.3nt.k ) ){
        
        # While there are seqnames
        # We also take advantage of this conditional and partition the sequences according 
        # to the 4 different types of data
        if( i <= tmp.sp.3nt.k ){
          
          tmp.line   <- strsplit( x = tmp.aln.3nt.k[i], split = "   " )
          tmp.spname <- tmp.line[[1]][1]
          index      <- which( names( tmp.list ) == tmp.spname )
          bool.check[i] <- all.equal( names( tmp.list )[index], tmp.spname )
          sp.check[i]   <- tmp.spname
          order.sp[i]   <- index
          # First, get partitioned data for species in pos "index" (in aln without 3nt)
          single.nucs <- strsplit( x = tmp.list[[ index ]], split = "" )[[1]]
          
          tmp.nt12cp.k[[ index ]] <- paste0( single.nucs[ tmp.parts.k[1,1]:tmp.parts.k[1,2] ], collapse = "" )
          nucs.nt12cp.k           <- length( single.nucs[ tmp.parts.k[1,1]:tmp.parts.k[1,2] ] )
          tmp.mt12cp.k[[ index ]] <- paste0( single.nucs[ tmp.parts.k[2,1]:tmp.parts.k[2,2] ], collapse = "" )
          nucs.mt12cp.k           <- length( single.nucs[ tmp.parts.k[2,1]:tmp.parts.k[2,2] ] )
          tmp.mt3cp.k[[ index ]]  <- paste0( single.nucs[ tmp.parts.k[3,1]:tmp.parts.k[3,2] ], collapse = "" )
          nucs.mt3cp.k            <- length( single.nucs[ tmp.parts.k[3,1]:tmp.parts.k[3,2] ] )
          tmp.mtrna.k[[ index ]]  <- paste0( single.nucs[ tmp.parts.k[4,1]:tmp.parts.k[4,2] ], collapse = "" )
          nucs.mtrna.k            <- length( single.nucs[ tmp.parts.k[4,1]:tmp.parts.k[4,2] ] )
          
          
          # Now, paste the 3nt sequence
          # SAC-20/03/30> Needed to update this as otherwise it was getting only spaces
          #               with the new format in some phylip files!
          #               This is same fixed that was done in previous lines with rest
          #               of alignment
          last.entry            <- length( tmp.line[[1]] )
          # Make sure it does not have spaces nor tabs !
          tmp.checked.aln.3nt <- gsub( pattern = " ", replacement = "", x = tmp.line[[1]][last.entry] )
          tmp.checked.aln.3nt <- gsub( pattern = "\t", replacement = "", x = tmp.checked.aln.3nt )
          
          tmp.list[[ index ]] <- paste( tmp.list[[ index ]],
                                        #tmp.line[[1]][2], sep = "" )   # SAC-20/03/30> Changed. Notes in comment above.
                                        tmp.checked.aln.3nt, sep = "" ) # SAC-20/03/30> Update for commented line above.
          
          # tmp.nt3cp.k[[ index ]] <- tmp.line[[1]][2] # SAC-20/03/30> Changed. Notes in comment above.
          tmp.nt3cp.k[[ index ]] <- tmp.checked.aln.3nt   # SAC-20/03/30> Update for commented line above.
          
          cat( "Species visited: ", tmp.line[[1]][1], "\n" )
          write( paste( "Species visited: ", tmp.line[[1]][1], sep = "" ),
                 paste( "Rout/log_concatenation/log_file_concatenate_seqs_", subt, ".txt", sep = "" ),
                 append = TRUE, sep = "\n" )
          
        }
        # If empty line, create j and add sequences to list
        else if( tmp.aln.3nt.k[i] == "" ){
          
          # This message is just to be printed for every blank line,
          # excluding the last line of the file
          if ( i < length( tmp.aln.3nt.k ) ){
            cat( "We are in a break line now! More sequence to be appended now...\n" )
            write( paste( "We are in a break line now! More sequence to be appended now...\n", sep = "" ),
                   paste( "Rout/log_concatenation/log_file_concatenate_seqs_", subt, ".txt", sep = "" ),
                   append = TRUE, sep = "\n"  )
          }
          j <- 0
          breaks <- breaks + 1
          
        }
        # In case old format, then append the sequences to the list in corresponding order
        else{
          j <- j+1
          cat( "Appending sequence to taxon ", j, ":",
               tmp.aln.k.names.sp[ order.sp[j] ], "\n" )
          write( paste( "Appending sequence to taxon ", j, ":",
                        tmp.aln.k.names.sp[ order.sp[j] ], sep = "" ),
                 paste( "Rout/log_concatenation/log_file_concatenate_seqs_", subt, ".txt", sep = "" ),
                 append = TRUE, sep = "\n" )
          tmp.list[[ order.sp[j] ]] <- paste( tmp.list[[ order.sp[j] ]],
                                              tmp.aln.3nt.k[i], sep = "" )
          tmp.nt3cp.k[[ order.sp[j] ]] <- paste( tmp.nt3cp.k[[ order.sp[j] ]],
                                                 tmp.aln.3nt.k[i], sep = "" )
        }
        
        # Print summary stats
        if( i == length( tmp.aln.3nt.k )  ){
          cat( "Total num of breaks = ", breaks, "\n" )
          write( paste( "Total num of breaks = ", breaks, "\n", sep = "" ),
                 paste( "Rout/log_concatenation/log_file_concatenate_seqs_", subt, ".txt", sep = "" ),
                 append = TRUE, sep = "\n" )
          nucs.nt3cp.k <- nchar( tmp.nt3cp.k[[ 1 ]] )
        }
        
      } 
    }
    
    
    # 4. Write down header of alignment 
    if( nt3 != FALSE ){
      total.nuc <- tmp.nuc.k + tmp.nuc.3nt.k
    }else{
      total.nuc <- tmp.nuc.k
    }
    
    write( x = paste( length( tmp.aln.k.names.sp ), total.nuc, sep = "  " ),
           file = paste( tmp.dirname.k, names( tmp.aln )[k], ".aln", sep = ""  ) )
    
    write( x = paste( length( tmp.aln.k.names.sp ), nucs.nt12cp.k, sep = "  " ),
           file = paste( tmp.dirname.k, names( tmp.aln )[k], "_nt12cp.aln", sep = ""  ) )
    write( x = paste( length( tmp.aln.k.names.sp ), nucs.mt12cp.k, sep = "  " ),
           file = paste( tmp.dirname.k, names( tmp.aln )[k], "_mt12cp.aln", sep = ""  ) )
    write( x = paste( length( tmp.aln.k.names.sp ), nucs.mt3cp.k, sep = "  " ),
           file = paste( tmp.dirname.k, names( tmp.aln )[k], "_mt3cp.aln", sep = ""  ) )
    write( x = paste( length( tmp.aln.k.names.sp ), nucs.mtrna.k, sep = "  " ),
           file = paste( tmp.dirname.k, names( tmp.aln )[k], "_mtrna.aln", sep = ""  ) )
    if( nt3 != FALSE ){
      write( x = paste( length( tmp.aln.k.names.sp ), nucs.nt3cp.k, sep = "  " ),
             file = paste( tmp.dirname.k, names( tmp.aln )[k], "_nt3cp.aln", sep = ""  ) )
    }
    
    
    # 5. Add sequences to file 
    for( i in 1:tmp.sp.k ){
      
      tmp.nameINlist <- names( tmp.list )[ i ]
      tmp.seqINlist  <- tmp.list[[ i ]]
      
      print( all.equal( nchar( tmp.list[[ i ]]), total.nuc ) )
      write( x = paste( tmp.nameINlist, tmp.seqINlist, sep = "   " ),
             file = paste( tmp.dirname.k, names( tmp.aln )[k], ".aln", sep = ""  ),
             append = T )
      
      tmp.nameINlist.nt12cp <- names( tmp.nt12cp.k )[ i ]
      tmp.seqINlist.nt12cp  <- tmp.nt12cp.k[[ i ]]
      tmp.nameINlist.mt12cp <- names( tmp.mt12cp.k )[ i ]
      tmp.seqINlist.mt12cp  <- tmp.mt12cp.k[[ i ]]
      tmp.nameINlist.mt3cp  <- names( tmp.mt3cp.k )[ i ]
      tmp.seqINlist.mt3cp   <- tmp.mt3cp.k[[ i ]]
      tmp.nameINlist.mtrna  <- names( tmp.mtrna.k )[ i ]
      tmp.seqINlist.mtrna   <- tmp.mtrna.k[[ i ]]
      if( nt3 != FALSE ){
        tmp.nameINlist.nt3cp  <- names( tmp.nt3cp.k )[ i ]
        tmp.seqINlist.nt3cp   <- tmp.nt3cp.k[[ i ]]
      }
      
      
      write( x = paste( tmp.nameINlist.nt12cp, tmp.seqINlist.nt12cp, sep = "   " ),
             file = paste( tmp.dirname.k, names( tmp.aln )[k], "_nt12cp.aln", sep = ""  ),
             append = T )
      write( x = paste( tmp.nameINlist.mt12cp, tmp.seqINlist.mt12cp, sep = "   " ),
             file = paste( tmp.dirname.k, names( tmp.aln )[k], "_mt12cp.aln", sep = ""  ),
             append = T )
      write( x = paste( tmp.nameINlist.mt3cp, tmp.seqINlist.mt3cp, sep = "   " ),
             file = paste( tmp.dirname.k, names( tmp.aln )[k], "_mt3cp.aln", sep = ""  ),
             append = T )
      write( x = paste( tmp.nameINlist.mtrna, tmp.seqINlist.mtrna, sep = "   " ),
             file = paste( tmp.dirname.k, names( tmp.aln )[k], "_mtrna.aln", sep = ""  ),
             append = T )
      if( nt3 != FALSE ){
        write( x = paste( tmp.nameINlist.nt3cp, tmp.seqINlist.nt3cp, sep = "   " ),
               file = paste( tmp.dirname.k, names( tmp.aln )[k], "_nt3cp.aln", sep = ""  ),
               append = T )
      }
      
      
    }
    
    # Write new partitions file including nt_3cp
    new.parts.k          <- matrix(0, nrow = 5, ncol = 2 )
    new.parts.k[c(1:4),] <- tmp.parts.k
    new.parts.k[5,1]     <- tmp.parts.k[4,2] + 1
    new.parts.k[5,2]     <- total.nuc
    
    part.DNAinfo <- paste( "DNA, ", c( rownames( tmp.parts.k ), "nt_3cp" ), " =", sep = "" )
    part.NUCinfo <- paste( new.parts.k[,1], new.parts.k[,2], sep = "-" )
    part.out     <- paste( part.DNAinfo, part.NUCinfo, sep = " " )
    write( x = part.out,
           file = paste( tmp.dirname.k, "partitions.txt", sep = ""  ) )
    
    # Last check: have we appended species properly from 3nt_CP to the rest of alignment?
    if( nt3 != FALSE ){
      tmp.mat.return             <- matrix( 0, nrow = length( order.sp ), ncol = 3 )
      colnames( tmp.mat.return ) <- c( "order.sp", "bool.check", "sp.check" )
      tmp.mat.return[,1]         <- order.sp
      tmp.mat.return[,2]         <- bool.check
      tmp.mat.return[,3]         <- sp.check
      names( mat.return )[k]     <- names( tmp.aln )[k]
      mat.return[[ k ]]          <- tmp.mat.return
      write.table( x = tmp.mat.return,
                   file = paste( "Rout/Rdata/check_appends_", subt, ".txt", sep = ""  ),
                   quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t" )
    }
  }
  
  # Return matrix with checks
  if( nt3 != FALSE ){
    return( mat.return )
  }
  
}

checks <- store.aln( x = aln.object, nt3 = TRUE )  



