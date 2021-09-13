#--------------------------#
# EXAMINE TAXONOMIC LEVELS #
#--------------------------#
# Compare if taxa belonging to another have genes that do 
# not share amongsy any of the others.
# If so, highlight them for posterior deletion
#
# There will be a log file "log_taxaNOTIN_order.txt" to
# keep track of the findings during filtering
#
# Arguments: 
#     unique.lin   character, vector with the name of all the entries for the taxonomical levels being checked. 
#     sp           character, vector with all the taxa names for the specific taxonomic level being checked.
#     genes        character, vector with all the lines inside the file "genes.txt".
#     name.lin     character, vector with the name of the taxonomic level being checked. E.g. "Family".
filter.lineages <- function( unique.lin, sp, genes, name.lin ){
  
  checked.lev          <- vector( mode = "list", length =  length( unique.lin ) )
  names( checked.lev ) <- unique.lin 
  count.lin            <- 0 
  for( lin in unique.lin ){
    
    count.lin <- count.lin + 1
    write( x = paste( "--------------------------------------------", sep = "" ),
           file = paste( "log_taxaNOTIN_", name.lin, ".txt", sep = "" ), append = TRUE )
    write( x = paste( "    ANALYSIS FOR ", toupper( name.lin ), " ", toupper( lin ), sep = "" ),
           file = paste( "log_taxaNOTIN_", name.lin, ".txt", sep = "" ), append = TRUE )
    write( x = paste( "--------------------------------------------\n\n", sep = "" ),
           file = paste( "log_taxaNOTIN_", name.lin, ".txt", sep = "" ), append = TRUE )
    # 1. Find index where order "lin" is in vector "lins"
    ind.lin <- which( lins == lin )
    # 2. Count num. species
    name.sp <- sp[ ind.lin ]
    # 3. Print sp 
    write( x =  paste( "For ", name.lin, " ", lin, ", there are ", length( name.sp ), " entries:\n", sep = "" ),
           file = paste( "log_taxaNOTIN_", name.lin, ".txt", sep = "" ), append = TRUE )
    write( x =  paste( name.sp, sep = "\t" ), file = paste( "log_taxaNOTIN_", name.lin, ".txt", sep = "" ), append = TRUE )
    write( x =  paste( "\n", sep = "" ), file = paste( "log_taxaNOTIN_", name.lin, ".txt", sep = "" ), append = TRUE )
    
    # 4. Find which genes these taxa have and if they share any of them
    list.genes          <- vector( mode = "list", length = length( name.sp ) )
    names( list.genes ) <- name.sp
    for( i in 1:length( name.sp ) ){
      tmp.sp     <- name.sp[ i ]
      matched.l1 <- grep( pattern = paste( tmp.sp, "[^_..*]", sep = "" ), x = genes, value = TRUE )
      matched.l2 <- grep( pattern = paste( tmp.sp, "$", sep = "" ), x = genes, value = TRUE )
      matched.l  <- unique( union( matched.l1, matched.l2 ) )
      tmp.genes    <- unique( sapply( stringr::str_split( string = matched.l, pattern = " " ),
                                      function( x ) x[2] ) )
      list.genes[[ i ]] <- tmp.genes
    }
    
    # 5. Compare 
    for( i in 1:length( list.genes ) ){
      
      # Start with entry i
      #cat( "> We start comparing ", names( list.genes )[i], "against the rest of entries ...\n" )
      write( x =  paste( "> We start comparing ", names( list.genes )[i],
                         " against the rest of entries ...", sep = "" ),
             file = paste( "log_taxaNOTIN_", name.lin, ".txt", sep = "" ), append = TRUE )
      # Now compare one by one:
      count   <- 0 
      for( j in 1:length( list.genes ) ){
        
        # Avoid comparing "i" with itself if it matches position "j"
        if( i != j ){
          tmp <- intersect( list.genes[[ i ]], list.genes[[ j ]] )
          #print( tmp )
          if( length( tmp ) == 0 ){
            count <- count + 1
            #cat( names( list.genes )[i], "does not have any gene in common with ",
            #     names( list.genes )[j], "\n" )
            write( x =  paste( names( list.genes )[i], " does not have any gene in common with ",
                               names( list.genes )[j], sep = "" ),
                   file = paste( "log_taxaNOTIN_", name.lin, ".txt", sep = "" ), append = TRUE )
          }
          
        }
        # Print sum. stats when we have visited all entries in "list.genes"
        if( j == length( list.genes ) && count != 0 ){
          # cat( "Total = ", count, "/", length( name.sp )-1, "possible taxa comparisons\n")
          write( x = paste( "Total = ", count, "/", length( name.sp )-1,
                            " possible taxa comparisons\n", sep = "" ),
                 file = paste( "log_taxaNOTIN_", name.lin, ".txt", sep = "" ), append = TRUE )
        }
        
        # Print out which taxa needs to be 100% sure deleted!
        if( i !=j & count == length( name.sp )-1 & length( name.sp ) != 1 ){
          # Find position of species that needs to be deleted in "list.genes"
          sp.eval_i.ind <- which( names( list.genes ) == names( list.genes )[i] )
          cat( toupper( name.lin ), " ", lin, " has ", length( name.sp ), " taxa. Species ", 
               names( list.genes )[i], " has ", length( list.genes[[ i ]] ), 
               " genes that are not shared by any of the other taxa. It should be deleted!\n" )
          write( x = paste( toupper( name.lin ), " ", lin, " has ", length( name.sp ), " taxa. Species ", 
                            names( list.genes )[i], " has ", length( list.genes[[ i ]] ), 
                            " genes that are not shared by any of the other taxa. It should be deleted!\n",
                            sep = "" ),
                 file = paste( "log_taxaNOTIN_", name.lin, ".txt", sep = "" ), append = TRUE )
        }
        
      }
      #cat( "\n ")
      
    }
    
    #cat( ">> END << \n\n" )
    write( x = paste( ">> END << \n", sep = "" ),
           file = paste( "log_taxaNOTIN_", name.lin, ".txt", sep = "" ), append = TRUE )
    
    # Return the list.genes for that specific "lin"
    checked.lev[[ count.lin ]] <- list.genes
    
  }
  
  # Return the "list.genes" var
  return( genes = checked.lev )
  
}