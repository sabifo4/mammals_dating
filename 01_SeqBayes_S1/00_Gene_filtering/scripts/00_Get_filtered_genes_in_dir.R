#----------------------#
# 0. CLEAN ENVIRONMENT #
#----------------------#
rm( list = ls( ) )

#--------------------------#
# 1. SET WORKING DIRECTORY #
#--------------------------#
library( rstudioapi ) 
# Get the path to current open R script and find main dir "00_Gene_filtering"
path_to_file <- getActiveDocumentContext()$path
script_wd <- paste( dirname( path_to_file ), "/", sep = "" )
wd <- gsub( pattern = "/scripts", replacement = "", x = script_wd )
setwd( wd )

#---------------------------------------------#
# 2. READ FILE WITH SUMM STATISTICS FOR GENES #
#---------------------------------------------#
# Read file with summary genes 
sum_genes <- read.csv( file = "out_data/genes.csv" )

#--------------------------------------#
# 3. APPLY RESTRICTIONS TO FILTER DATA #
#--------------------------------------#
# Find which genes do not have mouse sequence
no_mouse_ind  <- which( sum_genes$Mouse_presence == 0 )
length( no_mouse_ind )   # 2 

# Find which genes do not have human sequence 
no_human_ind <- which( sum_genes$Human_presence == 0 )
length( no_human_ind )   # 0

# Find which genes have under 100 codons 
less_cod_ind  <- which( sum_genes$Num_codons_hum < 100 ) 
length( less_cod_ind )   # 323

# Find which genes have less than 10 sequences 
less_10seq_ind  <- which( sum_genes$Num_seqs < 10 ) 
length( less_10seq_ind ) # 18

# Find out which genes do not meet at least one of the previous 
# restrictions 
restrictions <- unique( c( no_mouse_ind, no_human_ind, less_cod_ind, less_10seq_ind ) )
length( restrictions )   # 335

# >> CHECK | Double check everything is alright
diff.ab <- length( which( less_10seq_ind %in% less_cod_ind ) )
diff.ab # 8
diff.ac <- length( which( less_10seq_ind %in% no_human_ind ) )
diff.ac # 0
diff.bc <- length( which( less_cod_ind %in% no_human_ind ) )  
diff.bc # 0

all.equal( ( length( no_mouse_ind ) + length( no_human_ind ) + length( less_cod_ind ) +
               length( less_10seq_ind ) ),
           ( length( restrictions ) + (diff.ab + diff.ac + diff.bc ) ) )
#[1] TRUE # 343 (there are duplicates) == 335 (unique) + 8 (shared)
# >> END CHECK

# Get filtered genes
filt_genes      <- sum_genes[-c( restrictions ),]
dim( filt_genes )[1] # 15569 genes
name_filt_genes <- sum_genes[-c( restrictions ),1]

# >> CHECK | Check we have processed them alright 
all.equal( which( sum_genes$Gene %in% sum_genes$Gene[restrictions] ),
           sort( restrictions ) )
# [1] TRUE
# >> END CHECK

# Check if directory where to copy files exist. Otherwise, create it
if ( ! dir.exists( "filtered_genes/" ) ){
  dir.create( "filtered_genes/" )
}

# Copy filtered genes in a new dir "filtered_genes"
counter = 0 
for( i in seq( 1:length( name_filt_genes ) ) ){
  
  counter = counter + 1
  # Copy files
  file.copy( paste( "genes/", as.character( name_filt_genes[i] ), "/", sep = "" ),
             paste( "filtered_genes/", sep = ""),
             recursive = TRUE )
  cat( "Gene ", counter, " is: ", as.character( name_filt_genes[i] ), "\n" )
  write( paste( "Gene ", counter, " is: ",
                as.character( name_filt_genes[i] ), "\n", sep = "" ),
         "out_logs/log_01_copy_filtered_genes.txt",
         append = TRUE, sep = "\n" )
}

#-------------------------------------------------------------#
# 4. APPLY EXTRA RESTRICTION: GENES PRESENT IN ALL 72 SPECIES #
#-------------------------------------------------------------#
# Find which genes are present in the 72 species
# and meet the requirements
ind_all_72  <- which( sum_genes$Num_seqs == 72 )
all72sp     <- sum_genes[c(ind_all_72),] # 648 genes !

# >> CHECK | Check all the genes are present in 72 species!
all.equal( length( which( all72sp$Num_seqs == 72 ) ), length( ind_all_72 ) )
#[1] TRUE 
# >> CHECK END

# Find which genes do not have human sequence
no_human_72sp_ind  <- which( all72sp$Human_presence == 0 ) 
length( no_human_72sp_ind ) # 0 

# Find which genes do not have mouse sequence
no_mouse_72sp_ind  <- which( all72sp$Mouse_presence == 0 ) 
length( no_mouse_72sp_ind ) # 0 

# Find which genes have under 100 codons 
less_cod_72sp_ind  <- which( all72sp$Num_codons_hum < 100 ) 
length( less_cod_72sp_ind ) # 0 

# Find which genes have less than 10 sequences 
less_10seq_72sp_ind  <- which( all72sp$Num_seqs < 10 )
length( less_10seq_72sp_ind ) #0 

# Get genes names
name_genes_72sp <- all72sp[,1]
length( name_genes_72sp ) # 648

# Check if directory where to copy files exist. Otherwise, create it
if ( ! dir.exists( "filtered_genes_all72sp/" ) ){
  dir.create( "filtered_genes_all72sp/" )
}

# Copy genes that are within the 72sp in a new dir "filtered_genes_all72sp"
counter = 0 
for( i in seq( 1:length( name_genes_72sp ) ) ){
  
  counter = counter + 1
  # Copy files
  file.copy( paste( "genes/", as.character( name_genes_72sp[i] ), "/", sep = "" ),
             paste( "filtered_genes_all72sp/", sep = ""),
             recursive = TRUE )
  cat( "Gene ", counter, " is: ", as.character( name_genes_72sp[i] ), "\n" )
  write( paste( "Gene ", counter, " is: ",
                as.character( name_genes_72sp[i] ), "\n", sep = "" ),
         "out_logs/log_01_copy_filtered_all72sp_genes.txt",
         append = TRUE, sep = "\n" )
}
