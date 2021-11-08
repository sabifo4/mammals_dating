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
wd2 <- gsub( pattern = "02_SeqBayes_S2..*", replacement = "", x = wd )

#----------------#
# LOAD LIBRARIES #
#----------------#
library( sn )

#-------------------------#
# LOAD IN-HOUSE FUNCTIONS #
#-------------------------#
source( "../../../../../../src/Plots_check_dists.R" )

#------------------------#
# DEFINE ST CALIBRATIONS #
#------------------------#
# Fossil calibration information used.
# Node 432:  ST ( 1.6419, 0.4248, 12.6518, 1714.5649 ) --> MAMMALIA
# Node 433:  ST ( 0.5771, 0.0070, -0.6339, 7.5086 )    --> ARTIODACTYLA
# Node 434:  ST ( 0.5540, 0.0077, -0.9645, 7.3427 )    --> ARTIOFABULA
# Node 435:  ST ( 0.5070, 0.0092, -1.5132, 6.2939 )    --> CETRUMINANTIA
# Node 440:  ST ( 0.1617, 0.0311, 5.8630, 624.1483 )   --> BOVIDAE
# Node 446:  ST ( 0.0686, 0.0172, 3.2247, 169.6716 )   --> O.aries-C.hircus
#
# 0. Set dir name and get names of the calibrated nodes -- manually done
#    after checking excel file
subtree <- "L.cetartiodactyla"
ST.calib.nodes.names <- c( "MAMMALIA", "ARTIODACTYLA",
                           "ARTIOFABULA", "CETRUMINANTIA",
                           "BOVIDAE", "O.aries-C.hircus"
)

# Create dir to save RData files that will be later used 
# to plot final trees
if( dir.exists( "outRdata" ) != TRUE ){
  dir.create( "outRdata" )
}
save( ST.calib.nodes.names, file = "outRdata/ST.calib.nodes.names.RData" )

# 1. Use in-house function `define_ST` to extract info from `out.txt` output 
#    file by MCMCtree and generate objects with ST info 
ST.dists <- define_ST( log_ST = "00_Prior_onlyST/log.mammals.mcmc1.txt",
                       ST.fitted.dists.RData = paste( wd2, "01_SeqBayes_S1/03_Fit_ST_to_posteriors/00_fitST/Rdata/ST.fitted.dists.RData",
                                                      sep = "" ),
                       ST.calib.nodes.names = ST.calib.nodes.names )
attach( ST.dists )
save( ST.calib.nodes, file = "outRdata/ST.calib.nodes.RData" )

#---------------------------------#
# LOAD RESULTS WITH 72sp DATA SET #
#---------------------------------#==#
# A. Posterior ages (main analysis)  #
#====================================#
load( paste( wd2, "01_SeqBayes_S1/03_Fit_ST_to_posteriors/00_fitST/Rdata/post.divtimes.RData", sep = "" ) ) 
# Rename object and delete "divtimes"
divtimes.post72sp <- divtimes
remove( divtimes )

mean_est_divt.post72sp <- apply( X = divtimes.post72sp, MARGIN = 2, FUN = mean )
quant_est_divt.post72sp<- apply( X = divtimes.post72sp, MARGIN = 2, FUN = quantile, probs = c(0.025,0.975) )
quant_est_divt.post72sp<- t( quant_est_divt.post72sp)
all.equal( quant_est_divt.post72sp[1,], quantile( divtimes.post72sp[,1], probs = c( 0.025, 0.975 ) ) )
var_est_divt.post72sp <- apply( X = divtimes.post72sp, MARGIN = 2, FUN = var )
sd_est_divt.post72sp  <- sqrt( var_est_divt.post72sp )

#-------#
# PRIOR #
#-------#===================================#
# B. Load mcmc.txt from prior with ONLY ST  #
#===========================================#
# Load files and generate objects with sum stats using in-house function `load_data`
onlyST <- load_data( mcmc1 = "00_Prior_onlyST/mcmc1/mcmc.txt", mcmc2 = "00_Prior_onlyST/mcmc2/mcmc.txt")


#-------#
# PRIOR #
#-------#==================================================#
# C. Load mcmc.txt from prior with both ST and soft bounds #
#==========================================================#
# Load files and generate objects with sum stats using in-house function `load_data`
STnSB <- load_data( mcmc1 = "01_Prior_SBandST/mcmc1/mcmc.txt", mcmc2 = "01_Prior_SBandST/mcmc2/mcmc.txt")
var.STnSB <- apply( X = STnSB$divt, MARGIN = 2, FUN = var )
sd.STnSB  <- sqrt( var.STnSB )

#-------#
# PRIOR #
#-------#=============================================================================#
# D. Load mcmc.txt from prior with both ST and soft bounds after tweaking ST          #
#=====================================================================================#
# Load files and generate objects with sum stats using in-house function `load_data`
STnSB.2 <- load_data( mcmc1 = "02_Prior_SBandST_tweak1/mcmc1/mcmc.txt", mcmc2 = "02_Prior_SBandST_tweak1/mcmc2/mcmc.txt" )
var.STnSB.2 <- apply( X = STnSB.2$divt, MARGIN = 2, FUN = var )
sd.STnSB.2  <- sqrt( var.STnSB.2 )

#-------#
# PRIOR #
#-------#========================================================================#
# E. Load mcmc.txt from prior with both ST and soft bounds after tweaking ST     #
#================================================================================#
# Load files and generate objects with sum stats using in-house function `load_data`
STnSB.3 <- load_data( mcmc1 = "03_Prior_SBandST_tweak2/mcmc1/mcmc.txt", mcmc2 = "03_Prior_SBandST_tweak2/mcmc2/mcmc.txt" )
var.STnSB.3 <- apply( X = STnSB.3$divt, MARGIN = 2, FUN = var )
sd.STnSB.3  <- sqrt( var.STnSB.3 )

#-------#
# PRIOR #
#-------#======================================================================#
# F. Load mcmc.txt from prior with both ST and soft bounds after 3nd tweaking  #
#==============================================================================#
# Load files and generate objects with sum stats using in-house function `load_data`
# STnSB.4 <- load_data( mcmc1 = "04_Prior_SBandST_tweak3/mcmc1/mcmc.txt", mcmc2 = "04_Prior_SBandST_tweak3/mcmc2/mcmc.txt" )
# var.STnSB.4 <- apply( X = STnSB.4$divt, MARGIN = 2, FUN = var )
# sd.STnSB.4  <- sqrt( var.STnSB.4 )


#-------#
# PRIOR #
#-------#=====================================================================#
# G. Load mcmc.txt from prior with both ST and soft bounds after 4th tweaking #
#=============================================================================#
# Load files and generate objects with sum stats using in-house function `load_data`
# STnSB.5 <- load_data( mcmc1 = "05_Prior_SBandST_tweak4/mcmc1/mcmc.txt", mcmc2 = "05_Prior_SBandST_tweak4/mcmc2/mcmc.txt" )
# var.STnSB.5 <- apply( X = STnSB.5$divt, MARGIN = 2, FUN = var )
# sd.STnSB.5  <- sqrt( var.STnSB.5 )

#-------#
# PRIOR #
#-------#=====================================================================#
# G. Load mcmc.txt from prior with both ST and soft bounds after 5th tweaking #
#=============================================================================#
# STnSB.6 <- load_data( mcmc1 = "06_Prior_SBandST_tweak5/mcmc1/mcmc.txt", mcmc2 = "06_Prior_SBandST_tweak5/mcmc2/mcmc.txt" )
# var.STnSB.6 <- apply( X = STnSB.6$divt, MARGIN = 2, FUN = var )
# sd.STnSB.6  <- sqrt( var.STnSB.6 )

#---------------------#
# SET PLOT PARAMETERS #
#---------------------#
# Find indexes for Euarchonta post tree and 72sp post tree
ind.post.cal <- ind.post.72sp <- rep( 0, length( ST.calib.nodes ) )
names( ind.post.72sp ) <- names( ind.post.cal ) <- names( ST.calib.nodes ) 

for( i in 1:length( ind.post.cal ) ){
  ind.post.cal[i]  <- which( colnames( onlyST$divt ) == names( ST.calibs )[i] )
  ind.post.72sp[i] <- which( colnames( divtimes.post72sp ) == nodes.72sp[i] )
}
# Set nrows and ncol for plot
p.nrow  <- 2 
p.ncol  <- 3

#-------#
# PLOTS #
#-------#
# A1. Plot densities for ST calibrated nodes when only ST calibrations added
plot_func( post.cal = onlyST$divt, ind.post.cal = ind.post.cal,
           ST.cal = ST.calibs, ST.cal.nodes = names( ST.calibs ),
           post.72sp = divtimes.post72sp, ind.post.72sp = ind.post.72sp,
           names.calibs = ST.calib.nodes.names,
           p.nrow = p.nrow, p.ncol = p.ncol, out = TRUE, legend = TRUE,
           outname = paste( "00_Only_ST_", subtree, sep = "" ) )

# A2. Plot densities for ST calibrated nodes when only ST calibrations added (both runs)
plot_func_mcmc( post.cal.1 = onlyST$divt1, post.cal.2 = onlyST$divt2,
                ind.post.cal = ind.post.cal, ST.cal = ST.calibs,
                ST.cal.nodes = names( ST.calibs ),
                post.72sp = divtimes.post72sp, ind.post.72sp = ind.post.72sp,
                names.calibs = ST.calib.nodes.names,
                p.nrow = p.nrow, p.ncol = p.ncol, out = TRUE,
                legend = TRUE,
                outname = paste( "00_Only_ST_", subtree, "_MCMCruns", sep ="" ) )

# B1. Plot densities for ST calibrated nodes when both ST and SB calibrations added
plot_func( post.cal = STnSB$divt, ind.post.cal = ind.post.cal,
           ST.cal = ST.calibs, ST.cal.nodes = names( ST.calibs ),
           post.72sp = divtimes.post72sp, ind.post.72sp = ind.post.72sp,
           names.calibs = ST.calib.nodes.names,
           p.nrow = p.nrow, p.ncol = p.ncol, out = TRUE, legend = TRUE,
           outname = paste( "01_SBnST_", subtree, sep = "" ) )

# B2. Plot densities for ST calibrated nodes when both ST and SB calibrations added (both runs)
plot_func_mcmc( post.cal.1 = STnSB$divt1, post.cal.2 = STnSB$divt2,
                ind.post.cal = ind.post.cal, ST.cal = ST.calibs,
                ST.cal.nodes = names( ST.calibs ),
                post.72sp = divtimes.post72sp, ind.post.72sp = ind.post.72sp,
                names.calibs = ST.calib.nodes.names,
                p.nrow = p.nrow, p.ncol = p.ncol, out = TRUE,
                legend = TRUE,
                outname = paste( "01_SBnST_", subtree, "_MCMCruns", sep ="" ) )

# C1. Plot densities for ST calibrated nodes when both ST and SB calibrations added after tweaking
plot_func( post.cal = STnSB.2$divt, ind.post.cal = ind.post.cal,
           ST.cal = ST.calibs, ST.cal.nodes = names( ST.calibs ),
           post.72sp = divtimes.post72sp, ind.post.72sp = ind.post.72sp,
           names.calibs = ST.calib.nodes.names,
           p.nrow = p.nrow, p.ncol = p.ncol, out = TRUE, legend = TRUE,
           outname = paste( "02_SBnST_tweak1_", subtree, sep = "" ) )

# C2. Plot densities for ST calibrated nodes when both ST and SB calibrations added (both runs) after tweaking
plot_func_mcmc( post.cal.1 = STnSB.2$divt1, post.cal.2 = STnSB.2$divt2,
                ind.post.cal = ind.post.cal, ST.cal = ST.calibs,
                ST.cal.nodes = names( ST.calibs ),
                post.72sp = divtimes.post72sp, ind.post.72sp = ind.post.72sp,
                names.calibs = ST.calib.nodes.names,
                p.nrow = p.nrow, p.ncol = p.ncol, out = TRUE,
                legend = TRUE,
                outname = paste( "02_SBnST_tweak1_", subtree, "_MCMCruns", sep ="" ) )

# D1. Plot densities for ST calibrated nodes when both ST and SB calibrations added after tweaking
plot_func( post.cal = STnSB.3$divt, ind.post.cal = ind.post.cal,
           ST.cal = ST.calibs, ST.cal.nodes = names( ST.calibs ),
           post.72sp = divtimes.post72sp, ind.post.72sp = ind.post.72sp,
           names.calibs = ST.calib.nodes.names,
           p.nrow = p.nrow, p.ncol = p.ncol, out = TRUE, legend = TRUE,
           outname = paste( "03_SBandSTtweak2_", subtree, sep = "" ) )

# D2. Plot densities for ST calibrated nodes when both ST and SB calibrations added (both runs) after tweaking
plot_func_mcmc( post.cal.1 = STnSB.3$divt1, post.cal.2 = STnSB.3$divt2,
                ind.post.cal = ind.post.cal, ST.cal = ST.calibs,
                ST.cal.nodes = names( ST.calibs ),
                post.72sp = divtimes.post72sp, ind.post.72sp = ind.post.72sp,
                names.calibs = ST.calib.nodes.names,
                p.nrow = p.nrow, p.ncol = p.ncol, out = TRUE,
                legend = TRUE,
                outname = paste( "03_SBandSTtweak2_", subtree, "_MCMCruns", sep ="" ) )

# # E1. Plot densities for ST calibrated nodes when both ST and SB calibrations added after tweaking
# plot_func( post.cal = STnSB.4$divt, ind.post.cal = ind.post.cal,
#            ST.cal = ST.calibs, ST.cal.nodes = names( ST.calibs ),
#            post.72sp = divtimes.post72sp, ind.post.72sp = ind.post.72sp,
#            names.calibs = ST.calib.nodes.names,
#            p.nrow = p.nrow, p.ncol = p.ncol, out = TRUE, legend = TRUE,
#            outname = paste( "04_SBandSTtweak3_", subtree, sep = "" ) )
# 
# # E2. Plot densities for ST calibrated nodes when both ST and SB calibrations added (both runs) after tweaking
# plot_func_mcmc( post.cal.1 = STnSB.4$divt1, post.cal.2 = STnSB.4$divt2,
#                 ind.post.cal = ind.post.cal, ST.cal = ST.calibs,
#                 ST.cal.nodes = names( ST.calibs ),
#                 post.72sp = divtimes.post72sp, ind.post.72sp = ind.post.72sp,
#                 names.calibs = ST.calib.nodes.names,
#                 p.nrow = p.nrow, p.ncol = p.ncol, out = TRUE,
#                 legend = TRUE,
#                 outname = paste( "04_SBandSTtweak3_", subtree, "_MCMCruns", sep ="" ) )
# 
# # F1. Plot densities for ST calibrated nodes when both ST and SB calibrations added after tweaking
# plot_func( post.cal = STnSB.5$divt, ind.post.cal = ind.post.cal,
#            ST.cal = ST.calibs, ST.cal.nodes = names( ST.calibs ),
#            post.72sp = divtimes.post72sp, ind.post.72sp = ind.post.72sp,
#            names.calibs = ST.calib.nodes.names,
#            p.nrow = p.nrow, p.ncol = p.ncol, out = TRUE, legend = TRUE,
#            outname = paste( "05_SBandSTtweak4_", subtree, sep = "" ) )
# 
# # F2. Plot densities for ST calibrated nodes when both ST and SB calibrations added (both runs) after tweaking
# plot_func_mcmc( post.cal.1 = STnSB.5$divt1, post.cal.2 = STnSB.5$divt2,
#                 ind.post.cal = ind.post.cal, ST.cal = ST.calibs,
#                 ST.cal.nodes = names( ST.calibs ),
#                 post.72sp = divtimes.post72sp, ind.post.72sp = ind.post.72sp,
#                 names.calibs = ST.calib.nodes.names,
#                 p.nrow = p.nrow, p.ncol = p.ncol, out = TRUE,
#                 legend = TRUE,
#                 outname = paste( "05_SBandSTtweak4_", subtree, "_MCMCruns", sep ="" ) )
# 
# # G1. Plot densities for ST calibrated nodes when both ST and SB calibrations added after tweaking
# plot_func( post.cal = STnSB.6$divt, ind.post.cal = ind.post.cal,
#            ST.cal = ST.calibs, ST.cal.nodes = names( ST.calibs ),
#            post.72sp = divtimes.post72sp, ind.post.72sp = ind.post.72sp,
#            names.calibs = ST.calib.nodes.names,
#            p.nrow = p.nrow, p.ncol = p.ncol, out = TRUE, legend = TRUE,
#            outname = paste( "06_SBandSTtweak5_", subtree, sep = "" ) )
# 
# # G2. Plot densities for ST calibrated nodes when both ST and SB calibrations added (both runs) after tweaking
# plot_func_mcmc( post.cal.1 = STnSB.6$divt1, post.cal.2 = STnSB.6$divt2,
#                 ind.post.cal = ind.post.cal, ST.cal = ST.calibs,
#                 ST.cal.nodes = names( ST.calibs ),
#                 post.72sp = divtimes.post72sp, ind.post.72sp = ind.post.72sp,
#                 names.calibs = ST.calib.nodes.names,
#                 p.nrow = p.nrow, p.ncol = p.ncol, out = TRUE,
#                 legend = TRUE,
#                 outname = paste( "06_SBandSTtweak5_", subtree, "_MCMCruns", sep ="" ) )

###---- END PLOTS -----


#=================#
# PLOT MEAN/QUANT #
#=================#
# Extract mean/quant for results when both SB and ST are added
pdf( file = paste( "01_SBnST_", subtree, "_meanquant.pdf", sep = "" ),
     paper = "a4", width = 1024, height = 768 )
curr_plot <- dev.cur()
png( filename = paste( "01_SBnST_", subtree, "_meanquant.png", sep = "" ),
     width = 1024, height = 768 )
dev.control( "enable" )
par( mfrow = c( 1,3 ) )

## MEANS
x.mean <- mean_est_divt.post72sp[ind.post.72sp]*100
y.mean <- STnSB$mean_divt[ind.post.cal]*100
col.mean <- rep( "black", length( ST.calib.nodes ) )
cex.lab <- rep( 0.3, length( ST.calib.nodes ) )
dub.pos <- c(5,6)
col.mean[dub.pos] <- "red"
plot( x = x.mean, y = y.mean, col = col.mean, pch = 20,
      xlab = expression( 'Posterior mean time '*italic(t)*' (Ma)|72sp - data1' ),
      ylab = expression( 'Posterior mean time '*italic(t)*' (Ma)|431sp - artiodactyla_data2' ),
      xlim = c( min( x.mean )-5, max( x.mean )+ 90 ) )
coords.m <-  list( x = x.mean + 60, y = y.mean )
text( coords.m$x, coords.m$y, labels = ST.calib.nodes.names, cex= cex.lab, col = col.mean )
abline( lm( y.mean ~ 0 + x.mean ), col="darkgray", lty = 2 )
lm.m <- lm( y.mean ~ 0 + x.mean )
summary( lm.m )
text( x = 100, y = 180,
      labels = expression( atop( ''*italic( R^2 )*' =  0.9998' ) ) )
text( x = 100, y = 175,
      labels = expression( atop( ''*italic( w )*' = '*italic( 0.997349 * t )*'' ) ) )

## QUANTILE 2.5 
x.qlow <- quant_est_divt.post72sp[ind.post.72sp,1]*100
y.qlow <- STnSB$quant_divt[ind.post.cal,1]*100
plot( x = x.qlow, y = y.qlow, col = col.mean, pch = 20,
      main = subtree, 
      xlab = expression( 'Quant. 2.5% post. mean time '*italic(t)*' (Ma)|72sp - data1' ),
      ylab = expression( 'Quant. 2.5% post. mean time '*italic(t)*' (Ma)|431sp - artiodactyla_data2' ),
      xlim = c( min( x.qlow )-5, max( x.qlow )+ 90 ) )
coords.ql <-  list( x = x.qlow + 60, y = y.qlow )
text( coords.ql$x, coords.ql$y, labels = ST.calib.nodes.names, cex= 0.8, col = col.mean )
abline( lm( y.qlow ~ 0 + x.qlow ), col="darkgray", lty = 2 )
lm.qlow <- lm( y.qlow ~ 0 + x.qlow )
summary( lm.qlow )
text( x = 100, y = 150,
      labels = expression( atop( ''*italic( R^2 )*' = 0.9999' ) ) )
text( x = 100, y = 145,
      labels = expression( atop( ''*italic( w )*' = '*italic(  1.001158 * t )*'' ) ) )

## QUANTILE 97.5
x.qup <- quant_est_divt.post72sp[ind.post.72sp,2]*100
y.qup <- STnSB$quant_divt[ind.post.cal,2]*100
plot( x = x.qup, y = y.qup, col = col.mean, pch = 20,
      xlab = expression( 'Quant. 97.5% post. mean time '*italic(t)*' (Ma)|72sp - data1' ),
      ylab = expression( 'Quant. 97.5% post. mean time '*italic(t)*' (Ma)|431sp - artiodactyla_data2' ),
      xlim = c( min( x.qup )-5, max( x.qup )+ 115 ) )
coords.qu <-  list( x = x.qup + 75, y = y.qup )
text( coords.qu$x, coords.qu$y, labels = ST.calib.nodes.names, cex= 0.8, col = col.mean )
abline( lm( y.qup ~ 0 + x.qup ), col="darkgray", lty = 2 )
lm.qup <- lm( y.qup ~ 0 + x.qup )
summary( lm.qup )
text( x = 110, y = 200,
      labels = expression( atop( ''*italic( R^2 )*' = 0.9997' ) ) )
text( x = 110, y = 195,
      labels = expression( atop( ''*italic( w )*' = '*italic(  1.035801 * t )*'' ) ) )

# Stop graphics !
dev.copy( which = curr_plot )
dev.off()
dev.off()

# Check deviations
#(m_new - m_post) / m_post
mean.dev <- ( y.mean - x.mean )/x.mean
qup.dev  <- ( y.qup - x.qup )/x.qup
qlow.dev <- ( y.qlow - x.qlow )/x.qlow
var.dev  <- ( var.STnSB[ind.post.cal] - var_est_divt.post72sp[ind.post.72sp] )/ var_est_divt.post72sp[ind.post.72sp]
sd.dev   <- ( sd.STnSB[ind.post.cal] - sd_est_divt.post72sp[ind.post.72sp] )/ sd_est_divt.post72sp[ind.post.72sp]
sum.stats <- matrix( nrow = length( y.mean ), ncol = 5 )
sum.stats[,1] <- mean.dev
sum.stats[,2] <- qlow.dev
sum.stats[,3] <- qup.dev
sum.stats[,4] <- var.dev
sum.stats[,5] <- sd.dev
colnames( sum.stats ) <- c( "mean", "q-2.5%", "q-97.5%", "var", "sd" )
rownames( sum.stats ) <- paste( ST.calib.nodes, names( ST.calib.nodes), sep = "-")
sum.stats.ord <- round( sum.stats*100, digits = 1 )
tmp.col1 <- which( abs( sum.stats.ord[,1] ) > 5 )
tmp.col2 <- which( abs( sum.stats.ord[,2] ) > 5 )
tmp.col3 <- which( abs( sum.stats.ord[,3] ) > 5 )
tmp.dub  <- sort( unique( c( tmp.col1, tmp.col2, tmp.col3 ) ) )
rownames( sum.stats )[tmp.dub]
write.table( x = sum.stats.ord, file = paste( "01_SBnST_", subtree, "_sumstats.tsv", sep = "" ),
             sep = "\t", quote = FALSE )


#================#
# TWEAK ST DISTS #
#================#
# Values to add to ST for next tweaking (location ~ mean)
#   100                         ---> x.mean[dub.pos]
#   sum.stats.ord[dub.pos,1]    ---> add.to.STxi = ?
#
# A) E.g., 11.9% (mean deviation) of 18.64778 (My, x.mean_Bovidae) --> 2.219086
#    Let's subtract that deviation from x_i (ST.calibs[[ 5 ]][1]) after dividing into 100
#       0.1617 - 0.02219086 = 0.1395091
#===================================================================================================#
# Values to add to ST for next tweaking (scale ~ sd)
#   100                       ---> sd_est_divt.post72sp[ind.post.72sp[dub.pos]]
#   sum.stats.ord[dub.pos,5]  ---> add.to.STwi = ?
#
# A) E.g., 34.6% (mean deviation) of 0.01888672 (sd) --> 0.006534805
#    Let's subtract that deviation from w_i (ST.calibs[[ 5 ]][2])
#       0.0311 - 0.006534805 = 0.02456519
dub.pos              <- c( 5, 6 ) 
add.to.STxi          <- - sum.stats.ord[dub.pos,1] * x.mean[dub.pos] / 100
add.to.STwi          <- - sum.stats.ord[dub.pos,5] * sd_est_divt.post72sp[ind.post.72sp[dub.pos]] / 100
names( add.to.STwi ) <- rownames( sum.stats )[dub.pos]
node.labs.1     <- gsub( pattern = "-..*", replacement = "", x = names( add.to.STxi ) )
node.labs.2     <- gsub( pattern = "t_n", replacement = "", x = names( ST.calibs ) )
ST.calibs.tweak <- ST.calibs
tmp.tab             <- matrix( 0, nrow = length( add.to.STxi ), ncol = 4 )
rownames( tmp.tab ) <- names( add.to.STxi )
for( i in dub.pos ){
  
  # Match nodes and then replace accordingly in "xi" an "wi" for ST calibs
  index <- which( node.labs.1 == node.labs.2[i] )
  # x_i
  ST.calibs.tweak[[ i ]][1] <- ST.calibs[[ i ]][1] + ( add.to.STxi[ index ]/100 )
  cat( "Prior ST: ", ST.calibs[[ i ]][1], "| Mean post ST:",  y.mean[i]/100, "| Mean post ST with 72sp: ", x.mean[i]/100, "\n" )
  cat( "We will add", add.to.STxi[index], " to the xi parameter of the prior ST calib\n" )
  cat( ST.calibs[[ i ]][1], "+", add.to.STxi[index]/100, "=", ST.calibs[[ i ]][1] + ( add.to.STxi[ index ]/100 ), "\n\n" )
  # w_i
  ST.calibs.tweak[[ i ]][2] <- ST.calibs[[ i ]][2] + ( add.to.STwi[ index ] )
  cat( "Prior ST: ", ST.calibs[[ i ]][1], "| Mean post ST:",  sd.STnSB[ind.post.cal[i]],
       "| Mean post ST with 72sp: ", sd_est_divt.post72sp[ind.post.72sp[i]], "\n" )
  cat( "We will add", add.to.STwi[index], " to the w_i parameter of the prior ST calib\n" )
  cat( ST.calibs[[ i ]][2], "+", add.to.STwi[index], "=", ST.calibs[[ i ]][2] + ( add.to.STwi[ index ] ), "\n\n" )
  tmp.tab[index,] <- ST.calibs.tweak[[ i ]]
  
}

write.table( x = tmp.tab, file = paste( "outRdata/", subtree, "_tweaked_ST_1.tsv", sep = "" ),
             quote = FALSE, col.names = FALSE, sep = "\t" )

#=================================#

#=================================#
# PLOT MEAN/QUANT - STnSB2 OBJECT #
#=================================#
# Extract mean/quant for results when both SB and ST are added
pdf( file = paste( "02_SBnST_tweak1_", subtree, "_meanquant.pdf", sep = "" ),
     paper = "a4", width = 1024, height = 768 )
curr_plot <- dev.cur()
png( filename = paste( "02_SBnST_tweak1_", subtree, "_meanquant.png", sep = "" ),
     width = 1024, height = 768 )
dev.control( "enable" )
par( mfrow = c( 1,3 ) )

## MEANS
x.mean.2 <- mean_est_divt.post72sp[ind.post.72sp]*100
y.mean.2 <- STnSB.2$mean_divt[ind.post.cal]*100
col.mean <- rep( "black", length( ST.calib.nodes ) )
cex.lab <- rep( 0.3, length( ST.calib.nodes ) )
dub.pos.2 <- 5
col.mean[dub.pos.2] <- "red"
cex.lab[dub.pos.2] <- 0.7
plot( x = x.mean.2, y = y.mean.2, pch = 20, col = col.mean,
      xlab = expression( 'Posterior mean time '*italic(t)*' (Ma)|72sp - data1' ),
      ylab = expression( 'Posterior mean time '*italic(t)*' (Ma)|431sp - artiodactyla_data2' ),
      xlim = c( min( x.mean.2 )-5, max( x.mean.2 )+ 90 ) )
coords.m <-  list( x = x.mean.2 + 60, y = y.mean.2 )
text( coords.m$x, coords.m$y, labels = ST.calib.nodes.names, cex= cex.lab, col = col.mean )
abline( lm( y.mean.2 ~ 0 + x.mean.2 ), col="darkgray", lty = 2 )
lm.m <- lm( y.mean.2 ~ 0 + x.mean.2 )
summary( lm.m )
text( x = 100, y = 180,
      labels = expression( atop( ''*italic( R^2 )*' =  1' ) ) )
text( x = 100, y = 175,
      labels = expression( atop( ''*italic( w )*' = '*italic( 0.996535 * t )*'' ) ) )

## QUANTILE 2.5 
x.qlow.2 <- quant_est_divt.post72sp[ind.post.72sp,1]*100
y.qlow.2 <- STnSB.2$quant_divt[ind.post.cal,1]*100
plot( x = x.qlow.2, y = y.qlow.2, col = col.mean, pch = 20,
      main = subtree, 
      xlab = expression( 'Quant. 2.5% post. mean time '*italic(t)*' (Ma)|72sp - data1' ),
      ylab = expression( 'Quant. 2.5% post. mean time '*italic(t)*' (Ma)|431sp - artiodactyla_data2' ),
      xlim = c( min( x.qlow.2 )-5, max( x.qlow.2 )+ 90 ) )
coords.ql <-  list( x = x.qlow.2 + 60, y = y.qlow.2 )
text( coords.ql$x, coords.ql$y, labels = ST.calib.nodes.names, cex= cex.lab, col = col.mean )
abline( lm( y.qlow.2 ~ 0 + x.qlow.2 ), col="darkgray", lty = 2 )
lm.qlow <- lm( y.qlow.2 ~ 0 + x.qlow.2 )
summary( lm.qlow )
text( x = 100, y = 150,
      labels = expression( atop( ''*italic( R^2 )*' = 0.9999' ) ) )
text( x = 100, y = 145,
      labels = expression( atop( ''*italic( w )*' = '*italic(  1.000874 * t )*'' ) ) )

## QUANTILE 97.5
x.qup.2 <- quant_est_divt.post72sp[ind.post.72sp,2]*100
y.qup.2 <- STnSB.2$quant_divt[ind.post.cal,2]*100
plot( x = x.qup.2, y = y.qup.2, col = col.mean, pch = 20,
      xlab = expression( 'Quant. 97.5% post. mean time '*italic(t)*' (Ma)|72sp - data1' ),
      ylab = expression( 'Quant. 97.5% post. mean time '*italic(t)*' (Ma)|431sp - artiodactyla_data2' ),
      xlim = c( min( x.qup.2 )-5, max( x.qup.2 )+ 100 ) )
coords.qu <-  list( x = x.qup.2 + 70, y = y.qup.2 )
text( coords.qu$x, coords.qu$y, labels = ST.calib.nodes.names, cex= cex.lab, col = col.mean )
abline( lm( y.qup.2 ~ 0 + x.qup.2 ), col="darkgray", lty = 2 )
lm.qup <- lm( y.qup.2 ~ 0 + x.qup.2 )
summary( lm.qup )
text( x = 110, y = 200,
      labels = expression( atop( ''*italic( R^2 )*' = 0.9997' ) ) )
text( x = 110, y = 195,
      labels = expression( atop( ''*italic( w )*' = '*italic(  1.036259 * t )*'' ) ) )

# Stop graphics !
dev.copy( which = curr_plot )
dev.off()
dev.off()

# Check deviations
#  Error in the calibration: (m_new - m_post) / m_post
#   x.mean.2 (m_post) <- mean_est_divt.post72sp[ind.post.72sp]*100
#   y.mean.2 (m_new)  <- STnSB.2$mean_divt[ind.post.cal]*100
#
#  Error in the calibration: (qup_new - qup_post) / qup_post
#   x.qup.2 (m_post)  <- quant_est_divt.post72sp[ind.post.72sp,2]*100
#   y.qup.2 (m_new    <- STnSB.2$quant_divt[ind.post.cal,2]*100
#
#  Error in the calibration: (qlow_new - qlow_post) / qlow_post
#   x.qlow.2 (m_post) <- quant_est_divt.post72sp[ind.post.72sp,1]*100
#   y.qlow.2 (m_new)  <- STnSB.2$quant_divt[ind.post.cal,1]*100
mean.dev.2  <- ( y.mean.2 - x.mean.2 )/x.mean.2
qup.dev.2   <- ( y.qup.2 - x.qup.2 )/x.qup.2
qlow.dev.2  <- ( y.qlow.2 - x.qlow.2 )/x.qlow.2
var.dev.2   <- ( var.STnSB.2[ind.post.cal] - var_est_divt.post72sp[ind.post.72sp] )/ var_est_divt.post72sp[ind.post.72sp]
sd.dev.2    <- ( sd.STnSB.2[ind.post.cal] - sd_est_divt.post72sp[ind.post.72sp] )/ sd_est_divt.post72sp[ind.post.72sp]
sum.stats.2 <- matrix( nrow = length( y.mean.2 ), ncol = 5 )
sum.stats.2[,1] <- mean.dev.2
sum.stats.2[,2] <- qlow.dev.2
sum.stats.2[,3] <- qup.dev.2
sum.stats.2[,4] <- var.dev.2
sum.stats.2[,5] <- sd.dev.2
colnames( sum.stats.2 ) <- c( "mean", "q-2.5%", "q-97.5%", "var", "sd" )
rownames( sum.stats.2 ) <- paste( ST.calib.nodes, names( ST.calib.nodes), sep = "-")
sum.stats.ord.2 <- round( sum.stats.2*100, digits = 1 )
tmp.col1 <- which( abs( sum.stats.ord.2[,1] ) > 5 )
tmp.col2 <- which( abs( sum.stats.ord.2[,2] ) > 5 )
tmp.col3 <- which( abs( sum.stats.ord.2[,3] ) > 5 )
tmp.dub  <- sort( unique( c( tmp.col1, tmp.col2, tmp.col3 ) ) )
rownames( sum.stats )[tmp.dub]
write.table( x = sum.stats.ord.2, file = paste( "02_SBnSTtweak1_", subtree, "_sumstats.tsv", sep = "" ),
             sep = "\t", quote = FALSE )


#================#
# TWEAK ST DISTS #
#================#
# Values to add to ST for next tweaking (location ~ mean)
#   100                             ---> x.mean.2[dub.pos.2]
#   sum.stats.ord.2[dub.pos.2,1]    ---> add.to.STxi = ?
#
# A) E.g., -5% (mean deviation) of 18.64778 (My, x.mean) --> -0.932389
#    Let's subtract that deviation from x_i (ST.calibs.tweak[[ 5 ]])after dividing into 100
#       0.1395091 - (-0.00932389) = 0.148833
#===================================================================================================#
# Values to add to ST for next tweaking (scale ~ sd)
#   100                            ---> sd_est_divt.post72sp[ind.post.72sp[dub.pos.2]]
#   sum.stats.ord.2[dub.pos.2,5]  ---> add.to.STwi = ?
#
# A) E.g., 0.2% (mean deviation) of 0.01888672  (sd) --> 3.777344e-05
#    Let's subtract that deviation from w_i (ST.calibs.tweak)
#       0.0245652 - 3.777344e-05 = 0.02452743
add.to.STxi.2          <- - sum.stats.ord.2[dub.pos.2,1] * x.mean.2[dub.pos.2] / 100
add.to.STwi.2          <- - sum.stats.ord.2[dub.pos.2,5] * sd_est_divt.post72sp[ind.post.72sp[dub.pos.2]] / 100
names( add.to.STwi.2 ) <- rownames( sum.stats.2 )[dub.pos.2]
# Commented as they have been run above
node.labs.1     <- gsub( pattern = "-..*", replacement = "", x = rownames( sum.stats.ord.2 )[dub.pos.2] )
# node.labs.2     <- gsub( pattern = "t_n", replacement = "", x = names( ST.calibs ) )
ST.calibs.tweak.2   <- ST.calibs.tweak
tmp.tab             <- matrix( 0, nrow = length( add.to.STwi.2 ), ncol = 4 )
rownames( tmp.tab ) <- names( add.to.STwi.2 )
for( i in dub.pos.2 ){
  
  # Match nodes and then replace accordingly in "xi" an "wi" for ST calibs
  index <- which( node.labs.1 == node.labs.2[i] )
  # x_i
  ST.calibs.tweak.2[[ i ]][1] <- ST.calibs.tweak[[ i ]][1] + ( add.to.STxi.2[ index ]/100 )
  cat( "Prior ST: ", ST.calibs.tweak[[ i ]][1], "| Mean post ST:",  y.mean.2[i]/100, "| Mean post ST with 72sp: ", x.mean.2[i]/100, "\n" )
  cat( "We will add", add.to.STxi.2[index], " to the xi parameter of the prior ST calib\n" )
  cat( ST.calibs.tweak[[ i ]][1], "+", add.to.STxi.2[index]/100, "=", ST.calibs.tweak[[ i ]][1] + ( add.to.STxi.2[ index ]/100 ), "\n\n" )
  # w_i
  ST.calibs.tweak.2[[ i ]][2] <- ST.calibs.tweak[[ i ]][2] + ( add.to.STwi.2[ index ] )
  cat( "Prior ST: ", ST.calibs.tweak[[ i ]][2], "| Mean post SD.ST:",  sd.STnSB.2[ind.post.cal[i]],
       "| Mean post SD.ST with 72sp: ", sd_est_divt.post72sp[ind.post.72sp[i]], "\n" )
  cat( "We will add", add.to.STwi.2[index], " to the w_i parameter of the prior ST calib\n" )
  cat( ST.calibs.tweak[[ i ]][2], "+", add.to.STwi.2[index], "=", ST.calibs.tweak[[ i ]][2] + ( add.to.STwi.2[ index ] ), "\n\n" )
  tmp.tab[index,] <- ST.calibs.tweak.2[[ i ]]
  
}

write.table( x = tmp.tab, file = paste( "outRdata/", subtree, "_tweaked_ST_2.tsv", sep = "" ),
             quote = FALSE, col.names = FALSE, sep = "\t" )

#------------------------------

#==================================================#
# NEW PLOTS AFTER SECOND TWEAKING ST AS DONE ABOVE #
#==================================================#
#=================================#
# PLOT MEAN/QUANT - STnSB3 OBJECT #
#=================================#
# Extract mean/quant for results when both SB and ST are added
pdf( file = paste( "03_SBnST_tweak2_", subtree, "_meanquant.pdf", sep = "" ),
     paper = "a4", width = 1024, height = 768 )
curr_plot <- dev.cur()
png( filename = paste( "03_SBnST_tweak2_", subtree, "_meanquant.png", sep = "" ),
     width = 1024, height = 768 )
dev.control( "enable" )
par( mfrow = c( 1,3 ) )

## MEANS
x.mean.3 <- mean_est_divt.post72sp[ind.post.72sp]*100
y.mean.3 <- STnSB.3$mean_divt[ind.post.cal]*100
col.mean <- rep( "black", length( ST.calib.nodes ) )
cex.lab <- rep( 0.3, length( ST.calib.nodes ) )
dub.pos <- 5
col.mean[dub.pos] <- "red"
cex.lab[dub.pos] <- 0.7
plot( x = x.mean.3, y = y.mean.3, pch = 20, col = col.mean,
      xlab = expression( 'Posterior mean time '*italic(t)*' (Ma)|72sp - data1' ),
      ylab = expression( 'Posterior mean time '*italic(t)*' (Ma)|431sp - artiodactyla_data2' ),
      xlim = c( min( x.mean.3 )-5, max( x.mean.3 )+ 90 ) )
coords.m <-  list( x = x.mean.3 + 60, y = y.mean.3 )
text( coords.m$x, coords.m$y, labels = ST.calib.nodes.names, cex= cex.lab, col = col.mean )
abline( lm( y.mean.3 ~ 0 + x.mean.3 ), col="darkgray", lty = 2 )
lm.m <- lm( y.mean.3 ~ 0 + x.mean.3 )
summary( lm.m )
text( x = 100, y = 180,
      labels = expression( atop( ''*italic( R^2 )*' =  1 ' ) ) )
text( x = 100, y = 175,
      labels = expression( atop( ''*italic( w )*' = '*italic( 0.996103 * t )*'' ) ) )

## QUANTILE 2.5 
x.qlow.3 <- quant_est_divt.post72sp[ind.post.72sp,1]*100
y.qlow.3 <- STnSB.3$quant_divt[ind.post.cal,1]*100
plot( x = x.qlow.3, y = y.qlow.3, col = col.mean, pch = 20,
      main = subtree, 
      xlab = expression( 'Quant. 2.5% post. mean time '*italic(t)*' (Ma)|72sp - data1' ),
      ylab = expression( 'Quant. 2.5% post. mean time '*italic(t)*' (Ma)|431sp - artiodactyla_data2' ),
      xlim = c( min( x.qlow.3 )-5, max( x.qlow.3 )+ 90 ) )
coords.ql <-  list( x = x.qlow.3 + 60, y = y.qlow.3 )
text( coords.ql$x, coords.ql$y, labels = ST.calib.nodes.names, cex= cex.lab, col = col.mean )
abline( lm( y.qlow.3 ~ 0 + x.qlow.3 ), col="darkgray", lty = 2 )
lm.qlow <- lm( y.qlow.3 ~ 0 + x.qlow.3 )
summary( lm.qlow )
text( x = 100, y = 150,
      labels = expression( atop( ''*italic( R^2 )*' = 0.9999' ) ) )
text( x = 100, y = 145,
      labels = expression( atop( ''*italic( w )*' = '*italic(  1.000563 * t )*'' ) ) )

## QUANTILE 97.5
x.qup.3 <- quant_est_divt.post72sp[ind.post.72sp,2]*100
y.qup.3 <- STnSB.3$quant_divt[ind.post.cal,2]*100
plot( x = x.qup.3, y = y.qup.3, col = col.mean, pch = 20,
      xlab = expression( 'Quant. 97.5% post. mean time '*italic(t)*' (Ma)|72sp - data1' ),
      ylab = expression( 'Quant. 97.5% post. mean time '*italic(t)*' (Ma)|431sp - artiodactyla_data2' ),
      xlim = c( min( x.qup.3 )-5, max( x.qup.3 )+ 100 ) )
coords.qu <-  list( x = x.qup.3 + 70, y = y.qup.3 )
text( coords.qu$x, coords.qu$y, labels = ST.calib.nodes.names, cex= cex.lab, col = col.mean )
abline( lm( y.qup.3 ~ 0 + x.qup.3 ), col="darkgray", lty = 2 )
lm.qup <- lm( y.qup.3 ~ 0 + x.qup.3 )
summary( lm.qup )
text( x = 110, y = 200,
      labels = expression( atop( ''*italic( R^2 )*' = 0.9997' ) ) )
text( x = 110, y = 195,
      labels = expression( atop( ''*italic( w )*' = '*italic(  1.036816 * t )*'' ) ) )

# Stop graphics !
dev.copy( which = curr_plot )
dev.off()
dev.off()

# Check deviations
#  Error in the calibration: (m_new - m_post) / m_post
#   x.mean.3 (m_post) <- mean_est_divt.post72sp[ind.post.72sp]*100
#   y.mean.3 (m_new)  <- STnSB.3$mean_divt[ind.post.cal]*100
#
#  Error in the calibration: (qup_new - qup_post) / qup_post
#   x.qup.3 (m_post)  <- quant_est_divt.post72sp[ind.post.72sp,2]*100
#   y.qup.3 (m_new    <- STnSB.3$quant_divt[ind.post.cal,2]*100
#
#  Error in the calibration: (qlow_new - qlow_post) / qlow_post
#   x.qlow.3 (m_post) <- quant_est_divt.post72sp[ind.post.72sp,1]*100
#   y.qlow.3 (m_new)  <- STnSB.3$quant_divt[ind.post.cal,1]*100
mean.dev.3  <- ( y.mean.3 - x.mean.3 )/x.mean.3
qup.dev.3   <- ( y.qup.3 - x.qup.3 )/x.qup.3
qlow.dev.3  <- ( y.qlow.3 - x.qlow.3 )/x.qlow.3
var.dev.3   <- ( var.STnSB.3[ind.post.cal] - var_est_divt.post72sp[ind.post.72sp] )/ var_est_divt.post72sp[ind.post.72sp]
sd.dev.3    <- ( sd.STnSB.3[ind.post.cal] - sd_est_divt.post72sp[ind.post.72sp] )/ sd_est_divt.post72sp[ind.post.72sp]
sum.stats.3 <- matrix( nrow = length( y.mean.3 ), ncol = 5 )
sum.stats.3[,1] <- mean.dev.3
sum.stats.3[,2] <- qlow.dev.3
sum.stats.3[,3] <- qup.dev.3
sum.stats.3[,4] <- var.dev.3
sum.stats.3[,5] <- sd.dev.3
colnames( sum.stats.3 ) <- c( "mean", "q-2.5%", "q-97.5%", "var", "sd" )
rownames( sum.stats.3 ) <- paste( ST.calib.nodes, names( ST.calib.nodes), sep = "-")
sum.stats.ord.3 <- round( sum.stats.3*100, digits = 1 )
tmp.col1 <- which( abs( sum.stats.ord.3[,1] ) > 5 )
tmp.col2 <- which( abs( sum.stats.ord.3[,2] ) > 5 )
tmp.col3 <- which( abs( sum.stats.ord.3[,3] ) > 5 )
tmp.dub  <- sort( unique( c( tmp.col1, tmp.col2, tmp.col3 ) ) )
rownames( sum.stats )[tmp.dub]
write.table( x = sum.stats.ord.3, file = paste( "03_SBnSTtweak2_", subtree, "_sumstats.tsv", sep = "" ),
             sep = "\t", quote = FALSE )

detach( ST.dists )
