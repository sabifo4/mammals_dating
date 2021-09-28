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
source( "../../../../../../../src/Plots_check_dists.R" )

#------------------------#
# DEFINE ST CALIBRATIONS #
#------------------------#
# Which nodes were ST-calibrated when sampling from the prior?
#
# Node 631:  ST ( 1.6419, 0.4248, 12.6518, 1714.5649 )
# Node 632:  ST ( 0.5496, 0.0063, -0.1714, 6.7595 )
# Node 633:  ST ( 0.4732, 0.0175, -1.2816, 13.0445 )
# Node 635:  ST ( 0.3819, 0.0218, -0.6701, 36.8117 )
# Node 636:  ST ( 0.1744, 0.0351, 3.1168, 381.9704 )
# Node 638:  ST ( 0.1279, 0.0288, 3.5319, 404.4165 )
# Node 938:  ST ( 0.0880, 0.0217, 3.7564, 209.0654 )
# Node 1011:  ST ( 0.0889, 0.0241, 3.7938, 547.4454 )


# 0. Set dir name and get names of the calibrated nodes -- manually done
#    after checking excel file  (needs to follow same order as calibs in output file!)
subtree <- "Rod_subtree1"
ST.calib.nodes.names <- c( "Mammalia",
                           "dipodomys_ordi-MYOMORPHA",                
                           "MYOMORPHA",        
                           "nannospalax_galili-MURIDAE",                
                           "MURIDAE",
                           "Peromyscus-Microtus",    
                           "Mauratus-Cgriseus",
                           "Murinae"
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
onlyST <- load_data( mcmc1 = "00_Prior_onlyST/mcmc1/mcmc.txt", mcmc2 = "00_Prior_onlyST/mcmc2/mcmc.txt" )

#-------#
# PRIOR #
#-------#==================================================#
# C. Load mcmc.txt from prior with both ST and soft bounds #
#==========================================================#
# Load files and generate objects with sum stats using in-house function `load_data`
STnSB     <- load_data( mcmc1 = "01_Prior_SBandST/mcmc1/mcmc.txt", mcmc2 = "01_Prior_SBandST/mcmc2/mcmc.txt" )
var.STnSB <- apply( X = STnSB$divt, MARGIN = 2, FUN = var )
sd.STnSB  <- sqrt( var.STnSB )

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
p.ncol  <- 4

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

## ---- END PLOTS ----- ##


#================================#
# PLOT MEAN/QUANT - STnSB OBJECT #
#================================#
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
# dub.pos <- c( 7,8 )
# col.mean[dub.pos] <- "red"
# cex.lab[dub.pos] <- 0.7
plot( x = x.mean, y = y.mean, pch = 20, col = col.mean,
      xlab = expression( 'Posterior mean time '*italic(t)*' (Ma)|72sp - data1' ),
      ylab = expression( 'Posterior mean time '*italic(t)*' (Ma)|630sp - rod_subt1_data2' ),
      xlim = c( min( x.mean )-5, max( x.mean )+ 90 ) )
coords.m <-  list( x = x.mean + 60, y = y.mean )
text( coords.m$x, coords.m$y, labels = ST.calib.nodes.names, cex= cex.lab, col = col.mean )
abline( lm( y.mean ~ 0 + x.mean ), col="darkgray", lty = 2 )
lm.m <- lm( y.mean ~ 0 + x.mean )
summary( lm.m )
text( x = 100, y = 180,
      labels = expression( atop( ''*italic( R^2 )*' =  1' ) ) )
text( x = 100, y = 175,
      labels = expression( atop( ''*italic( w )*' = '*italic( 0.996113 * t )*'' ) ) )

## QUANTILE 2.5 
x.qlow <- quant_est_divt.post72sp[ind.post.72sp,1]*100
y.qlow <- STnSB$quant_divt[ind.post.cal,1]*100
plot( x = x.qlow, y = y.qlow, col = col.mean, pch = 20,
      main = subtree, 
      xlab = expression( 'Quant. 2.5% post. mean time '*italic(t)*' (Ma)|72sp - data1' ),
      ylab = expression( 'Quant. 2.5% post. mean time '*italic(t)*' (Ma)|630sp - rod_subt1_data2' ),
      xlim = c( min( x.qlow )-5, max( x.qlow )+ 90 ) )
coords.ql <-  list( x = x.qlow + 60, y = y.qlow )
text( coords.ql$x, coords.ql$y, labels = ST.calib.nodes.names, cex= cex.lab, col = col.mean )
abline( lm( y.qlow ~ 0 + x.qlow ), col="darkgray", lty = 2 )
lm.qlow <- lm( y.qlow ~ 0 + x.qlow )
summary( lm.qlow )
text( x = 100, y = 150,
      labels = expression( atop( ''*italic( R^2 )*' = 1' ) ) )
text( x = 100, y = 145,
      labels = expression( atop( ''*italic( w )*' = '*italic(  0.999748 * t )*'' ) ) )

## QUANTILE 97.5
x.qup <- quant_est_divt.post72sp[ind.post.72sp,2]*100
y.qup <- STnSB$quant_divt[ind.post.cal,2]*100
plot( x = x.qup, y = y.qup, col = col.mean, pch = 20,
      xlab = expression( 'Quant. 97.5% post. mean time '*italic(t)*' (Ma)|72sp - data1' ),
      ylab = expression( 'Quant. 97.5% post. mean time '*italic(t)*' (Ma)|630sp - rod_subt1_data2' ),
      xlim = c( min( x.qup )-5, max( x.qup )+ 100 ) )
coords.qu <-  list( x = x.qup + 70, y = y.qup )
text( coords.qu$x, coords.qu$y, labels = ST.calib.nodes.names, cex= cex.lab, col = col.mean )
abline( lm( y.qup ~ 0 + x.qup ), col="darkgray", lty = 2 )
lm.qup <- lm( y.qup ~ 0 + x.qup )
summary( lm.qup )
text( x = 110, y = 200,
      labels = expression( atop( ''*italic( R^2 )*' = 0.9998' ) ) )
text( x = 110, y = 195,
      labels = expression( atop( ''*italic( w )*' = '*italic(  1.040474 * t )*'' ) ) )

# Stop graphics !
dev.copy( which = curr_plot )
dev.off()
dev.off()

# Check deviations
#  Error in the calibration: (m_new - m_post) / m_post
#   x.mean (m_post) <- mean_est_divt.post72sp[ind.post.72sp]*100
#   y.mean (m_new)  <- STnSB$mean_divt[ind.post.cal]*100
#
#  Error in the calibration: (qup_new - qup_post) / qup_post
#   x.qup (m_post)  <- quant_est_divt.post72sp[ind.post.72sp,2]*100
#   y.qup (m_new    <- STnSB$quant_divt[ind.post.cal,2]*100
#
#  Error in the calibration: (qlow_new - qlow_post) / qlow_post
#   x.qlow (m_post) <- quant_est_divt.post72sp[ind.post.72sp,1]*100
#   y.qlow (m_new)  <- STnSB$quant_divt[ind.post.cal,1]*100
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

write.table( x = sum.stats.ord, file = paste( "01_SBnST_", subtree, "_sumstats.tsv", sep = "" ),
             sep = "\t", quote = FALSE )


## PLOT TWO ST DENSITIES ##
## Density 1: M.auratus- C. griseus
curve( dst( x, xi = 0.091, omega = 0.024,  # manually adjusted numbers
            alpha = 3.756, nu = 209.065),  # same numbers as in original ST
       from=0, to=0.3, n = 5e2, col = "red", 
       main = "Manually adjusted (red) VS unadjusted (black) - M.aur-C.gri", adj = 1 )
curve( dst( x, xi = 0.088, omega = 0.022,  
            alpha = 3.756, nu = 209.065),  # same numbers as in original ST
       from=0, to=0.3, n = 5e2, add = TRUE, col = "black", adj = 1 )
# Add neighbouring node 
## P. maniculatus - M. ochrogaster
curve( dst( x, xi = 0.128, omega = 0.029,  
            alpha = 3.532, nu = 404.416),  
       from=0, to=0.3, n = 5e2, add = TRUE, col = "green", adj = 1 )
## Muridae
curve( dst( x, xi = 0.174, omega = 0.035,  
            alpha = 3.117, nu = 381.970),  
       from=0, to=0.3, n = 5e2, add = TRUE, col = "purple", adj = 1 )


## Density 2: Murinae
curve( dst( x, xi = 0.089, omega = 0.024,  # manually adjusted numbers
            alpha = 3.794, nu = 547.445),  # same numbers as in original ST
       from=0, to=0.5, n = 5e2, col = "black",
       main = "Manually adjusted (red) VS unadjusted (black) - Murinae", adj = 1 )
curve( dst( x, xi = 0.104, omega = 0.0361,  # manually adjusted numbers
            alpha = 3.794, nu = 547.445),   # same numbers as in original ST
       from=0, to=0.5, n = 5e2, col = "red", add = TRUE, adj = 1 )
# Add neighbouring node 
## Nannospalax galili - Muridae
curve( dst( x, xi = 0.382, omega = 0.022,  
            alpha = -0.670, nu = 36.812),  
       from=0, to=0.5, n = 5e2, add = TRUE, col = "green", adj = 1 )
## Muridae
curve( dst( x, xi = 0.174, omega = 0.035,  
            alpha = 3.117, nu = 381.970),  
       from=0, to=0.5, n = 5e2, add = TRUE, col = "purple", adj = 1 )