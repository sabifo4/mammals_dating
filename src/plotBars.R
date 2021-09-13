# Helper function to reorder tips 
# Code taken from StackOverflow to reoder the tips 
# so they are ladderized properly
# https://stackoverflow.com/questions/34364660/how-to-get-correct-order-of-tip-labels-in-ape-after-calling-ladderize-function
reorder_tips <- function( phy ){
  
  # First step is to filter out internal nodes from the the second
  # column of the edge matrix
  phy    <- ladderize(phy, right = TRUE)
  is_tip <- phy$edge[,2] <= length(phy$tip.label)
  ordered_tips <- phy$edge[is_tip, 2]
  
  #Then you can use this vector to extract the tips in the right order:
  phy$tip.label[ordered_tips]
  
  # Return ladderized phy object
  return( phy )
}

## Function 2
# New function that uses the plotBars function from phytools
# as a basis to now plot much cooler stuff !
# New args:
#
# tree: phylo, object of class "phylo" with the tree.
#   CI: list, length equal to the number of models that you
#       want to plot. Each element is an object of class "data.frame"
#       where the columns correspond to the "mean", "low.bounds", 
#       and "upper.bounds" calculated from the values sampled during
#       the MCMC. Bear in mind that the first element needs to correspond
#       to the data of the model you want to plot in the nodes (main model),
#       the second model to the one you want to use to get the nodes in which
#       you had the calibrations (prior), and then add as many data frames as 
#       other models you have tested.
#       
# models: integer, number that indicates how many models you have 
#         ## note: redundant, I will eventually remove it.
# transp: numeric, vector which length corresponds to the number of models 
#         tested. Each number is used as "alpha" to control the transparency
#         of the color in the error bars.
# cex.dots: numeric, number to control the size of the dots and bars
# reduce.y: numeric, number that controls how much "y" needs to decrease to plot 
#           the error bars of all the models.
# calibs: data.frame, it contains the values with the calibrations.
#         First column defins which calibration is used ("ST" for skew-T
#         and "B" for soft bounds). The rest of columns contain the parameters
#         for these two calibrations specified in MCMCtree.
# ladder.tt: boolean, TRUE if you want to ladderize your tree. FALSE, otherwise.
new_plotBars <- function (tree, CI, models, transp, cex.dots,
                          reduce.y, calibs, ladder.tt, 
                          text.calibs = NULL, name.nodes = NULL, 
                          reduce.y.calibs = NULL, 
                          cex.y.node = 1, cex.x.node = 0.05, 
                          cex.x.calib = 1, cex.y.calib = 1, 
                          plot.calibs = TRUE, ...) 
{
  ## SAC: Get only L and U columns of CI
  ## as well as the mean ages for first model 
  mean.ages.1 <- CI[[1]][,1]
  for( i in 1:models ){
    CI[[ i ]] <- CI[[ i ]][,c(2:3)]
  }
  
  args <- list(...)
  if (!is.null(args$gridlines)) {
    gridlines <- args$gridlines
    args$gridlines <- NULL
  }
  else gridlines <- TRUE
  if (is.null(args$mar)) 
    args$mar <- c(4.1, 1.1, 1.1, 1.1)
  if (is.null(args$ftype)) 
    args$ftype <- "i"
  ## SAC note: "ftype" is the parameter 
  ## that actualliy takes size for tip labels !
  fsize <- if (!is.null(args$fsize)) 
    args$fsize
  else 1
  if (is.null(args$direction)) 
    args$direction <- "leftwards"
  if (!is.null(args$bar.width)) {
    bar.width <- args$bar.width
    args$bar.width <- NULL
  }
  else bar.width <- 11
  if (!is.null(args$cex)) {
    cex <- args$cex
    args$cex <- NULL
  }
  else cex <- 0.5
  print( cex )
  if (!is.null(args$bar.col)) {
    bar.col <- args$bar.col
    args$bar.col <- NULL
  }
  else bar.col <- "blue"
  par(mar = args$mar)
  plot.new()
  th <- max(nodeHeights(tree))
  ## ------------------------------------------------- ##
  ## SAC: Here we will use the maximum CI of the first ##
  ## model, assuming this is the one we want to have   ##
  ## as the "main" model                               ##
  ## ------------------------------------------------- ##
  # h <- max(th, max(CI))         ## before SAC
  CI.for.lims <- CI[[ 1 ]]        ## after SAC
  h <- max(th, max(CI.for.lims))  ## after SAC
  if (is.null(args$xlim)) {
    # m <- min(min(nodeHeights(tree)), min(CI))         ## before SAC
    m <- min(min(nodeHeights(tree)), min(CI.for.lims))  ## after SAC
    d <- diff(c(m, h))
    pp <- par("pin")[1]
    sw <- fsize * (max(strwidth(tree$tip.label, units = "inches"))) + 
      1.37 * fsize * strwidth("W", units = "inches")
    alp <- optimize(function(a, d, sw, pp) (a * 1.04 * d + 
                                              sw - pp)^2, d = d, sw = sw, pp = pp, interval = c(0, 
                                                                                                1e+06))$minimum
    args$xlim <- if (args$direction == "leftwards") 
      c(h, m - sw/alp)
    else c(m, h + sw/alp)
  }
  if (is.null(args$at)){ 
    at <- seq(0, h, by = h/5)
    ## SAC: Add K-Pg point
    if( ! 0.66 %in% at ){
      at <- sort( c( at, 0.66 ) )
    }
    ## SAC: Add root age point
    if( ! mean.ages.1[1] %in% at ){
      at <- sort( c( at, mean.ages.1[1] ) )
      print( at )
    }
  }
  else {
    at <- args$at
    args$at <- NULL
  }
  if( ladder.tt == TRUE ){
    args$tree <- ladderize( tree )
  }else{
    args$tree <- tree
  }
  args$add <- TRUE
  print(args)
  do.call(plotTree, args = args)
  if (gridlines) 
    abline(v = at, lty = "dashed", col = make.transparent("grey", 
                                                          0.5))
    ind.KPg <- which( at == 0.66 )
    abline(v = at[ind.KPg], lty = "dashed", col = make.transparent("black", 
                                                        0.8))
  axis(1, at = at, labels = signif(at, 3),
       cex.axis = cex, las = 2 )
  obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)

  ## SAC: Run this for every model passed to the function
  ##      Note that the first model is the main one, then 
  ##      the second model is the one with the calibrations
  ##      and then other models can be added
  for( mod in 1:models ){ 
    
    ## SAC: Set vars in for loop
    CI.mod       <- CI[[ mod ]]
    bar.col.mod  <- bar.col[ mod ]
    alpha.transp <- transp[ mod ]
    count        <- 0
    for (i in 1:tree$Nnode + Ntip(tree)){
      
      ## SAC: mod == 2 is the one with calibrations
      if( mod == 2 && plot.calibs == TRUE ){
        ind.calib       <- which( rownames(CI.mod) %in% rownames( calibs ) )
        tmp.calib.nodes <- as.numeric( gsub( "t_n", "", rownames(CI.mod)[ind.calib] ) )
        if( i %in% tmp.calib.nodes ){
          if( is.null( reduce.y.calibs ) == TRUE ){
            reduce.y.calibs <- reduce.y
          }
          lines(x = c(CI.mod[i - Ntip(tree), 1], CI.mod[i - Ntip(tree), 2]),
                y = rep(obj$yy[i] - reduce.y.calibs, 2),
                lwd = 0.8, lend = 0, # col = make.transparent(bar.col.mod, 0.4))
                col = bar.col.mod )
          ## -- SAC: my old code not needed to plot this now -- ##
          # points( x = (CI.mod[i - Ntip(tree), 1] + CI.mod[i - Ntip(tree), 2])/2,
          #         y = obj$yy[i] - reduce.y.calibs, pch = 19, col = bar.col.mod,
          #         cex = cex.dots )
          ## -------------------------------------------------- ##
          points( x = c(CI.mod[i - Ntip(tree), 1], CI.mod[i - Ntip(tree), 2]),
                  y = rep(obj$yy[i] - reduce.y.calibs, 2),
                  pch = "|", col = bar.col.mod, cex = cex.dots )
        }
        
      } ## SAC: mod == 1 is the main model 
      if( mod == 1 ){
        lines(x = c(CI.mod[i - Ntip(tree), 1], CI.mod[i - Ntip(tree), 2]),
              y = rep(obj$yy[i], 2),
              lwd = bar.width, #lend = 0, # col = make.transparent(bar.col.mod, 0.4))
              col = make.transparent(bar.col.mod, alpha.transp))
        ## SAC: Print node labels if calibration is there
        ind.calib       <- which( rownames( CI.mod ) %in% rownames( calibs ) )
        tmp.calib.nodes <- as.numeric( gsub( "t_n", "", rownames(CI.mod)[ind.calib] ) )
        if( i %in% tmp.calib.nodes ){
          count <- count + 1
          text( obj$xx[i]+cex.x.node,
                obj$yy[i]+cex.y.node, 
                labels = paste( i ), cex = 0.8 )
          if( is.null( text.calibs ) == FALSE ){
            text( obj$xx[i]+cex.x.calib,
                  obj$yy[i]-cex.y.calib, 
                  labels = text.calibs[count], cex = 0.4 )
          }
          if( is.null( name.nodes ) == FALSE ){
            text( obj$xx[i]+cex.x.node,
                  obj$yy[i]-cex.y.node, 
                  labels = names( name.nodes )[count], cex = 0.4 )
          }
          
          
        }
      } ## SAC: mod > 2 --> the rest of models ! 
      else if( mod > 2 ){
        lines(x = c(CI.mod[i - Ntip(tree), 1], CI.mod[i - Ntip(tree), 2]),
              y = rep(obj$yy[i] - reduce.y, 2),
              lwd = bar.width, lend = 0, # col = make.transparent(bar.col.mod, 0.4))
              col = make.transparent(bar.col.mod, alpha.transp))
        points( x = (CI.mod[i - Ntip(tree), 1] + CI.mod[i - Ntip(tree), 2])/2,
                y = obj$yy[i] - reduce.y, pch = 19, col = bar.col.mod,
                cex = cex.dots )
      }else if( mod == 2 && plot.calibs == FALSE  ){
        lines(x = c(CI.mod[i - Ntip(tree), 1], CI.mod[i - Ntip(tree), 2]),
              y = rep(obj$yy[i] - reduce.y, 2),
              lwd = bar.width, lend = 0, # col = make.transparent(bar.col.mod, 0.4))
              col = make.transparent(bar.col.mod, alpha.transp))
        points( x = (CI.mod[i - Ntip(tree), 1] + CI.mod[i - Ntip(tree), 2])/2,
                y = obj$yy[i] - reduce.y, pch = 19, col = bar.col.mod,
                cex = cex.dots )
      }
      
    }
    ## SAC: Add points for model == 1
    if( mod == 1 ){
      points(obj$xx[1:tree$Nnode + Ntip(tree)],
             obj$yy[1:tree$Nnode + Ntip(tree)], pch = 19, col = bar.col.mod,
             cex = cex.dots )
    }
    
    ## SAC: Update reduce.y y axis (decrease) for next model
    reduce.y <- reduce.y + reduce.y 
  }
  
  ## SAC: Return object
  return( obj )
  
}
