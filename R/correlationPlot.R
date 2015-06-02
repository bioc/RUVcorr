#' Correlation plot to compare estimated correlations with true correlations.
#'
#' \code{correlationPlot} produces a correlation plot to compare true and estimated
#   correlations.
#'
#' @param true A matrix of true gene-gene correlation values.
#' @param est A matrix of estimated gene expression values.
#' @param plot.genes A vector of indices of genes used in plotting; 
#' the suggested length of this vector is 18.
#' @param boxes A logical scalar to indicate whether boxes 
#' are drawn around sets of 6 genes; only available if \code{plot.genes} has length 18.
#' @param title A character string describing the title of the plot.
#' @return \code{correlationPlot} returns a plot.
#'
#' @details
#' The upper triangle of the correlation plot shows the true gene-gene correlation values,
#' while the lower triangle of the correlation plot shows the gene-gene correlation values
#' calculated from the estimated gene expression values. This is possible because correlation
#' matrices are symmetric.
#' @inheritParams graphics::mtext
#' @examples
#' Y<-simulateGEdata(500, 500, 10, 2, 5, g=NULL, Sigma.eps=0.1, 
#' 250, 100, intercept=FALSE, check.input=FALSE)
#' correlationPlot(Y$Sigma, Y$Y, title="Raw", 
#' plot.genes=c(sample(1:100, 6), sample(101:250, 6), sample(251:500, 6)))
#' @seealso \code{\link[corrplot]{corrplot}}
#' @author Saskia Freytag
#' @export
correlationPlot<-function(
      true, ##matrix 
      est,
      plot.genes=sample(1:dim(true)[1], 18), ## index of which genes to plot
      boxes=TRUE, ## draw boxes FALSE/TRUE
      title, 
      line=-1)
{
  
  if(is.matrix(true)==FALSE) stop("true needs to be a matrix object.")
  if(is.matrix(est)==FALSE) stop("est needs to be a matrix object.")
  
  if(length(plot.genes)!=18 & boxes==TRUE){
    warning("No boxes can be drawn when not exactly 18 genes are specified.")
    boxes=FALSE
  }
  
  est.true<-mashUp(true, est, plot.genes)
  ## lower triangle estimated values, upper triangle true values
  
  corrplot(est.true, method = c("color"), order = "original", tl.col = "black")
  ## draw object
  mtext(paste(title), side=3, line=line, cex = 1)
  if(boxes){
    rect(0.5, 12.5, 6.5, 18.5, border = "black", lwd = 2)
    rect(6.5, 6.5, 12.5, 12.5, border="black", lwd=2)
    rect(12.5, 0.5, 18.5, 6.5, border="black", lwd=2)
  }
}

