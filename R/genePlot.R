#' Plot of means and inter-quantile ranges of all genes.
#' 
#' \code{genePlot} plots the means vs. the inter-quantile ranges of the gene
#' 	expression values of all genes with the possibility to highlight interesting sets of genes.
#'
#' @param Y A matrix of gene expression values or an object of the class \code{simualteGEdata}.
#' @param index A vector of indices of genes of interest to be displayed in a different color, if \code{index=NULL} no genes are highlighted.
#' @param legend A character string describing the highlighted genes.
#' @param col.h The color of the highlighted genes.
#' @param title A character string describing the title of the plot.
#' @return \code{genePlot} returns a plot.
#'
#' @examples
#' Y<-simulateGEdata(500, 500, 10, 2, 5, g=NULL, Sigma.eps=0.1, 
#' 250, 100, intercept=FALSE, check.input=TRUE)
#' try(dev.off(), silent=TRUE)
#' par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
#' genePlot(Y, index=1:100, legend="Expressed genes", title="IQR-Mean Plot")
#' @author Saskia Freytag
#' @export
genePlot<-function(
      Y, ## simulateGEdata object or matrix 
      index=NULL, ## index of genes that should be highlighted
      legend=NULL, # description
      col.h="red",
      title
){
  
  if(is.simulateGEdata(Y)) Y<-Y$Y
  
  tmpMeans<-colMeans(Y)
  tmpIQR<-apply(Y, 2, function(x) findIQR(x))
  
  plot(tmpMeans, tmpIQR, col="gray", pch=20, xlab="Means", ylab="IQR", 
  main=paste(title), bty="l")
  if(is.null(index)==FALSE){
    points(tmpMeans[index], tmpIQR[index], col=col.h, pch=20)
  }
  if(is.null(legend)==FALSE){
    legend("topright", legend=legend, text.col=c(col.h), bty="n")
  }		
}