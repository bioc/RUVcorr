#' Plots an object of class \code{Threshold}.
#'
#' \code{plotThreshold} plots the objects of class \code{Threshold}.
#'
#' @param x An object of class \code{Threshold} or a list of objects of class \code{Threshold}.
#' @param main A character string describing the title of the plot.
#' @param legend A vector of character strings decribing the different \code{Threshold} objects in
#' \code{x}; only applicable when x is a list.
#' @param col A vector giving the colors, if \code{NULL} colors
#' are generated automatically.
#' @param ... Further arguments passed to or from other methods.
#' @return \code{plotThreshold} returns a plot.
#' 
#' @examples
#' Y<-simulateGEdata(500, 500, 10, 2, 5, g=NULL, Sigma.eps=0.1, 
#' 250, 100, intercept=FALSE, check.input=FALSE)
#' anno<-as.matrix(sample(1:4, dim(Y$Y)[1], replace=TRUE))
#' colnames(anno)<-"Factor"
#' weights<-findWeights(Y$Y, anno, "Factor")
#' Thresh<-calculateThreshold(Y$Y, exclude=1:100, index.ref=1:10, 
#' Weights=weights, anno=anno, Factor="Factor")
#' plotThreshold(Thresh)
#' @author Saskia Freytag
#' @seealso \code{\link{calculateThreshold}}
#' @export
plotThreshold<-function(x,
      main="", 
      legend,
      col=NULL, 
      ...)
{
  
  if(class(x)!="Threshold" & class(x)!="list"){
    stop("x must be of class Threshold or a list of Threshold objects.")
  }
  if(class(x)=="Threshold"){
    plot(x$Thresholds, x$Prop.values, type="b", lwd=2, xlim=c(0, 1), 
    ylim=c(0, 1), col="red", ylab="Proportion Prioritised", xlab="Threshold", bty="l", ...)
  }
  if(class(x)=="list"){
    if(any(sapply(x, function(i) class(i)=="Threshold"))==FALSE) {
      stop("x must be of class Threshold or a list of Threshold objects.")
    }
    plot(x[[1]]$Thresholds, x[[1]]$Prop.values, type="n", lwd=2, xlim=c(0, 1), 
    ylim=c(0, 1), ylab="Proportion Priortised", xlab="Threshold", bty="l", main=main, ...)
    if(is.null(col)){
      col<-hcl(h = seq(0, 360, round(360/length(x), 2)), c=45, l=70)
    }
    lapply(1:length(x), function(i) points(x[[i]]$Thresholds, x[[i]]$Prop.values, 
    type="b", col=col[i], lwd=2))
     
     legend("bottomleft", legend, col=col[1:length(x)], bty="n", lwd=2)
  }
  	
}
