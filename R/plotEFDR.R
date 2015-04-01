#' Plots an object of class \code{EFDR}.
#'
#' \code{plotEFDR} plots the EDFR curves stored in the objects of class \code{EDFR}.
#'
#' @param x An object of class \code{EFDR} or a list of objects of class \code{EFDR}.
#' @param main A character string describing the title of the plot.
#' @param legend A vector of character strings decribing the different EDFR objects in
#' \code{x}; only applicable when x is a list.
#' @param col A vector giving the colors, if \code{NULL} colors
#' are generated automatically.
#' @param ... Further arguments passed to or from other methods.
#' @return \code{plot.EDFR} returns a plot.
#' 
#' @examples
#' Y<-simulateGEdata(500, 500, 10, 2, 5, g=NULL, Sigma.eps=0.1, 
#' 250, 100, intercept=FALSE, check.input=FALSE)
#' anno<-as.matrix(sample(1:4, dim(Y$Y)[1], replace=TRUE))
#' colnames(anno)<-"Factor"
#' weights<-findWeights(Y$Y, anno, "Factor")
#' efdr<-calculateEFDR(Y$Y, exclude=1:100, index.ref=1:10, 
#' Weights=weights, anno=anno, Factor="Factor")
#' plotEFDR(efdr)
#' @author Saskia Freytag
#' @export
plotEFDR<-function(x,
      main="", 
      legend,
      col=NULL, 
      ...)
{
  
  if(class(x)!="EFDR" & class(x)!="list"){
    stop("x must be of class EFDR or a list of EFDR objects.")
  }
  if(class(x)=="EFDR"){
    plot(x$Thresholds, x$EFDR.values, type="b", lwd=2, xlim=c(0, 1), 
    ylim=c(0, 1), col="red", ylab="Empirical FDR", xlab="Threshold", bty="l", ...)
  }
  if(class(x)=="list"){
    if(any(sapply(x, function(i) class(i)=="EFDR"))==FALSE) {
      stop("x must be of class EFDR or a list of EFDR objects.")
    }
    plot(x[[1]]$Thresholds, x[[1]]$EFDR.values, type="n", lwd=2, xlim=c(0, 1), 
    ylim=c(0, 1), ylab="Empirical FDR", xlab="Threshold", bty="l", main=main, ...)
    if(is.null(col)){
      col<-hcl(h = seq(0, 360, round(360/length(x), 2)), c=45, l=70)
    }
    lapply(1:length(x), function(i) points(x[[i]]$Thresholds, x[[i]]$EFDR.values, 
    type="b", col=col[i], lwd=2))
     
     legend("bottomleft", legend, col=col[1:length(x)], bty="n", lwd=2)
  }
  	
}
