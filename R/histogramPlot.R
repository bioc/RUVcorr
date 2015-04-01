#' Plot histogram of correlations.
#'
#' \code{histogramPlot} plots histograms of correlation values in expression data and
#'	its reference.
#'
#' @param X A matrix or a list of matrices of estimated gene-gene correlations.
#' @param Y A matrix of reference gene-gene correlations (i.e. known underlying correlation structure).
#' @param legend A vector of character strings describing the data contained in \code{X} and \code{Y}.
#' @param title A character string describing title.
#' @param col.X A vector or character string defining the color/colors associated with the data contained in \code{X}.
#' @param col.Y The color associated with the data in \code{Y}.
#' @param line A vector giving the line type.
#' @return \code{histogramPlot} returns a plot.
#'
#' @inheritParams graphics::hist
#' @details
#' The default for breaks is \code{"Sturges"}.
#' Other names for which algorithms are supplied are \code{"Scott"} and \code{"FD"} / \code{"Freedman-Diaconis"} 
#' Case is ignored and partial 
#' matching is used. Alternatively, a function can be supplied which will compute the 
#' intended number of breaks or the actual breakpoints as a function of \code{x}.
#' @examples
#' Y<-simulateGEdata(500, 500, 10, 2, 5, g=NULL, Sigma.eps=0.1, 
#' 250, 100, intercept=FALSE, check.input=FALSE)
#' Y.hat<-RUVNaiveRidge(Y, center=TRUE, nc_index=251:500, 0, 10, check.input=FALSE)
#' Y.hat.cor<-cor(Y.hat[,1:100])
#' try(dev.off(), silent=TRUE)
#' par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0, mfrow=c(1, 1))
#' histogramPlot(Y.hat.cor, Y$Sigma[1:100, 1:100], title="Simulated data", 
#' legend=c("RUV", "Truth"))
#' try(dev.off(), silent=TRUE)
#' histogramPlot(list(Y.hat.cor, cor(Y$Y[, 1:100])), Y$Sigma[1:100, 1:100],
#' title="Simulated data", col.Y="black", legend=c("RUV", "Raw", "Truth"))
#' @author Saskia Freytag
#' @export
histogramPlot<-function(
      X, ##Matrix or list of matrices of estimated correlations.
      Y, ##Matrix of reference correlations.
      legend, ##Vector of characters describing input matrices
      breaks=40, 
      title,
      col.X="red",
      col.Y="black",
      line=NULL
){
  
  if(class(X)=="matrix"){
    X<-X[lower.tri(X)]
    Y<-Y[lower.tri(Y)]
    max.val<-max(c(max(density(X)$y),max(density(Y)$y)))
    hist(Y, freq=FALSE, xlim=c(-1, 1), breaks=breaks, ylim=c(0,max.val),
    main=paste(title), xlab="Correlation Size", border=col.Y)
    if(is.null(line)) line<-1
    lines(density(X), col=col.X, lwd=2, lty=line) 
    legend("topleft", paste(legend), bty="n" ,lty=1, lwd=3, cex=0.95, ncol=1, col=c(col.X,col.Y))
  }
  
  if(class(X)=="list"){
    if(length(col.X)!=length(X)){
      warning("Specified colors are no longer valid.")
      
      col.X<-hcl(h = seq(0,360,round(360/length(X),2)), c=45, l=70)[1:length(X)]
    }
    X<-lapply(X, function(x) x[lower.tri(x)])
    Y<-Y[lower.tri(Y)]
    max.val<-max(c(max(unlist(lapply(X, function(x) max(density(x)$y)))), max(density(Y)$y)))
    hist(Y, freq=FALSE, xlim=c(-1, 1), breaks=breaks, ylim=c(0, max.val),
    main=paste(title), xlab="Correlation Size", border=col.Y)
    if(is.null(line)) line<-rep(1, length(X))
    
    for(i in 1:length(X)){
      lines(density(X[[i]]), col=col.X[i], lwd=2, lty=line[i])
    }
    
    legend("topleft", paste(legend), bty="n" , lty=line, lwd=3, cex=0.95, ncol=1, col=c(col.X, col.Y))
  }
  
}
