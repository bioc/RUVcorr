#' Plot empirical cumulative distribution function for correlations.
#'
#' \code{ECDFPlot} generates empirical cumulative distribution 
#' functions (ECDF) for gene-gene correlation values.
#'
#' @param X A matrix or list of matrices of estimated gene-gene correlations.
#' @param Y A matrix of reference gene-gene correlations (i.e. 
#' underlying known correlation structure).
#' @param index A vector of indicies of genes of interest.
#' @param col.X The color or colors for ECDF as estimated from \code{X}.
#' @param col.Y The color for ECDF as estimated from \code{Y}.
#' @param title A character string describing title of plot.
#' @param legend A vector describing \code{X} and \code{Y}.
#' @return \code{ECDFPlot} returns a plot.
#'
#' @examples
#' Y<-simulateGEdata(500, 500, 10, 2, 5, g=NULL, Sigma.eps=0.1, 
#' 250, 100, intercept=FALSE, check.input=FALSE)
#' Y.hat<-RUVNaiveRidge(Y, center=TRUE, nc_index=251:500, 0, 10, check.input=TRUE)
#' Y.hat.cor<-cor(Y.hat)
#' par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0, mfrow=c(1, 1))
#' ECDFPlot(Y.hat.cor, Y$Sigma, index=1:100, title="Simulated data", 
#' legend=c("RUV", "Truth"))
#' ECDFPlot(list(Y.hat.cor, cor(Y$Y)), Y$Sigma, index=1:100, 
#' title="Simulated data", legend=c("RUV", "Raw", "Truth"), col.Y="black")
#' @author Saskia Freytag
#' @export
ECDFPlot<-function(
      X, ## matrix estimated correlations or list of matrices
      Y, ## matrix reference correlations
      index="all", ##vector
      col.X="red",
      col.Y="black",
      title,
      legend
){ 
  
  if(class(X)=="matrix"){
    if(index[1]!="all"){
      X<-X[index, index]
      Y<-Y[index, index]
    }
    
    X<-X[lower.tri(X)]
    Y<-Y[lower.tri(Y)]
    
    ECDF.X<-ecdf(abs(X))
    ECDF.Y<-ecdf(abs(Y))
    
    plot(ECDF.Y, xlim=c(0, 1), ylim=c(0, 1), col=col.Y, 
    xlab="|Correlation|", ylab="Proportion of |Correlation|", lwd=2,
    main=paste(title), verticals = TRUE, do.p=FALSE, bty="l")
    lines(ECDF.X, col = col.X, lwd=2, verticals=TRUE, do.p=FALSE)
    legend("bottomright", legend=paste(legend), col = c(col.X, col.Y), 
    bty="n", lty=1, lwd=3, cex=0.95)
  }
  
  if(class(X)=="list"){
    if(index[1]!="all"){
      X<-lapply(X, function(x) x[index, index])
      Y<-Y[index, index]
    }
    
    X<-lapply(X, function(x) x[lower.tri(x)])
    Y<-Y[lower.tri(Y)]
    
    ECDF.X<-lapply(X, function(x) ecdf(abs(x)))	
    ECDF.Y<-ecdf(abs(Y))
    
    if(length(col.X)!=length(X)){
      warning("Specified colors are no longer valid.")
      
      col.X<-hcl(h = seq(0, 360, round(360/length(X), 2)), c=45, l=70)[1:length(X)]
    }
    plot(ECDF.Y, xlim=c(0, 1), ylim=c(0, 1), col=col.Y, 
    xlab="|Correlation|", ylab="Proportion of |Correlation|", lwd=2,
    main=paste(title), verticals=TRUE, do.p=FALSE, bty="l")
    for (i in 1:length(X)){
      lines(ECDF.X[[i]], col=col.X[i], lwd=2, verticals=TRUE, do.p=FALSE)
    }
    legend("bottomright", legend=paste(legend), col=c(col.X, col.Y), 
    bty="n", lty=1, lwd=3, cex=0.95)
    
  }
  
}