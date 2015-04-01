#' Plot eigenvalues of SVD of the negtaive controls.
#'
#' \code{eigenvaluePlot} plots the ratio of the i$^{th}$ eigenvalue 
#' of the SVD of the negative controls to the eigenvalue total.
#'
#' @param Y A matrix of gene expressions.
#' @param nc_index A vector of indices for the negative controls.
#' @param k A numeric value giving the number of eigenvalues that should be displayed.
#' @param center A logical character to indicate whether centering is needed.
#' @param title A character string describing title.
#' @return \code{eigenvaluePlot} returns a plot.
#'
#' @examples
#' Y<-simulateGEdata(500, 500, 10, 2, 5, g=NULL, Sigma.eps=0.1, 
#' 250, 100, intercept=FALSE, check.input=FALSE)
#' eigenvaluePlot(Y$Y, nc_index=251:500, k=20, center=TRUE)
#' @author Saskia Freytag
#' @export

eigenvaluePlot<-function(Y, nc_index, k=10, center=TRUE, title="Eigenvalue Plot"){
  Yc<-Y[, nc_index]
  Yc<-scale(Yc, center=TRUE, scale=FALSE)
  
  out<-svd(Yc)$d
  out<-out/sum(out)
  
  plot(out[1:k], xlab="Index", ylab="Proportion of Eigenvalue", 
  main=paste(title), bty="l")
  
}