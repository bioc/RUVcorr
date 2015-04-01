#' Randomly choose background genes.
#'
#' \code{background} returns background genes for judging the quality of the cleaning. 
#'   These genes are supposed to represent the majority of genes. The positive control and 
#'   negative control genes should be excluded.
#' 
#' @param Y A matrix of gene expression values or an object of 
#' the class \code{simulateGEdata}.
#' @param nBG An integer setting the number of background genes.
#' @param exclude A vector of indices of genes to exclude.
#' @param nc_index A vector of indices of negative controls (also
#'  excluded from being background genes).
#' @return \code{background} returns a vector of randomly chosen indices.
#'
#' @examples
#' Y<-simulateGEdata(500, 500, 10, 2, 5, g=NULL, Sigma.eps=0.1, 
#' 250, 100, intercept=FALSE, check.input=FALSE)  
#' background(Y, nBG=20, exclude=1:100, nc_index=251:500)
#' @author Saskia Freytag
#' @export  
background<-function(Y,
      nBG, ##number of background genes
      exclude,
      nc_index){
  
  if(is.simulateGEdata(Y)){
    Y<-Y$Y
    colnames(Y)<-1:dim(Y)[2]	
  }	
  Ytmp<-Y[,-c(exclude,nc_index)]
  
  if(dim(Ytmp)[2]<=nBG) stop("Trying to pick too many backgrounds.")
  BG<-sample(colnames(Ytmp), nBG)
  
  return(which(is.element(colnames(Y),BG)))
}