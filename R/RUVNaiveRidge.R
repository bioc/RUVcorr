#' Removal of unwanted variation for gene correlations.
#' 
#' \code{RUVNaiveRidge} applies the ridged version of global removal of unwanted variation 
#' to simulated or real gene expression data. 
#'
#' @param Y A matrix of gene expression values or an object of 
#' class \code{simulateGEdata}. 
#' @param center A logical scalar; if \code{TRUE} the data is centered, 
#' if \code{FALSE} data is assumed to be already centered.
#' @param nc_index A vector of indices of negative controls.
#' @param nu A numeric scalar value of \code{nu} \eqn{\geq 0}.
#' @param kW An integer setting the number of dimensions for the estimated noise.
#' @param check.input A logical scalar; if \code{TRUE} all input is 
#' checked (not advisable for large simulations).
#' @return \code{RUVNaiveRidge} returns a matrix of the cleaned 
#' (RUV-treated) centered gene expression values.
#' @details 
#' The parameter \code{kW} controls how much noise is cleaned, whereas the 
#' parameter \code{nu} controls the amount of ridging to deal with possible dependence of 
#' the noise and the factor of interest.
#' @examples 
#' Y<-simulateGEdata(500, 500, 10, 2, 5, g=NULL, Sigma.eps=0.1, 
#' 250, 100, intercept=TRUE, check.input=FALSE)
#' Y
#' Y.hat<-RUVNaiveRidge(Y, center=TRUE, nc_index=251:500, 0, 9, check.input=TRUE)
#' cor(Y.hat[,1:5])
#' Y$Sigma[1:5,1:5]
#' Y.hat<-RUVNaiveRidge(Y, center=FALSE, nc_index=251:500, 0, 10, check.input=TRUE)
#' cor(Y.hat[,1:5])
#' Y$Sigma[1:5,1:5]
#' @references Jacob L., Gagnon-Bartsch J., Speed T. Correcting gene expression 
#'   data when neither the unwanted variation nor the factor of interest are observed. 
#'   Berkley Technical Reports (2012).
#' @author Saskia Freytag, Laurent Jacob
#' @exportMethod RUVNaiveRidge
#' @export
RUVNaiveRidge<-function(
      Y, 
      center=TRUE, ##set equal to FALSE in case of centered data
      nc_index, ##column index for negative controls 
      nu, ## Ridge factor
      kW, ## number of noise dimensions
      check.input=FALSE) UseMethod("RUVNaiveRidge")

#' \code{RUVNaiveRidge.default} applies the ridged version of global 
#' removal of unwanted variation to matrices.
#'
#' @rdname RUVNaiveRidge
#' @export
RUVNaiveRidge.default<-function(
      Y, ##matrix of gene expression data 
      center=TRUE, ##set equal to FALSE in case of centered data
      nc_index, ##column index for negative controls 
      nu, ## Ridge factor
      kW, ## number of noise dimensions
      check.input=FALSE)
{
  
  if(check.input){
    if(is.matrix(Y)==FALSE){stop("Y needs to be a matrix.")}
    if(nu<0){stop("nu has to be positive or 0.")}
    if(kW>dim(Y)[1]){ stop("kW is too big.") }
   }
  
  if(center){
    Y<-scale(Y, center=TRUE, scale=FALSE)
    ## center data
  }
  
  Yc<-Y[, nc_index]
  ## subset negative controls
  
  tmp<-svd(Yc, nu=kW, nv=kW)
  S.d<-diag(tmp$d[1:kW], nrow=kW, ncol=kW)
  ## SVD of negative controls
  
  W.hat<-tmp$u%*%S.d
  ## estimate W.hat
  
  alpha.hat<-solve(t(W.hat)%*%W.hat+nu*diag(dim(W.hat)[2]))%*%t(W.hat)%*%Y
  ## estimate alpha.hat
  
  return(Y-W.hat%*%alpha.hat)
  ## calculate Y.hat (note that this is the mean centered Y)
}

#' \code{RUVNaiveRidge.simulateGEdata} applies the ridged version of 
#' removal of unwanted variation to objects of class \code{simulateGEdata}.
#'
#' @rdname RUVNaiveRidge 
#' @export
RUVNaiveRidge.simulateGEdata<-function(
      Y, ##object of the class simulateGEdata
      center=TRUE, ##set equal to FALSE in case of centered data
      nc_index, ##column index for negative controls 
      nu, ## Ridge factor
      kW, ## number of noise dimensions
      check.input=FALSE){
  
  Y<-Y$Y 
  ## extract right expression values
  
  RUVNaiveRidge.default(Y, center, nc_index, nu, kW, check.input)
}
