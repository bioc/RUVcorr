#' Quality assessment for cleaning procedures.
#'
#' \code{assessQuality} allows to assess the quality of cleaning procedures in the context
#'   of correlations when the true underlying correlation structure is known. 
#' 
#' @param est A matrix of estimated gene expression values.
#' @param true A matrix of true correlations.
#' @param index A vector of indices of genes to be included in 
#' the assessment; if \code{index="all"} all genes are considered.
#' @param methods The method used for quality assessment; 
#' if \code{method="fnorm"} the squared Frobenius norm is used;
#' if \code{method="wrong.sign"} the percentage of wrongly 
#' estimated signs is calculated if \code{method="all"}
#' both are calculated.
#' @return \code{assessQuality} returns a vector of the requested quality assessments.
#'
#' @details
#' The squared Frobenius norm used for \code{assessQuality} has the following structure
#'   \deqn{F=\frac{\|E-T\|^2}{s}}
#' Here, the parameter \eqn{E} and the parameter \eqn{T} denote the lower triangles of the 
#' estimated and true Fisher transformed correlation matrices, respectively.
#' The parameter \eqn{s} denotes the number of elements in \eqn{E} and \eqn{T}.
#' @examples
#' Y<-simulateGEdata(500, 500, 10, 2, 5, g=NULL, 
#' Sigma.eps=0.1, 250, 100, intercept=FALSE, check.input=FALSE)  
#' assessQuality(Y$Y, Y$Sigma, index=1:100, methods="wrong.sign")
#' assessQuality(Y$Y, Y$Sigma, index=1:100, method="fnorm")
#' @author Saskia Freytag
#' @export
assessQuality<-function(
      est, # estimated expression values in matrix
      true, # true correlation values in matrix
     index="all", # index of genes wanted for analysis
     methods=c("all", "fnorm", "wrong.sign") # methods
){
  if (index[1]!="all"){
    est<-est[, index]
    true<-true[index, index]
  }
  
  est<-cor(est)
  res<-vector(length=0)
  ## initialise result vector
  
  if(any(methods=="fnorm"| methods=="all")){
    true1<-fisherz(true[lower.tri(true)])
    est1<-fisherz(est[lower.tri(est)])
    res<-sum((true1-est1)^2)/length(true1)
    ## calculation of squared Frobenius Norm based on lower triangle of matrices
  }
  
  if(any(methods=="wrong.sign"| methods=="all")){
    neg<-which(true<0, arr.ind=TRUE)
    ## index of truly negative values
    pos<-which(true>0, arr.ind=TRUE)
    ## index of truly positive values
    res<-c(res, (sum(est[pos]<0) + sum(est[neg]>0))/length(true))
    ## calculation of percentage of correlations with wrong sign
  }
  if(length(res)==1){
    if(methods=="fnorm") {
      names(res)<-"fnorm"
    } else {names(res)<-"wrong.sign"}
  } else names(res)<-c("fnorm", "wrong.sign")
  
  return(res)
}
