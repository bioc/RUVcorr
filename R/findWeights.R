#' Finds weights of each level of a factor.
#' 
#' \code{findWeights} returns a list of variances and weights based on the correlation
#' between genes for each level of a factor found in the annotation. This function is 
#' typically used to find the weights of each individual in the data set.
#'
#' @param X A matrix of gene expression values.
#' @param anno A dataframe or a matrix containing the annotation of arrays in \code{X}.
#' @param Factor A character string corresponding to a column name of \code{anno}. For all
#' levels of this factor corresponding weights will be calculated.
#' @return \code{findWeights} returns output of the class \code{Weights}.
#' An object of class \code{Weights} is a \code{list} with the following components:
#'     \itemize{
#'        \item{\code{Weights}}{ A list containing the weights of each level of \code{Factor}.}
#'        \item{\code{Inv.Sigma}}{ A list containing the inverse variances of each level of \code{Factor}.}
#'     }
#'
#' @details
#' Note that because calculations of weights include finding correlations between all genes,
#' this function might take some time. Hence, recalculation of weights is not advisable and 
#' should be avoided. However often the inverse variances can be used to calculate new weights. 
#' In particlular, when \eqn{W_i} denotes the weight of the \eqn{i^{th}} level and \eqn{V_i}
#' the variance as calculated from the gene-gene correlations:
#' \deqn{W_i=\frac{\frac{1}{V_i}}{\sum_{i=1}^{n}\frac{1}{V_i}}}
#'
#' @examples
#' Y<-simulateGEdata(500, 500, 10, 2, 5, g=NULL, Sigma.eps=0.1, 
#' 250, 100, intercept=FALSE, check.input=FALSE)
#' anno<-as.matrix(sample(1:4, dim(Y$Y)[1], replace=TRUE))
#' colnames(anno)<-"Factor"
#' findWeights(Y$Y, anno, "Factor")
#' @author Saskia Freytag
#' @export
findWeights<-function(X, 
      anno, 
      Factor
){
  
  if(any(colnames(anno)==Factor)==FALSE) stop("Factor has to be a column of anno.")
  if(dim(anno)[1]!=dim(X)[1]) stop("Annotation does not match.")
  
  expression<-splitByFactor(X, anno, Factor)
  expression.cor<-lapply(expression, function(x) cor(x)[lower.tri(cor(x))])
  
  Inv.Sigma<-lapply(expression.cor, function(x) 1/var(x))
  Weights<-lapply(Inv.Sigma, function(x) x/Reduce("+", Inv.Sigma))
  res<-list(Weights=Weights, Inv.Sigma=Inv.Sigma)
  class(res)<-"Weights"
  
  return(res)
}
