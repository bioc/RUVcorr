#' Calculate weighted correlations.
#' 
#' \code{wcor} returns correlations weighted according to a provided object of
#' class \code{Weights}.
#'
#' @param X A matrix of gene expression values.
#' @param anno A dataframe or a matrix containing the annotation of arrays in \code{X}.
#' @param Factor A character string corresponding to a column name of \code{anno}; this
#' should be the same used to generate \code{Weights}.
#' @param Weights An object of class \code{Weights} or a list of weights.
#' @return \code{wcor} returns a matrix.
#'
#' @examples
#' Y<-simulateGEdata(500, 500, 10, 2, 5, g=NULL, Sigma.eps=0.1, 
#' 250, 100, intercept=FALSE, check.input=FALSE)
#' anno<-as.matrix(sample(1:5, dim(Y$Y)[1], replace=TRUE))
#' colnames(anno)<-"Factor"
#' weights<-findWeights(Y$Y, anno, "Factor")
#' wcor(Y$Y[,1:5], anno, "Factor", weights)
#' @author Saskia Freytag
#' @export
wcor<-function(
      X, 
      anno, 
      Factor, 
      Weights)
{
  
  if(class(Weights)=="Weights") Weights<-Weights$Weights
  if(any(colnames(anno)==Factor)==FALSE){ 
    stop("Factor has to be a column of anno.")
  }
  if(dim(anno)[1]!=dim(X)[1]) {
    stop("Annotation does not match.")
  }
  
  col.num<-which(colnames(anno)==Factor)
  categories<-sort(unique(anno[, col.num]))
  if(length(categories)!=length(Weights)) {
     stop("Need to have weights for all levels of Factor!")
  }
  
  expression<-splitByFactor(X, anno, Factor)
  all.cor<-lapply(1:length(expression), function(x) 
  cor(expression[[x]])*as.numeric(Weights[x]))
  return(Reduce("+", all.cor))
  
}
