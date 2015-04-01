#' Prioritising candidate genes.
#' 
#' \code{prioritise} returns a set of genes from a candidate set of genes that are
#' correlated above a provided threshold with at least one of the provided reference
#' genes.
#' 
#' @param X A matrix of gene expression values.
#' @param ref_index A vector of indices of reference genes.
#' @param cand_index A vector of indices of candidate genes. 
#' @param anno A dataframe or a matrix containing the annotation of arrays in \code{X}.
#' @param Factor A character string corresponding to a column name of \code{anno}; this
#' should be the same used to generate \code{Weights}.
#' @param Weights An object of class \code{Weights} or a list of weights. If \code{NULL}
#' the unweighted correlation is used.
#' @param threshold A value in the range \eqn{[0,1]}. 
#' @return \code{prioritise} returns a matrix with three columns. The first column gives 
#' the names of the genes that were prioiritised, while the second column gives the 
#' number of correlations above the threshold for the gene in question. The
#' columns gives the sum of the absolute value of all correlations with reference genes
#' above the threshold. 
#' @examples
#' Y<-simulateGEdata(500, 500, 10, 2, 5, g=NULL, Sigma.eps=0.1, 
#' 250, 100, intercept=FALSE, check.input=TRUE)
#' colnames(Y$Y)<-1:dim(Y$Y)[2]
#' anno<-as.matrix(sample(1:5, dim(Y$Y)[1], replace=TRUE))
#' colnames(anno)<-"Factor"
#' weights<-findWeights(Y$Y, anno, "Factor")
#' prioritise(Y$Y, 1:10, 51:150, anno, "Factor", weights, 0.6)
#' @author Saskia Freytag
#' @export
prioritise<-function(X, 
      ref_index, 
      cand_index, 
      anno, 
      Factor, 
      Weights, 
      threshold)
{ 
  if(is.null(Weights)!=TRUE){
    tmp<-wcor(X[, c(ref_index, cand_index)], anno, Factor, Weights)
  } else {
    tmp<-cor(X[, c(ref_index, cand_index)])
  }
  tmp<-tmp[1:length(ref_index), (length(ref_index)+1):dim(tmp)[2]]
  index<-apply(tmp, 2, function(z) any(abs(z)>=threshold))
  prioritisedGenes<-colnames(X)[cand_index[index]]
  strength<-apply(tmp, 2, function(z) length(which(abs(z)>=threshold)))[index]
  strength2<-apply(tmp, 2, function(z) sum(abs(z[which(abs(z)>=threshold)])))[index]
  
  return(cbind(prioritisedGenes, strength, strength2))
}