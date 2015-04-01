#' Calculates the empricial false discovery rate.
#' 
#' \code{calculateEFDR} returns the empirical false dicovery rate (EDFR) for supplied
#' thresholds. This function also fits a loess curve to the estimated points. This 
#' allows the calculation of a threshold for priortisation of genes.
#'  
#' @param X A matrix of gene expression values.
#' @param exclude A vector of indices of genes to exclude. 
#' @param index.ref A vector of indices of reference genes used for prioritisation.
#' @param set.size A interger giving the size of the set of genes that are to be
#' prioritised. 
#' @param Weights A object of class \code{Weights} or a list of weights. The weights
#' should correspond to \code{Factor}. If \code{NULL} the unweighted correlations are 
#' used.
#' @param thresholds A vector of thresholds; values should be in the range \eqn{[0,1]}. 
#' @param anno A dataframe or a matrix containing the annotation of arrays in \code{X}.
#' @param Factor A character string corresponding to a column name of \code{anno}.
#' @return \code{calculateEFDR} returns an object of class \code{EFDR}.
#' An object of class \code{EFDR} is a \code{list} with the following components:
#'   \itemize{
#'     \item{\code{EFDR.values}}{ A vector of EDFRs.}
#'     \item{\code{Thresholds}}{ A vector containing the values in \code{threshold}.}
#'     \item{\code{loess.estimate}}{ An object of class \code{loess}.} 
#'   }
#'
#' @details
#' The empirical false discovery rate is estimated by drawing 1000 random sets of genes
#' and calculating how many would be prioritised at every given threshold. A gene is
#' is prioritised if at least one correlation with a known reference gene is above the
#' given threshold.
#'
#' @examples
#' Y<-simulateGEdata(500, 500, 10, 2, 5, g=NULL, Sigma.eps=0.1, 
#' 250, 100, intercept=FALSE, check.input=FALSE)
#' anno<-as.matrix(sample(1:4, dim(Y$Y)[1], replace=TRUE))
#' colnames(anno)<-"Factor"
#' weights<-findWeights(Y$Y, anno, "Factor")
#' calculateEFDR(Y$Y, exclude=251:500, index.ref=1:10, 
#' Weights=weights, anno=anno, Factor="Factor")
#' @author Saskia Freytag
#' @export
#' @exportMethod EFDR
calculateEFDR<-function(X, 
      exclude, 
      index.ref,
      set.size=length(index.ref),
      Weights, 
      thresholds=seq(0.05,1,0.05),
      anno,
      Factor
){

  if(class(Weights)=="Weights") Weights<-Weights$Weights
  if(is.null(Weights)!=TRUE){
    if(length(unique(anno[, 
    which(colnames(anno)==Factor)]))!=length(Weights)){
      stop("The number of Weights should be equal 
      to the number of unique levels in Factor.")
    }
  }
  if(any(thresholds<0)|any(thresholds>1)) { 
    stop("All thresholds need to be in the range [0,1].")
  }
  
  all.genes<-setdiff(1:dim(X)[2], c(exclude,index.ref))
  genes<-lapply(1:1000, function(x) c(index.ref, sample(all.genes, set.size)))
  names(genes)<-1:1000
  
  if(is.null(Weights)!=TRUE){
    data.cor<-lapply(genes, function(x) wcor(X[, x], anno, Factor, Weights))
  } else {
    data.cor<-lapply(genes, function(x) cor(X[, x]))
  }
  
  data.cor<-lapply(1:length(data.cor), function(x) 
  data.cor[[x]][-(1:length(index.ref)), (1:length(index.ref))])
  max.cor<-lapply(data.cor, function(y) apply(y, 1, function(x) max(abs(x))))
  res<-lapply(max.cor, function(y) sapply(thresholds, 
  function(x) length(which(y>x))/set.size))
  
  EFDRs<-Reduce("+", res)/length(genes)
  out<-list(EFDR.values=EFDRs, Thresholds=thresholds, 
  loess.estimate=loess(thresholds~EFDRs))
  class(out)<-"EFDR"
  return(out)
}