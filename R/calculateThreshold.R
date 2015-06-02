#' Calculates the correlation threshold.
#' 
#' \code{calculateThreshold} returns the proportion of prioritised genes from a random selection 
#' for supplied threshold. Furthermore, this function also fits a loess curve to the estimated points. 
#' This allows the calculation of a threshold for priortisation of genes.
#'  
#' @param X A matrix of gene expression values.
#' @param exclude A vector of indices of genes to exclude. 
#' @param index.ref A vector of indices of reference genes used for prioritisation.
#' @param set.size An integer giving the size of the set of genes that are to be
#' prioritised. 
#' @param Weights A object of class \code{Weights} or a list of weights. The weights
#' should correspond to \code{Factor}. If \code{NULL} the unweighted correlations are 
#' used.
#' @param thresholds A vector of thresholds; values should be in the range \eqn{[0,1]}. 
#' @param anno A dataframe or a matrix containing the annotation of arrays in \code{X}.
#' @param Factor A character string corresponding to a column name of \code{anno}.
#' @param cpus An integer giving the number of cores that are supposed to be used.
#' @param parallel A logical value indicating whether parallel comuting should be used.
#' @return \code{calculateThreshold} returns an object of class \code{Threshold}.
#' An object of class \code{Threshold} is a \code{list} with the following components:
#'   \itemize{
#'     \item{\code{Prop.values}}{ A vector of the proportion of prioritized genes.}
#'     \item{\code{Thresholds}}{ A vector containing the values in \code{threshold}.}
#'     \item{\code{loess.estimate}}{ An object of class \code{loess}.} 
#'   }
#'
#' @details
#' The proportion of prioritized random genes is estimated by drawing 1000 random sets of genes
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
#' calculateThreshold(Y$Y, exclude=251:500, index.ref=1:10, 
#' Weights=weights, anno=anno, Factor="Factor")
#' @author Saskia Freytag
#' @seealso \code{\link{funcThresh}}
#' @export
#' @exportClass Threshold
calculateThreshold<-function(X, 
      exclude, 
      index.ref,
      set.size=length(index.ref),
      Weights=NULL, 
      thresholds=seq(0.05,1,0.05),
      anno=NULL,
      Factor=NULL,
      cpus=1,
      parallel=FALSE
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
  
  if(parallel){
    multicoreParam <- MulticoreParam(workers = cpus)
    res<-bplapply(genes, funcThresh, BPPARAM = multicoreParam, Y=X, 
    Weights=Weights, Factor=Factor, anno=anno, index.ref=index.ref, 
    thresholds=thresholds, set.size=set.size)
  } else {
    res<-lapply(genes, function(x) funcThresh (x, X, Weights, Factor, anno, index.ref, thresholds, set.size))
  }
  
  EFDRs<-Reduce("+", res)/length(genes)
  out<-list(Prop.values=EFDRs, Thresholds=thresholds, 
  loess.estimate=loess(thresholds~EFDRs))
  class(out)<-"Threshold"
  return(out)
}