#' Compare ranking of known reference gene pairs.
#'
#' \code{compareRanks} allows to calculate the difference of the ranks of known reference 
#' gene pairs from two versions of the same data.
#' 
#' @param Y A matrix of raw gene expression values.
#' @param Y.hat A matrix of cleaned gene expression values.
#' @param ref_index A vector of indices that are referrring to genes of interest.
#' @param no.random An integer giving the number of random genes.
#' @param exclude_index A vector of indices to be excluded from the selection of random 
#' genes.
#' @return \code{compareRanks} returns a vector of the differences in ranks of 
#' the correlations of reference gene pairs estimated using raw or cleaned data.
#'
#' @details
#' The correlations between all random genes and reference genes is calculated 
#' (including correlations between random and reference) using the two versions of
#' the data. The correlations are then ranked according to their absolute value (highest
#' to lowest). The ranks of the reference gene pairs are extracted. For a paticular 
#' reference gene pair, the difference in the ranks between the two versions of the data 
#' is calculated:
#' Rank in \code{Y} - Rank in \code{Y.hat}
#' @examples
#' Y<-simulateGEdata(500, 500, 10, 2, 5, g=NULL, 
#' Sigma.eps=0.1, 250, 100, intercept=FALSE, check.input=FALSE)
#' Y.hat<-RUVNaiveRidge(Y, center=TRUE, nu=0, kW=10)
#' compareRanks(Y$Y, Y.hat, ref_index=1:30, no.random=100, exclude_index=c(31:100,251:500))
#' @author Saskia Freytag
#' @export
compareRanks<-function(
      Y,
      Y.hat,
      ref_index,
      no.random=1000,
      exclude_index
      )
{
  
  if(is.null(colnames(Y))){
    colnames(Y)<-1:dim(Y)[2]
  }
  if(is.null(colnames(Y.hat))){
    colnames(Y.hat)<-1:dim(Y.hat)[2]
  }
  
  indices<-(1:dim(Y)[2])[-c(exclude_index, ref_index)]
  ran<-sample(indices, no.random)
  
  cor.Y<-cor(Y[, c(ref_index, ran)])
  tmp.Y<-makeRankedList(cor.Y)
  
  cor.Y.hat<-cor(Y.hat[, c(ref_index, ran)])
  tmp.Y.hat<-makeRankedList(cor.Y.hat)
  
  ref<-colnames(Y)[ref_index]
  rank.Y<-which(is.element(tmp.Y[, 1], ref) & is.element(tmp.Y[, 2], ref))
  
  rank.Y.hat<-which(is.element(tmp.Y.hat[, 1], ref) & is.element(tmp.Y.hat[, 2], ref))
  
  return(as.numeric(rank.Y-rank.Y.hat))
}

