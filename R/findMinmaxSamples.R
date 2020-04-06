#' Find minimum and maximum samples in gene expression data.
#' 
#' Internal function that returns 5 samples with smallest inter-quantile range
#' and 5 samples with highest inter-quantile range.
#'
#' @param x Matrix of gene expression values.
#' @return Vector of indices.
#' @keywords internal
#' @author Saskia Freytag
findMinmaxSamples<-function(
      x
){
  
  Summary<-apply(x, 1, summary)
  Summary<-Summary[5, ]-Summary[2, ]
  Summary<-rank(Summary)
  output<-c(which(Summary<=5),which(Summary>(dim(x)[1]-5)))
  Summary<-Summary[output]
  output<-output[order(Summary)]
  return(output)
}