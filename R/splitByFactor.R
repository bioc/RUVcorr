#' Splitting a data set by a particular factor.
#'
#' Internal function that splits a data set according to a particular factor.
#' 
#' @param X  A matrix containing gene expressions.
#' @param anno A dataframe or a matrix containing the annotation of arrays in \code{X}.
#' @param Factor A character string corresponding to a column name of 
#' \code{anno} to be used for splitting. 
#' @return \code{splitByFactor} returns a list object.
#' @keywords internal
#' @author Saskia Freytag
splitByFactor<-function(
      X,
      anno, 
      Factor
      )
{
  
  expression<-list()
  col.num<-which(colnames(anno)==Factor)
  categories<-sort(unique(anno[, col.num]))
  ind<-lapply(categories, function(x) which(anno[, col.num]==x))
  expression<-lapply(ind, function(x) X[x, ])
  names(expression)<-categories
  return(expression)
  
}
