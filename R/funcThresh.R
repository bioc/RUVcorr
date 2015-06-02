#' Function to calculate correlation threshold in parallel.
#'
#' Internal function for parallel computing.
#' 
#' @param .x Vector.
#' @param Y Matrix.
#' @param Weights A object of class \code{Weights} or a list of weights.
#' @param Factor Character string.
#' @param anno Dataframe.
#' @param index.ref Vector.
#' @param thresholds Vector.
#' @param set.size Integer.
#' @return Matrix.
#' @keywords internal
#' @author Saskia Freytag
funcThresh<-function(.x, Y, Weights, Factor, anno, index.ref, thresholds, set.size)
{
  
  if(is.null(Weights)!=TRUE){
    data.cor<-wcor(Y[, .x], anno, Factor, Weights)
  } else {
    data.cor<-cor(Y[, .x])
  }
  
  data.cor<-data.cor[-(1:length(index.ref)), (1:length(index.ref))]
  max.cor<-apply(data.cor, 1, function(a) max(abs(a)))
  res<-sapply(thresholds, function(a) length(which(max.cor>a))/set.size)
  return(res)
}