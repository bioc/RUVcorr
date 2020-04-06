#' Find the inter quantile range.
#'
#' Internal function to find the inter quantile range.
#'
#' @param Vector of gene expression values.
#' @return Numeric value.
#' @keywords internal
#' @author Saskia Freytag
findIQR<-function(
      x ##vector of expression data
){
  
  tmp<-summary(x)
  return(tmp[5]-tmp[2])
}