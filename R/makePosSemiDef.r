#' Makes square matrices positive semi-definite.
#' 
#' Internal function which returns closest positive semi-definite matrix to input matrix.
#' 
#' @param a Square matrix.
#' @param offset Offset.
#' @return Positive semi-definite matrix.
#' @keywords internal
#' @author Saskia Freytag
makePosSemiDef<-function(
      a, # matrix 
      offset=0)
{
  
  lambda<-abs(min(eigen(a, only.values=TRUE)$values))
  return((a + (diag(dim(a)[1]) * (lambda+offset))) /(1+lambda+offset))
}
