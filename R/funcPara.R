#' Function to optimize parameters in parallel.
#'
#' Internal function for parallel computing.
#' 
#' @param x Vector.
#' @param Y \code{simulateGE} object.
#' @param nc_index Vector.
#' @param center Logical.
#' @param index Vector.
#' @param methods Vector.
#' @return List.
#' @keywords internal
#' @author Saskia Freytag
funcPara<-function(x, 
      Y, 
      nc_index,
      center=TRUE,
      index,
      methods)
{
  Y.hat.tmp<-RUVNaiveRidge(Y, nc_index, center=TRUE, x[2], x[1], check.input=FALSE)
  assessQuality(Y.hat.tmp, Y$Sigma, index=index, methods=methods)
}