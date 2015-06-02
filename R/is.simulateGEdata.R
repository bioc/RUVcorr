#' Checking \code{simulateGEdata} class.
#'
#' \code{is.simulateGEdata} checks if object is of \code{simulateGEdata} class.
#'
#' @param x An object.
#' @return \code{is.simulateGEdata} returns a logical scaler; 
#' \code{TRUE} if the object is of the class \code{simulateGEdata}.
#' @examples
#' \donttest{
#' Y<-simulateGEdata(500, 500, 10, 2, 5, g=NULL, Sigma.eps=0.1, 
#' 250, 100, intercept=TRUE, check.input=TRUE)
#' is.simulateGEdata(Y)
#' }
#' @seealso \code{\link{simulateGEdata}}
#' @author Saskia Freytag
#' @export
is.simulateGEdata<-function(
      x ##object
 ){
  
return(class(x)=="simulateGEdata")
}