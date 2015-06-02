#' Checking \code{optimizeParameters} class.
#'
#' \code{is.optimizeParameters} checks if object is of \code{optimizeParameters} class.
#'
#' @param x An object.
#' @return \code{is.optimizeParameters} returns a logical scalar; 
#' \code{TRUE} if the object is of the class \code{optimizeParameters}.
#' @examples
#' \donttest{
#' Y<-simulateGEdata(500, 500, 10, 2, 5, g=NULL, Sigma.eps=0.1, 
#' 250, 100, intercept=FALSE, check.input=FALSE)
#' opt<-optimizeParameters(Y, kW.hat=c(1,5,10), nu.hat=c(100,1000), 
#' nc_index=251:500, methods=c("fnorm"), cpus=1, parallel=FALSE)
#' opt
#' is.optimizeParameters(opt)
#' }
#' @seealso \code{\link{optimizeParameters}}
#' @author Saskia Freytag
#' @export
is.optimizeParameters<-function(
       x ##object
){
  
  return(class(x)=="optimizeParameters")
}