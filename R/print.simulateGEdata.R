#' Print an object of class \code{simualteGEdata}.
#'
#' \code{print.simualteGEdata} is the \code{print} generic 
#' for object so f the class \code{simulateGEdata}.
#'
#' @param x An object of the class \code{simulateGEdata}.
#' @param ... Further arguments passed to or from other methods.
#' @return \code{print.simualteGEdata} returns the information about simulation and
#' the first 5 rows and 5 columns of all matrices.
#' @examples
#' \donttest{
#' Y<-simulateGEdata(500, 500, 10, 2, 5, g=NULL, Sigma.eps=0.1, 
#' 250, 100, intercept=TRUE, check.input=FALSE)
#' Y
#' }
#' @export
#' @author Saskia Freytag
#' @seealso \code{\link{simulateGEdata}}
#' @method print simulateGEdata
print.simulateGEdata<-function(
      x,
      ... 
)
{
  
  cat("Simulated Data:\n")
  cat("Number of samples: ")
  print(dim(x$Y)[1], ...)
  cat("\nNumber of genes: ")
  print(dim(x$Y)[2], ...)
  cat("\nInfo: ")
  print(x$Info)
  
  cat("\n\n Truth\n")
  print(x$Truth[1:5, 1:5], ...)
  cat("\n\n Y\n")
  print(x$Y[1:5, 1:5],...)
  cat("\n\n Noise\n")
  print(x$Noise[1:5, 1:5], ...)
  cat("\n\n Sigma\n")
  print(x$Sigma[1:5, 1:5], ...)
  
}