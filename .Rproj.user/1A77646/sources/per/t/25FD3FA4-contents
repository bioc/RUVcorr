#' Optimize parameters of removal of unwanted variation.
#'
#' \code{optimizeParameters} returns the optimal parameters to be 
#' used in the removal of unwanted variation procedure when using simulated data.
#'
#' @param Y An object of the class \code{simualteGEdata}.
#' @param kW.hat A vector of integers for \code{kW} in \code{RUVNaiveRidge}.
#' @param nu.hat A vector of values for \code{nu} in \code{RUVNaiveRidge}.
#' @param nc_index A vector of indices of the negative controls 
#' used in \code{RUVNaiveRidge}.
#' @param check.input Logical; if \code{TRUE} all input is checked; 
#' not advisable for large simulations.
#' @inheritParams assessQuality
#' @param cpus A number specifiying how many workers to use for parallel computing.
#' @param parallel Logical: if \code{TRUE} parallel computing is used.
#'
#' @return \code{optimizeParameters} returns output of the class 
#' \code{optimizeParameters}.
#' An object of class \code{optimizeParameters} is a list containing the
#'  following components:
#'   \describe{
#'     \item{\code{All.results}}{A matrix of output of the quality assessment for all combinations of input parameters.}
#'     \item{\code{Compare.raw}}{A vector of the quality assessment for the uncorrected data.}
#'     \item{\code{Optimal.parameter}}{A matrix or a vector giving the optimal parameter combination.}
#'   }
#'
#' @details
#' The simulated data is cleaned using removal of unwanted variation with all 
#' combinations of the input parameters. The quality of each cleaning is judged by the
#' Frobenius Norm of the correlation as estimated from the cleaned data and the known data
#' or the percentage of correlations with estimated to have the wrong sign. 
#' @examples
#' Y<-simulateGEdata(500, 500, 10, 2, 5, g=NULL, Sigma.eps=0.1, 
#' 250, 100, intercept=FALSE, check.input=FALSE)
#' opt<-optimizeParameters(Y, kW.hat=c(1,5,10), nu.hat=c(100,1000), nc_index=251:500, 
#' methods=c("fnorm"), cpus=1, parallel=FALSE, check.input=TRUE)
#' opt
#' @seealso \code{\link{assessQuality}}, \code{\link{RUVNaiveRidge}}, 
#' \code{\link{funcPara}}
#' @author Saskia Freytag
#' @export
#' @exportClass optimizeParameters
optimizeParameters<-function(
      Y, ##object of the class simulateGEdata
      kW.hat=seq(5,25,5), ##possible kW.hats
      nu.hat=c(0,10,100,1000,10000), ##possible nu.hats
      nc_index, ## index of negative controls for RUV
      methods=c("all", "fnorm", "wrong.sign"), ##methods for quality assessment
      cpus=1,
     parallel=FALSE, 
     check.input=FALSE)
{
  
  if(check.input){
    if(is.simulateGEdata(Y)==FALSE) stop("Y has to be of the class simulateGEdata.")
    if(any(c(kW.hat,nu.hat)<0)) stop("All parameters have to be greater tor equal to 0.")
  }
  
  kW.nu.matrix<-cbind(rep(kW.hat,length=length(nu.hat)*length(kW.hat)),
  rep(nu.hat, each=length(kW.hat)))
  kW.nu.list<-as.list(as.data.frame(t(kW.nu.matrix)))
  index<-1:(dim(Y$Y)[2]-length(nc_index))
  
  if(parallel){
    multicoreParam <- MulticoreParam(workers = cpus)
    res1<-bplapply(kW.nu.list, funcPara, BPPARAM = multicoreParam, Y=Y, 
    nc_index=nc_index, center=TRUE, index=index, methods=methods)
    if(length(methods)==2|any(methods=="all")){
    res<-rbind(sapply(res1, function(x) x[1]), sapply(res1, function(x) x[2]))
    } else {
      res<-as.vector(unlist(res1))
    }
  } else {
    Y.hat<-lapply(kW.nu.list, function(x) RUVNaiveRidge(Y, nc_index, 
    center=TRUE, x[2], x[1], check.input=FALSE))
    res<-sapply(Y.hat, function(x) assessQuality(x, Y$Sigma, index=index, 
    methods=methods))
  }
  
  if(length(methods)==2|any(methods=="all")){
    res<-cbind(kW.nu.matrix, t(res))
    colnames(res)[1:2]<-c("kW", "nu")
  } else {
    res<-cbind(kW.nu.matrix, res)
    colnames(res)<-c("kW", "nu", methods)
  }
  res<-as.data.frame(res)
  
  optimal.parameters<-matrix(0, ncol=2, nrow=2)
  colnames(optimal.parameters)<-c("kW", "nu")
  if(any(methods=="fnorm"| methods=="all")){
    optimal.parameters<-res[which(res$fnorm==min(res$fnorm)), 1:2]
  }
  if(any(methods=="wrong.sign"| methods=="all")){	
    optimal.parameters<-rbind(optimal.parameters, 
    res[which(res$wrong.sign==min(res$wrong.sign)), 1:2])
  }
  if(dim(optimal.parameters)[1]==2){
    rownames(optimal.parameters)<-c("fnorm", "wrong.sign")
  }
  
  compare.raw<-assessQuality(Y$Y, Y$Sigma, index=index, methods=methods)
  
  out<-list(All.results=res, Compare.raw=compare.raw, 
  Optimal.parameters=optimal.parameters)
  class(out)<-"optimizeParameters"
  
  return(out)
}