#' Plots an object of class \code{optimizeParameters}.
#'
#' \code{plot.optimizeParameters} generates a heatmap of the quality assessment values
#' stored in the object of class \code{optimizeParameters} .
#'
#' @param x An object of the class \code{optimizeParameters}.
#' @param main A character string describing title of plot.
#' @param ... Further arguments passed to or from other methods.
#' @return \code{plot.optimizeParameters} returns a plot.
#'
#' @details
#' The black point in the heatmap denotes the optimal parameter combination.
#' @examples
#' Y<-simulateGEdata(500, 500, 10, 2, 5, g=2, Sigma.eps=0.1, 
#' 250, 100, intercept=FALSE, check.input=FALSE)
#' opt<-optimizeParameters(Y, kW.hat=c(1,5,10), nu.hat=c(100,100000), 
#' nc_index=251:500, methods=c("fnorm"), cpus=1, parallel=FALSE)
#' try(dev.off(), silent=TRUE)
#' plot(opt, main="Heatmap Plot")
#' @seealso \code{\link{optimizeParameters}}
#' @author Saskia Freytag
#' @export
#' @method plot optimizeParameters
plot.optimizeParameters<-function(
      x, ## an object of the class optimizeParameter
      main=colnames(opt$All.results)[3:dim(opt$All.results)[2]],
      ...
){
  
  var.names<-c("kW", "nu")
  if(is.optimizeParameters(x)==FALSE) {stop("x must be of the class optimizeParameters.")}
  
  if(dim(x$All.results)[2]==4){
    opt1<-x$All.results[, 1:3]
    colnames(opt1)<-c(var.names, "values")
    opt2<-x$All.results[, c(1:2,4)]
    colnames(opt2)<-c(var.names, "values")
    opt.f<-list(opt1,opt2)
  } else{
    opt.f<-x$All.results
    colnames(opt.f)<-c(var.names, "values")
    opt.f<-list(opt.f)
  }
  
  opt<-lapply(opt.f, function(i)
  xtabs(eval(paste("values~", var.names[1], "+", var.names[2], sep="")), data=i))
  
  min.value<-lapply(opt, function(i) which(i==min(i), arr.ind=TRUE)[1,])
  f<-colorRampPalette(c("dodgerblue4", "dodgerblue1", "cornflowerblue", 
  "lightblue2", "white", "lightcoral", "indianred1", "firebrick", "indianred4")) 
  plotL<-lapply(1:length(opt), function(i) levelplot(opt[[i]], 
  scales=list(x=list(rot=90)), col.regions=f(100), main=main[i], colorkey=list(width=0.8),
  panel=function(...){
    panel.levelplot(...)
     grid.points(min.value[[i]][1], min.value[[i]][2], pch=19) 
  }))
  if(length(opt)==2){
    grid.arrange(plotL[[1]], plotL[[2]], ncol=2)} else {
      print(plotL)
    }
}