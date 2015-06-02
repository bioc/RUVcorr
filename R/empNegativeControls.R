#' Empirically choose negative control genes. 
#' 
#' \code{empNegativeControls} finds suitable negative controls in real or simulated data.
#'
#' @param Y A matrix of gene expression values or an 
#'object of the class \code{simulateGEdata}.
#' @param exclude A vector of indices to be excluded from being chosen
#'  as negative controls.
#' @param smoothing A numerical scalar determining the amount of smoothing to be applied.
#' @param nc An integer setting the number of negative controls.
#' @return \code{empNegativeControls} returns a vector of indicies of 
#' empirically chosen negative controls.
#'
#' @details
#' First the mean of all genes (except the excluded genes) is calculated 
#' and genes are accordingly assigned to bins. The bins have the size
#' of the smoothing parameter. In each bin the function picks a number of negative control 
#' genes proportional to the total number of genes in the bin. 
#' The picked genes in each bin have the lowest inter-quantile ranges of all 
#' genes in the respective bin.
#' @examples
#' Y<-simulateGEdata(500, 500, 10, 2, 5, g=NULL, Sigma.eps=0.1, 
#' 250, 100, intercept=FALSE, check.input=TRUE)
#' empNegativeControls(Y, exclude=1:100, nc=100)
#' @section Warning:
#' For simulated data it is advisable to use the known negative controls or 
#' restrict the empirical
#' choice to the known negative controls by excluding all other genes.
#' @author Saskia Freytag
#' @exportMethod empNegativeControls
#' @export
empNegativeControls<-function(Y, ##matrix of expression data
      exclude, ##index of genes to be excluded from negative controls
      smoothing=0.1, ## smoothing factor
      nc ##number of negative controls
      ) UseMethod("empNegativeControls")
  
#' \code{empNegativeControls.default} empirically chooses negative 
#' control genes for matrix input.
#'
#' @rdname empNegativeControls
#' @export 
empNegativeControls.default<-function(
      Y, ## matrix of expression data
      exclude, ##index of genes to be excluded from negative controls
      smoothing=0.1, ## smoothing factor
      nc ##approx number of negative controls
){
  if(is.null(colnames(Y))) colnames(Y)<-1:dim(Y)[2]
  tmpY<-Y[, -exclude]
  prop<-nc/dim(tmpY)[2]
  ## find the proportion of negative controls to all genes
  
  tmpMeans<-colMeans(tmpY)
  ## find means for every gene
  
  tmpIQR<-apply(tmpY, 2, function(x) findIQR(x))
  ## find IQR for every gene
  
  empNC<-vector(length=0)
  ## initialise vector of empirical negative controls
  
  partition<-seq(floor(min(tmpMeans)), ceiling(max(tmpMeans)), smoothing)
  ## set partition of data
  
  for(i in 1:(length(partition)-1)){
    start<-partition[i]
    end<-partition[i+1]
    index<-which(tmpMeans<end & tmpMeans>start)
    ## find genes in partition
    if(length(index)==0) next
    ## if no genes are in the partition go to next
    size<-ceiling(length(index)*prop)
    ## find how many genes to pick
    tmpIQRsorted<-sort(tmpIQR[index], decreasing=FALSE)
    empNC<-c(empNC, names(tmpIQRsorted)[1:size])
  }
  
  empNC<-empNC[-sample(1:length(empNC), (length(empNC)-nc))]
  # randomly remove too many controls
  
  return(which(is.element(colnames(Y), empNC)))
}

#' \code{empNegativeControls.simulateGEdata} empircially chooses negative control genes for \code{simulateGEdata} object. 
#'
#' @rdname empNegativeControls
#' @export
empNegativeControls.simulateGEdata<-function(
      Y, ##object of class simualteGEdata
      exclude, ##index of genes to be excluded from negative controls
      smoothing=0.1, ## smoothing factor
      nc ##approx number of negative controls
){
  
  if(is.simulateGEdata(Y)==FALSE) stop("Y needs to be of class simulateGEdata.")	
  Y<-Y$Y
  colnames(Y)<-1:dim(Y)[2]
  
  empNegativeControls.default(Y, exclude, smoothing, nc)
  
}