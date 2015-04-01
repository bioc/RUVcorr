#' Plot nested design structure.
#' 
#' \code{plotDesign} returns a plot with different color strips representing
#' different factors relating to the study design.
#' genes.
#' 
#' @param anno A dataframe or matrix containing the annotation of the study.
#' @param Factors A vector of factors that should be plotted.
#' @param anno.names A vector containing the names, the default \code{Factors}. 
#' @param orderby A character describing an element in \code{Factor} by which the 
#' data should be ordered.
#' @return \code{plotDesign} returns a plot.
#' 
#' @examples
#' library(bladderbatch)
#' data(bladderdata)
#' expr.meta <- pData(bladderEset)
#' plotDesign(expr.meta, c("cancer", "outcome", "batch"), 
#' c("Diagnosis", "Outcome", "Batch"), orderby="batch")
#' @author Saskia Freytag
#' @export
plotDesign<-function(
      anno, 
      Factors, 
      anno.names=Factors, 
      orderby=NULL)
{
  
  ##testing
  if(any(is.element(Factors,colnames(anno))==FALSE)) stop("Factors need to be in anno.")
  if(is.null(orderby)==FALSE){
    if(is.element(orderby, colnames(anno))==FALSE) stop("orderby need to be in anno.")
  }
  if(length(Factors)!=length(anno.names)){ 
    warning("Not enough names specified in anno.names.")
  }
  
  if(is.null(orderby)==FALSE){
    tmp<-which(is.element(colnames(anno), orderby))
     Order<-order(as.factor(anno[, tmp]))
    anno<-anno[Order, ]
  }
  extract<-vector()
  for(i in Factors){
    extract<-c(extract, which(colnames(anno)==i))
  }
  anno<-anno[, extract]
  a<-dim(anno)[1]
  
  par(mfrow=c(dim(anno)[2], 1), mar=c(1, 2, 0.5, 1), mgp=c(0.5, 1, 0))
  
  for(i in 1:length(anno.names)){
    x<-anno[, i]
    x.name<-anno.names[i]
    x<-as.numeric(as.factor(x))
    n.category<-length(unique(x))
    colors<-hcl(h = seq(0, 360, round(360/n.category, 2))[1:n.category], c=45, l=70)
    
    image(seq(0,a*0.1, 0.1)[1:a], y=1, as.matrix(x), col=colors,
    axes=FALSE, ylab=paste(x.name), xlab="", cex=3)
  }
  
}
