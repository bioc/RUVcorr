#' Plot principle component analysis for gene expression data.
#'
#' \code{PCAPlot} generates principle component plots for with the possibility
#' to color arrays according to a known factor.
#'
#' @param Y A matrix of gene expression values or an object of class \code{prcomp}.
#' @param comp A vector of length 2 specifying which principle components to be used.
#' @param anno A dataframe or a matrix containing the annotation of the arrays.
#' @param Factor A character string describing the column name of 
#' \code{anno} used for coloring.
#' @param numeric A logical scalar indicating whether \code{Factor} is numerical.
#' @param new.legend A vector describing the names used for labelling; if \code{NULL} 
#' labels in \code{Factor} are used. 
#' @param title A character string giving the title. 
#' @return \code{PCAPlot} returns a plot.
#' @examples
#' Y<-simulateGEdata(500, 500, 10, 2, 5, g=NULL, Sigma.eps=0.1, 
#' 250, 100, intercept=FALSE, check.input=FALSE)
#' PCAPlot(Y$Y, title="")
#' 
#' ## Create random annotation file
#' anno<-as.matrix(sample(1:4, dim(Y$Y)[1], replace=TRUE))
#' colnames(anno)<-"Factor"
#' try(dev.off(), silent=TRUE)
#' par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0, mfrow=c(1, 1))
#' PCAPlot(Y$Y, anno=anno, Factor="Factor", numeric=TRUE, title="")
#' @seealso \code{\link[stats]{prcomp}}
#' @author Saskia Freytag
#' @export
PCAPlot<-function(
      Y, ## gene expression data or principle components
      comp=c(1, 2), ## Vector of components to plot
      anno=NULL, 
      Factor=NULL, 
      numeric=FALSE, 
      new.legend=NULL,
      title) 
{
  
  if(any(is.null(Factor)|is.null(anno))==FALSE){
    if(is.null(anno)|is.null(Factor)) stop("Both anno and factor have to be specified.")
    if(any(colnames(anno)==Factor)==FALSE) stop("Factor has to be a column of anno.")
    if(class(Y)!="prcomp"){
      if(dim(anno)[1]!=dim(Y)[1]) stop("Annotation does not match.")
    } else {
      if(dim(anno)[1]!=dim(Y$x)[1]) stop("Annotation does not match.")
    }
  }
  
  if(class(Y)!="prcomp"){
    pc<-prcomp(Y, center=TRUE)
  } else{
    pc<-Y
  }
  
  print("Calculation of principle components finished. Start plotting...")
  
  if(is.null(anno)==FALSE){
    category<-cbind(as.matrix(as.character(anno[, which(colnames(anno)==Factor)])), 
    1:dim(pc$x)[1])
    colnames(category)<-c("x", "ord")
    n.category<-length(unique(anno[, which(colnames(anno)==Factor)]))
    
    colours<-hcl(h = seq(0,360, round(360/n.category, 2)), c=45, l=70)
    if (numeric==TRUE){
      colour.code<-cbind(unique(as.character(category[, 1]))[order(
      unique(as.numeric(category[, 1])))], colours[1:n.category])
     } else {
       colour.code<-cbind(sort(unique(as.character(category[,1])), na.last=TRUE), 
       colours[1:n.category])
    }
    colnames(colour.code)<-c("x", "colour")
    new.colours<-merge(category, colour.code, by="x")
    new.colours<-new.colours[, 3]
    
    layout(matrix(c(1, 2), nrow=1), widths=c(1, 3))
    par(mar = c(3.5, 2.2, 2, 1), mgp = c(1.3, 0.45, 0))
    plot(0, type='n', axes=FALSE, ann=FALSE)
    if(is.null(new.legend)){
    legend("topleft", legend=colour.code[, 1], col=colour.code[, 2], 
    ncol=1, bty="n", pch=16, cex=0.75)
    } else {
      legend("topleft", legend=new.legend, col=colour.code[, 2], 
      ncol=1, bty="n", pch=16, cex=0.75)
    }
    plot(pc$x[, comp[1]], pc$x[, comp[2]], col=as.character(new.colours), 
    xlab=paste("PC", comp[1]), ylab=paste("PC", comp[2]), bty="l")
    mtext(text=paste(title), outer=TRUE, side=3, line=-1.5, font=2, cex=1.25)	
    mtext(paste("St. Dev PC", comp[1], "=", round(pc$sdev[comp[1]], 2), 
    " St. Dev PC", comp[2], "=", round(pc$sdev[comp[2]], 2), sep=""),
    outer=TRUE, side=1, line=-1, font=2, cex=0.75)
  }
  
  if(is.null(anno)){
    plot(pc$x[, comp[1]],  pc$x[, comp[2]], xlab=paste("PC", comp[1]),
    ylab=paste("PC", comp[2]), bty="l")
    mtext(text=paste(title), outer=TRUE, side=3, line=-2, font=2, cex=1.25)
    mtext(paste("St. Dev PC", comp[1], "=", round(pc$sdev[comp[1]],2), " St. Dev PC", 
    comp[2], "=", round(pc$sdev[comp[2]],2), sep=""), outer=TRUE, side=1, 
    line=-1, font=2, cex=0.75)
  }
}