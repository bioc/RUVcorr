#' Plots different versions of relative log expression plots
#'
#' \code{RLEPlot} generates three different types of 
#' relative log expression plots for high-dimensional data.
#'
#' @param X A matrix of gene expression values.
#' @param Y A matrix of gene expression values.
#' @param center A logical scalar; \code{TRUE} if centering should be applied.
#' @param name A vector of characters describing the data contained in 
#' \code{X} and \code{Y}.
#' @param title A character string describing the title of the plot.
#' @param method The type of RLE plot to be displayed; possible inputs are 
#' \code{"IQR.points"}, \code{"IQR.boxplots"} and \code{"minmax"} 
#' (for information see details).
#' @param anno A dataframe or a matrix containing the annotation of 
#' arrays in \code{X} and \code{Y} (only applicable for \code{method="IQR.points"});
#' if \code{anno=NULL} data points are not colored.
#' @param Factor A character string corresponding to a column name of 
#' \code{anno} to be used for coloring.
#' @param numeric A logical scalar indicating whether \code{Factor} is numerical.
#' @param new.legend A vector describing the names used for labelling; if \code{NULL} 
#' labels in \code{Factor} are used. 
#' @param outlier A logical indicating whether outliers should be plotted; only 
#' applicable when \code{method="minmax"}.
#' @return \code{RLEPlot} returns a plot.
#' 
#' @details
#' There are three different RLE plots that can be generated using \code{RLEPlot}:
#'    \describe{
#'      \item{\code{"IQR.points"}}{Median expression vs. inter-quantile range of every array.}
#'      \item{\code{"IQR.boxplots"}}{Boxplots of the 25\% and 75\% quantile of all arrays.}
#'      \item{\code{"Minmax"}}{Ordinary RLE plots for the 5 arrays with the smallest and largest inter-quantile ranges.}
#'    }
#' Note that normal RLE plots are not supplied as they 
#' are not very suitable for high-dimensional data.
#' @examples
#' Y<-simulateGEdata(500, 500, 10, 2, 5, g=NULL, Sigma.eps=0.1, 
#' 250, 100, intercept=FALSE, check.input=FALSE)
#' Y.hat<-RUVNaiveRidge(Y, center=TRUE, nc_index=251:500, 0, 10, check.input=TRUE)
#' try(dev.off(), silent=TRUE)
#' par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
#' RLEPlot(Y$Y, Y.hat, name=c("Raw", "RUV"), title="", method="IQR.points")
#' try(dev.off(), silent=TRUE)
#' par(mfrow=c(1, 1))
#' RLEPlot(Y$Y, Y.hat, name=c("Raw", "RUV"), title="", method="IQR.boxplots")
#' try(dev.off(), silent=TRUE)
#' RLEPlot(Y$Y, Y.hat, name=c("Raw", "RUV"), title="", method="minmax")
#'
#' #Create a random annotation file
#' anno<-as.matrix(sample(1:4, dim(Y.hat)[1], replace=TRUE))
#' colnames(anno)<-"Factor"
#' try(dev.off(), silent=TRUE)
#' RLEPlot(Y$Y, Y.hat, name=c("Raw", "RUV"), title="", method="IQR.points",
#' anno=anno, Factor="Factor", numeric=TRUE)
#' @author Saskia Freytag, Terry Speed
#' @export
RLEPlot<-function( 
      X, ## old matrix of gene expression
      Y, ## new matrix of gene expression
      center=TRUE,
      name,
      title,
      method=c("IQR.points", "IQR.boxplots", "minmax"),
      anno=NULL,
      Factor=NULL,
      numeric=FALSE,
      new.legend=NULL,
      outlier=FALSE
){
  
  if(any(is.null(Factor)|is.null(anno))==FALSE){
    if(is.null(anno)|is.null(Factor)) { 
      stop("Both anno and factor have to be specified.")
    }
    if(any(colnames(anno)==Factor)==FALSE) {
      stop("Factor has to be a column of anno.")
    }
    if(dim(anno)[1]!=dim(Y)[1]) {
      stop("Annotation does not match.")
    }
  }
  
  if(center){
    Y<-scale(Y, center=TRUE, scale=FALSE)
    X<-scale(X, center=TRUE, scale=FALSE)
  }
  
  if(method[1]!="IQR.points"){
    Y.summary<-apply(Y, 1, summary)
    X.summary<-apply(X, 1, summary)
  }
 
  if(method[1]=="IQR.points"){
    Y.IQR<-apply(Y, 1, function(x) findIQR(x))
    X.IQR<-apply(X, 1, function(x) findIQR(x))
    
    Y.Median<-apply(Y,1, function(x) median(x))
    X.Median<-apply(X,1, function(x) median(x))
    
    if((is.null(anno)|is.null(factor))==FALSE){
      category<-as.matrix(as.character(anno[, which(colnames(anno)==Factor)]))
      category<-cbind(category, 1:length(category))
      colnames(category)<-c("x", "Ind")
      n.category<-length(unique(anno[, which(colnames(anno)==Factor)]))
      
      colours<-hcl(h = seq(0, 360,round(360/n.category, 2)), c=45, l=70)
      if (numeric==TRUE){
        colour.code<-cbind(unique(as.character(
        category[, 1]))[order(unique(as.numeric(category[, 1])))], colours[1:n.category])
      } else {colour.code<-cbind(sort(unique(as.character(category[,1])), 
      na.last=TRUE), colours[1:n.category])}
      colnames(colour.code)<-c("x", "colour")
      new.colours<-merge(category, colour.code, by="x")
      new.colours<-new.colours[order(as.numeric
      (as.character(new.colours[ ,2]))),]
      new.colours<-new.colours[ ,3]
      
      layout(matrix(c(1, 2, 3), nrow=1), widths=c(2, 4, 4))
      par(mar = c(3, 2.5, 2, 1), mgp = c(1.3, 0.45, 0))
      plot.new()
      if(is.null(new.legend)){
      legend("topleft", legend=colour.code[,1], col=colour.code[,2], ncol=1, 
      bty="n", pch=16, cex=0.75)
      } else{
        legend("topleft", legend=new.legend, col=colour.code[,2], 
        ncol=1, bty="n", pch=16, cex=0.75)
      }
      mtext(paste(title), side=3, line=0, cex = 1)
      xlim<-c(min(c(X.Median, Y.Median)),max(c(X.Median,Y.Median)))
      ylim<-c(min(c(X.IQR, Y.IQR)),max(c(X.IQR,Y.IQR)))
      plot(X.Median, X.IQR , pch=1, main=paste(name[1]), 
      col=as.character(new.colours), 
      xlab="Median", ylab="IQR", bty="l", xlim=xlim, ylim=ylim)
      plot(Y.Median, Y.IQR ,pch=1, main=paste(name[2]), 
      col=as.character(new.colours), 
      xlab="Median", ylab="IQR", bty="l", xlim=xlim, ylim=ylim)
      
    } else {
      par(mfrow=c(1,2))
      plot(X.Median, X.IQR ,pch=19, main=paste(name[1]), col="gray", 
      xlab="Median", ylab="IQR", bty="l")
      plot(Y.Median, Y.IQR ,pch=19, main=paste(name[2]), col="gray", 
      xlab="Median", ylab="IQR", bty="l")
      mtext(paste(title), side=3, line=-1, cex = 1)
    }
    
  }
  
  if(method[1]=="IQR.boxplots"){
    boxplot(t(X.summary[c(2,5), ]), border="red", main=paste(title), las=2)
    abline(h=0, lwd=2, lty=2, col="grey", bty="l")
    boxplot(t(Y.summary[c(2, 5), ]), border="black", add=TRUE, las=2, boxwex=0.5)
    legend("bottomright", legend=name ,text.col=c("red", "black"), bty="n")
  }
  
  if(method[1]=="minmax"){
    index.Y<-findMinmaxSamples(X)
    boxplot(t(X[index.Y, ]), border=c(rep("tomato", 5), rep("red", 5)), 
    las=2, outline=outlier, main=title, bty="l")
    boxplot(t(Y[index.Y, ]), border=c(rep("gray30",5), rep("black", 5)), 
    main=paste(title), las=2, outline=outlier, boxwex=0.5, add=TRUE)
    
    legend("topleft", legend=name ,text.col=c("red", "black"), bty="n")
    abline(h=0, lwd=2, lty=2, col="grey")
  }
  
}