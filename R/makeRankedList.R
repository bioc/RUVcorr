#' Make ranked list of correlations.
#'
#' Internal function.
#'
#' @param Data matrix of gene-gene correlations.
#' @return Matrix.
#' @keywords internal
#' @author Saskia Freytag
makeRankedList<-function(Data){
	reorder<-order(rownames(Data))
	Data<-Data[reorder, reorder]

	Data[lower.tri(Data, diag=TRUE)]<-NA
	a<-melt(Data)
	del<-which(is.na(a[ ,3]))
	a<-a[-del, ]
	a<-a[order(-abs(a[ ,3])), ]
	return(a)
}