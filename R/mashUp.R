#' Joining two correlation matrices by diagonal.
#'
#' Internal function that joins two matrices at their diagonal.
#' 
#' @param true Matrix.
#' @param est Matrix.
#' @param plot.genes Vector of indices.
#' @return Matrix.
#' @keywords internal
#' @author Saskia Freytag
mashUp<-function(
      true, ## correlation values
      est, ## matrix with expression values
      plot.genes)
{
  
  est<-cor(est[, plot.genes])
  true<-true[plot.genes, plot.genes]
  est<-est[lower.tri(est)]
  mashup<-true
  mashup[lower.tri(mashup)]<-est
  return(mashup)
}