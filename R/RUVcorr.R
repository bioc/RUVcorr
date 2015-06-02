#' Removal of unwanted variation for gene-gene correlations.
#'
#' \pkg{RUVcorr} allows to apply global removal of 
#' unwanted variation (ridged version of RUV) to real and simulated gene expression data. 
#'
#' @details
#' All gene expression data are assumed to be in the following format:
#'   \itemize{
#'     \item Rows correspond to arrays. 
#'     \item Columns correspond to genes.
#'   }
#' @docType package
#' @name RUVcorr
#' @author Saskia Freytag
#' @import corrplot
#' @importFrom MASS mvrnorm 
#' @import grDevices
#' @import psych
#' @import snowfall
#' @import BiocParallel
#' @import bladderbatch
#' @import lattice
#' @importFrom gridExtra grid.arrange
#' @importFrom grid grid.points
#' @importFrom stats ecdf 
#' @importFrom reshape2 melt
NULL