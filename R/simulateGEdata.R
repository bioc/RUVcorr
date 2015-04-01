#' Simulate gene expression data.
#'
#' \code{simulateGEdata} returns simulated noisy gene expression values of specified size 
#' and its underlying gene-gene correlation.
#'
#' @param n An integer setting the number of genes.
#' @param m An integer setting the number of arrays.
#' @param k An integer setting number of dimensions of noise term, 
#' controls dimension of \eqn{W} and \eqn{\alpha}.
#' @param size.alpha A numeric scalar giving the maximal and 
#' minimal absolute value of \eqn{\alpha}.
#' @param g An integer value between [1, min(\code{k}, \code{corr.strength})) giving the 
#' correlation between \eqn{X} and \eqn{W} or \code{NULL} for independence.
#' @param corr.strength An integer controlling the dimension of \eqn{X} and \eqn{\beta}.
#' @param Sigma.eps A numeric scalar setting the amount of random variation in 
#' \eqn{\epsilon}; \code{Sigma.eps} \eqn{>0}.
#' @param nc An integer setting the number of negative controls.
#' @param ne An integer setting the number of strongly expressed genes.
#' @param intercept An logical value indicating whether the systematic noise has an intercept.
#' @param check.input A logical scalar; if \code{TRUE} all input is checked 
#' (not advisable for large simulations).
#' @return \code{simulateGEdata} returns output of the class \code{simulateGEdata}.
#' An object of class \code{simulateGEdata} is a \code{list} with the 
#' following components:
#'     \itemize{
#'       \item{\code{Truth}}{ A matrix containing the values of \eqn{X\beta}.}
#'       \item{\code{Y}}{ A matrix containing the values in \eqn{Y}.}
#'       \item{\code{Noise}}{ A matrix containing the values in \eqn{W\alpha}.} 
#'       \item{\code{Sigma}}{ A matrix containing the true gene-gene correlations, as defined by \eqn{X\beta}.}
#'       \item{\code{Info}}{ A matrix containing some of the general information about the simulation.}
#'     }
#'
#' @details
#' This function generates log2-transformed expression values of \code{n} genes in 
#' \code{m} arrays. The expression values consist of true expression and noise:
#' \deqn{Y=X\beta+W\alpha+\epsilon}
#' The dimensions of the matrices \eqn{X} and \eqn{\beta} are used to control the size of
#' the correlation between the genes. It is possible to simualte three different classes
#' of genes:
#'     \itemize{
#'       \item correlated genes expressed with true log2-transformed values from 0 to 16
#'       \item correlated genes expressed with true log2-transformed values with mean 0
#'       \item uncorrelated genes with true log2-transformed expression equal to 0 (negative controls)
#'    }
#' The negative control are always the last \code{nc} genes in the data, 
#' whereas the strongly expressed genes are always the first \code{ne} genes in the data.
#' The parameter \code{intercept} controls whether the systematic noise has an 
#' offset or not. Note that the intercept is one dimension of \eqn{W}.
#' It is possible to either simulate data where \eqn{W} and \eqn{X} are independent by
#' setting \code{g} to NULL, or increasing correlation \eqn{bWX} between 
#' \eqn{W} and \eqn{X} by increasing \code{g}.
#' @examples 
#' Y<-simulateGEdata(500, 500, 10, 2, 5, g=NULL, Sigma.eps=0.1, 
#' 250, 100, intercept=TRUE, check.input=TRUE)
#' Y
#' Y<-simulateGEdata(500, 500, 10, 2, 5, g=3, Sigma.eps=0.1, 
#' 250, 100, intercept=TRUE, check.input=TRUE)
#' Y
#' @references Jacbo L., Gagnon-Bartsch J., Speed T. Correcting gene expression 
#'   data when neither the unwanted variation nor the factor of interest are observed. 
#'   Berkley Technical Reports (2012).
#' @author Saskia Freytag, Johann Gagnon-Bartsch 
#' @exportMethod simulateGEdata
#' @exportClass simulateGEdata
#' @export
 simulateGEdata<-function (
       n, 
       m, 
       k, 
       size.alpha, 
       corr.strength, 
       g = NULL, 
       Sigma.eps = 0.1,
       nc, 
       ne, 
       intercept=TRUE, 
       check.input = FALSE)
{
  if (check.input) {
    if (any(c(n, m, k, size.alpha, corr.strength, Sigma.eps,
    nc, ne) < 0)) {
      stop("All input variables have to be positive!")
    }
    if (nc + ne > n) {
      stop("The total number of genes has to be greater or equal\n\t\tthan 
      the number negative controls plus the number of expressed genes.")
    }
    if (corr.strength%%1 != 0) {
       stop("The variable corr.strength needs to be\n\t\tan integer!")
    }
    if (k > m) {
      warning("The number of dimensions included in the noise 
      is\n\t\tgreater than the number of arrays.")
    }
  }
  
  
  if (is.null(g)) {
    beta.e <- matrix(runif(corr.strength * (n - nc), -2/sqrt(corr.strength),   
    2/sqrt(corr.strength)), nrow = corr.strength, ncol = n - nc)
    beta <- cbind(beta.e, matrix(0, nrow = corr.strength,
    ncol = nc))
    X <- matrix(rnorm(m * corr.strength, 0, (2/corr.strength)),
    nrow = m, ncol = corr.strength)
    X.beta <- X %*% beta
    X.beta <- X.beta + matrix(rep(c(runif(ne, 0, 16), rep(0,
    n - nc - ne), rep(0, (nc))), each = m), nrow = m,
    ncol = n)
    
    if(intercept){
      W <- matrix(rnorm(m * (k-1), 0, 1), nrow = m, ncol = (k-1))
      W<-cbind(1, W)
     } else {
      W <- matrix(rnorm(m * k, 0, 1), nrow = m, ncol = k)
     }
     
  } else {
    if (check.input) {
      if (g<=0 |g>=min(corr.strength, k)| g%%1 !=0) { 
       stop("g is not an integer or in the range 0<g<=min(corr.strength, k)")
      }
    }
    beta.e <- matrix(runif(corr.strength * (n - nc), -2/sqrt(corr.strength),
    2/sqrt(corr.strength)), nrow = corr.strength, ncol = n - nc)
    beta <- cbind(beta.e, matrix(0, nrow = corr.strength,
     ncol = nc))
     
    Lg1<-cbind(diag(g), matrix(0,nrow=g, ncol=(k-g)))
    Lg2<-matrix(0, nrow=corr.strength-g, ncol=k)
    Lg<-rbind(Lg1, Lg2)
    
    Sigma.XW.1<-cbind(diag(corr.strength), Lg)
    Sigma.XW.2<-cbind(t(Lg), diag(k))
    Sigma.XW<-rbind(Sigma.XW.1, Sigma.XW.2)
    
    if (min(eigen(Sigma.XW, only.values=TRUE)$values) <= 0) {
      print("Need to make positive semi-definite!")
      Sigma.XW <- makePosSemiDef(Sigma.XW, offset = 5e-04)
    }
    X.W <- mvrnorm(m, rep(0, (corr.strength + k)), Sigma.XW)
    
    X.beta <- X.W[, 1:corr.strength] %*% beta
    X.beta <- X.beta + matrix(rep(c(runif(ne, 0, 16), rep(0,
    n - nc - ne), rep(0, (nc))), each = m), nrow = m,
     ncol = n)
    if(intercept){
      W <- X.W[, (corr.strength + 1):(corr.strength + (k-1))]
      W<-cbind(1, W)
    } else {
      W <- X.W[, (corr.strength + 1):(corr.strength + k)]
    }
  }
  
  eps <- matrix(rnorm(m * n, 0, Sigma.eps), nrow = m, ncol = n)
  Sigma <- diag(n)
  Sigma.tmp <- cor(X.beta[, 1:(n - nc)])
  Sigma[1:(n - nc), 1:(n - nc)] <- Sigma.tmp
  alpha <- matrix(runif(n * k, -2*size.alpha/sqrt(k), 2*size.alpha/sqrt(k)), ncol = n,
   nrow = k)
  noise <- W %*% alpha
  if (is.null(g)) {
    info <- cbind(c("k", "Mean correlation", "Size alpha", "Intercept"),
    c(k, round(mean(abs(Sigma.tmp[lower.tri(Sigma.tmp)])),
    5), size.alpha, as.logical(intercept)))
  }
  else {
    suppressWarnings(bWX.tmp<-cor(X.beta,noise))
    suppressWarnings(info <- cbind(c("k", "Mean correlation", "bWX", 
    "Size alpha", "Intercept"),
    c(k, round(mean(abs(Sigma.tmp[lower.tri(Sigma.tmp)])),
    5), round(mean(abs(bWX.tmp[lower.tri(bWX.tmp)]), na.rm=TRUE),5), 
    size.alpha, as.logical(intercept))))
  }
  res <- list(Truth = X.beta, Y = X.beta + noise + eps, Noise = noise,
  Sigma = Sigma, Info = info)
  class(res) <- "simulateGEdata"
  
  return(res)
}
