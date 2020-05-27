#' Variance component estimation 
#'
#' Maximum Likelihood estimation of variance components using the eigenvalues and eigenvectors, 
#'   which are derived from the eigendecomposition of a relationship matrix. 
#' 
#' @param y A vector including the phenotypes of n individuals.  
#' @param Evector A matrix  (n x n) of eigenvectors.
#' @param Evalue A vector of n eigenvalues.
#' 
#' @return A list contains variance components
#' \describe{
#'   \item{$Ve}{An estimate of residual variance.}
#'   \item{$Vu}{An estimate of additive genetic variance.}
#' }
#' 
#' @author Haipeng Yu and Gota Morota 
#' 
#' Maintainer: Haipeng Yu \email{haipengyu@@vt.edu}
#' 
#' @example man/examples/varcomp.R
#' 
#' @export
# The following code is modified from Dr. Gustavo de los Campos: https://github.com/gdlc 
varcomp <- function(y, Evector, Evalue){
  startVal <- log(c(0.5, 0.5))
  var.opt <- optim(fn = log.Lik, y=y, Evector = Evector, Evalue = Evalue, par = startVal,
                   hessian=FALSE) 
  var.est <- list(Ve = exp(var.opt$par[1]) , Vu = exp(var.opt$par[2]))
  return(var.est)
}


# loglikelihood func
log.Lik <- function(y, Evector, Evalue, startVar){
  y <- y - mean(y)
  n <- length(y)
  V <- Evector
  D <- Evalue
  V_y <- crossprod(V,y) 
  V_2y <- as.vector(V_y)^2 
  sigma2e <- exp(startVar[1])
  sigma2a <- exp(startVar[2])
  lambda <- sigma2a/sigma2e 
  DStar <- (D * lambda + 1)
  sumLogD <- sum(log(DStar))
  part1 <- ( n * log(sigma2e) + sumLogD ) 
  part2 <- (sum(V_2y/DStar)) / sigma2e
  LogLik <- part1 + part2
  return(LogLik)
}
