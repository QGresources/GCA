#' Compute numerator relationship matrix.
#'
#' Use sorted pedigree to calculate a numerator relationship matrix based on a tabular method, where the individuals in \code{Progeny} must be ordered genealogically. 
#'  The missing \code{Sire} and \code{Dam} are coded as 0.
#' 
#' @param Progeny A numeric vector of sorted progenies.
#' @param Sire A numeric vector of sires according to progenies in \code{Progeny} column.
#' @param Dam A numeric vector of dams according to progenies in \code{Progeny} column.
#' @return A n by n numerator relationship matrix, where n refers to total number of individuals in \code{Progeny} column.
#' 
#' @author Haipeng Yu and Gota Morota 
#' 
#' Maintainer: Haipeng Yu \email{haipengyu@@vt.edu}
#' 
#' @example man/examples/computeA.R
#' 
#' @export
computeA <- function(Progeny, Sire, Dam){
  if(!is.numeric(Progeny)) stop('Progeny is not numeric')
  if(!is.numeric(Sire)) stop('Sire is not numeric')
  if(!is.numeric(Dam)) stop('Dam is not numeric')
  if(is.unsorted(Progeny)) stop('Pedigree must be sorted from old to young individuals')
  if (any(duplicated(Progeny))) stop("Progeny must be unique")
  n <- length(Progeny)
  A <- diag(n)
  for (i in 1 : n){
    if (Sire[i] == 0 && Dam[i] != 0){
      for (j in 1 : i-1){
        A[i, j] <- A[j, i] <- 0.5 *(A[j, Dam[i]]) 
      } 
    } else if (Sire[i] != 0 && Dam[i] == 0) {
      for (j in 1 : i-1){
        A[i, j] <- A[j, i] <- 0.5 *(A[j, Sire[i]])
      }
    } else if (Sire[i] != 0 && Dam[i] != 0){
      for (j in 1 : i -1){
        A[i, j] <- A[j, i] <- 0.5 *(A[j, Sire[i]] + A[j, Dam[i]])
      }
      A[i, i] <- A[i, i] + 0.5 * (A[Sire[i], Dam[i]])
    }
  }
  return(A)
}


