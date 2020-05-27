#' Compute genomic relationship matrix
#'
#' Use single nucleotide polymorphisms markers to derive an additive genomic relationship matrix. Missing markers are allowed, but should be coded as NA. 
#' 
#' @param snpmatrix A marker matrix with the dimension of n by m, where the elements are coded as 0, 1, 2, or NA, 
#'   where n and m are the total number of individuals and markers, respectively.  
#' @param maf A minor allele frequency cutoff for quality control. The default minor allele frequency is 0.05.
#' @param impute Perform genotype imputation for missing markers if applicable. Two methods of 'mean' and 'rbinom' are available, 
#'   where the 'mean' imputes missing markers using mean and 'rbinom' imputes the missing markers by random sampling from 
#'   a binomial distribution. The default method is 'rbinom'. This argument will be ignored if the \code{snpmatrix} does not include any missing markers. 
#' @param method A type of genomic relationship matrix including 'G1' and 'G2' (VanRaden 2008). The default method is 'G1'.  
#' @return An n by n additive genomic relationship matrix. 
#' 
#' @import stats
#' @author Haipeng Yu and Gota Morota 
#' 
#' Maintainer: Haipeng Yu \email{haipengyu@@vt.edu}
#' 
#' @example man/examples/computeG.R
#' 
#' @references \emph{VanRaden, P.M., 2008. Efficient methods to compute genomic predictions. Journal of dairy science, 91(11), pp.4414-4423.}
#' @export
computeG <- function(snpmatrix, maf = 0.05, impute = 'rbinom', method = 'G1') {
  if(anyNA(snpmatrix)) {
    nullid <- unname(which(apply(snpmatrix, 2, anyNA)))
    p_temp <- (colMeans(snpmatrix, na.rm = T) / 2)
    for(j in nullid) {
      if(impute == 'mean') {
        snpmatrix[, j] <- ifelse(is.na(snpmatrix[, j]), 2 * p_temp[, j], snpmatrix[, j])
      } else if(impute == 'rbinom') {
        set.seed(007) # reproducible imputation
        snpmatrix[, j] <- ifelse(is.na(snpmatrix[, j]), rbinom(1, 2, p_temp[j]), snpmatrix[, j])
      }
    }
  }
  p <- colMeans(snpmatrix) / 2
  maf_df <- pmin(p, 1-p)
  maf.index <- which(maf_df < maf)
  ifelse(length(maf.index) < 1, W <- snpmatrix, W <- snpmatrix[, -maf.index])
  if (method == 'G1') { 
    W_c <- scale(W, center = TRUE, scale = FALSE)
    G1 <- tcrossprod(W_c) / sum(2 * p * (1 - p))
    diag(G1) <- diag(G1) + 0.0001 
    cat('Genomic relationship matrix has been computed. Number of SNPs removed:',
        length(maf.index), sep = c(' '))
    return(G1)
  } else if (method == 'G2') {
    W_cs <- scale(W, center = TRUE, scale = TRUE)
    G2 <- tcrossprod(W_cs) / ncol(W_cs)
    diag(G2) <- diag(G2) + 0.0001
    cat('Genomic relationship matrix has been computed. Number of SNPs removed:',
        length(maf.index), sep = c(' '))
    return(G2)
  }
}


