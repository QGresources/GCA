#' Cattle dataset 
#'
#' The cattle dataset is simulated with QMSim software (Sargolzaei and Schenkel, 2009).   
#' This dataset includes 2,500 individuals across six generations (from founder to generation 5), 
#' each with 10,000 single nucleotide polymorphisms spread over 29 autosomes. Single phenotype with heritability of 0.6 was simulated.
#' Two fixed covariates of sex and unit are available.
#'
#' @docType data
#' @name GCcattle
#' @keywords datasets
#' @usage 
#' data(GCcattle)
#' 
#' @format 
#' \itemize{
#'   \item cattle.pheno: A matrix with a dimension of 2500 x 6, which includes one phenotype, two fixed covariates and pedigree information.
#'   \item catlle.W: A 2500 by 10,000 matrix, which contains marker information.
#' }
#'
#' @author Haipeng Yu and Gota Morota 
#' 
#' Maintainer: Haipeng Yu \email{haipengyu@@vt.edu}
#'
#' @example man/examples/GCcattle.R 
#' @references Sargolzaei, M., and F. S. Schenkel. 2009. Qmsim: a large-scale
#'   genome simulator for livestock. Bioinformatics 25:680–681. doi:10.1093/bioinformatics/btp045
NULL

