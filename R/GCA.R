#' \strong{GCA}: Genetic connectedness analysis
#' 
#' @description 
#' An R package for genetic connectedness analysis across units using pedigree and genomic data.
#' 
#' @details 
#'The GCA package encompasses numerous connectedness statistics, which could be labeled as two groups, 
#'by reference to connectedness based on prediction error variance (PEV) and variance of unit effect estimates (VE).
#'The PEV-derived metrics include prediction error variance of differences (PEVD), coefficient of determination (CD), 
#'and prediction error correlation (r). These PEV-derived metrics can be summarized at the unit level as the average PEV within
#'and across units (GrpAve), average PEV of all pairwise differences between individuals across units (IdAve), 
#'or using a contrast vector (Contrast). VE-derived metrics comprise variance of differences in management unit effects (VED), 
#'coefficient of determination of VED (CDVED), and connectedness rating (CR). Three correction factors accounting for the number 
#'of fixed effects can be applied for each VE-derived metric. These include non-correction (0), correction of unit effect (1), 
#'and correction of two or more fixed effects (2). The core function of GCA is integrated with C++ to improve computational 
#'efficiency using the Rcpp package (Eddelbuettel and François 2011). The details of these connectedness statistics can be found
#'in Yu and Morota 2019.
#' 
#' @section Available functions in GCA package:
#' \itemize{
#'   \item computeA(): Computation of numerator relationship matrix. 
#'   \item computeG(): Computation of genomic relationship matrix.
#'   \item gca(): Measures of genetic connectedness.
#'   \item varcomp(): Estimates of variance components using eigenvalues and eigenvectors. 
#' }
#'
#' @useDynLib GCA, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' 
#' @author Haipeng Yu and Gota Morota 
#' 
#' Maintainer: Haipeng Yu \email{haipengyu@@vt.edu}
#' @references \emph{Eddelbuettel D, François R (2011). Rcpp: Seamless R and C++ Integration. Journal of Statistical Software, 40(8), 1–18.}
#' @references \emph{Yu H and Morota G. 2019. GCA: An R Package for Genetic Connectedness Analysis Using Pedigree and Genomic Data.}
#' 
#'@keywords internal
"_PACKAGE"

#' Genetic connectedness analysis
#'
#' Estimates genetic connectedness across units using pedigree or genomic data.
#' 
#' @param Kmatrix A relationship matrix with a dimension of n by n, where n refers to the total number of individuals.  
#' @param Xmatrix A design matrix which associates fixed effects with phenotypes and the intercept is excluded. 
#'   The first column of \code{Xmatrix} should start with design matrix of unit effects, following with other fixed effects if applicable. 
#' @param sigma2a Additive genetic variance.
#' @param sigma2e Residual variance.
#' @param MUScenario A vector of units which will be treatd as a factor. 
#' @param statistic A statistic measures genetic connectedness, which includes 
#' \itemize{
#'   \item 'PEVD_IdAve' : Individual average PEVD, the optional argument of 'scale' is available. 
#'   \item 'PEVD_GrpAve' : Groupd average PEVD, the optional arguments of 'scale' and 'diag' are available.
#'   \item 'PEVD_contrast' : Contrast PEVD, the optional argument of 'scale' is available.
#'   \item 'CD_IdAve' : Individual average CD.
#'   \item 'CD_GrpAve' : Group average CD, the optional argument of 'diag' is available.
#'   \item 'CD_contrast' : Contrast CD.
#'   \item 'r_IdAve' : Individual average r.
#'   \item 'r_GrpAve' : Group average r, the optional argument of 'diag' is available.
#'   \item 'r_contrast' : Contrast r. 
#'   \item 'VED0' : Variance of estimate of unit effects differences. The optional argument of 'scale' is available.
#'   \item 'VED1' : Variance of estimate of unit effects differences with the correction of unit effect. The optional argument of 'scale' is available.
#'   \item 'VED2' : Variance of estimate of unit effects differences with the correction of unit effect and additional fixed effects. 
#'     The additional argument of 'Uidx' is required and the optional argument of 'scale' is available. 
#'   \item 'CDVED0' : Coefficient of determination of VED, the optional argument of 'diag' is available.
#'   \item 'CDVED1' : Coefficient of determination of VED with the correction of unit effect. The optional argument of 'diag' is available.
#'   \item 'CDVED2' : Coefficient of determination of VED with the correction of unit effect and additional fixed effects. 
#'     The additional argument of 'Uidx' is required and the optional argument of 'diag' is available.
#'   \item 'CR0' : Connectedness rating. 
#'   \item 'CR1' : Connectedness rating with the correction of unit effect. 
#'   \item 'CR2' : Connectedness rating with the correction of unit effect and additional fixed effects. The additional argument of 'Uidx' is required. 
#' }
#' @param NumofMU The number of management unit to summarize connectedness. The available options include 'Pairwise' and 'Overall', where the prior calculates 
#'   the connectedness across all pairwise units, and the later averages all pairwise connectedness across units. 
#' @param Uidx An interger to indicate the last column of unit effects in \code{Xmatrix}. This Uidx is required for statistics VED2, CDVED2 and CR2. The default is NULL. 
#' @param scale Logical argument. Should \code{sigma2a} be used to scale \code{statistic} (e.g., PEVD_IdAve, PEVD_GrpAve, PEVD_contrast, VED0, VED1, and VED2) 
#'   to remove units effects? Default is TRUE.
#' @param diag Logical argument. Should diagonal elements of PEV matrix (e.g., PEVD_GrpAve, CD_GrpAve, and r_GrpAve) or K matrix (CDVED0, CDVED1, and CDVED2) be included? Default is TRUE.  
#' 
#' @return 
#' A value of overall connectedness measurements across units when \code{NumofMU} is set as 'Overall'.
#' A matrix of connectedness measurments with diagonal as NA when \code{NumofMU} is set as 'Pairwise'.
#' 
#' @author Haipeng Yu and Gota Morota 
#' 
#' Maintainer: Haipeng Yu \email{haipengyu@@vt.edu}
#'
#' @example man/examples/gca.R
#' 
#' @export 
gca <- function(Kmatrix, Xmatrix, sigma2a, sigma2e, MUScenario, statistic, NumofMU, Uidx = NULL, scale = TRUE, diag = TRUE) {
  if(is.factor(MUScenario) != TRUE) stop("Management unit is not factor!")
  Z <- diag(x = 1, nrow = nrow(Kmatrix), ncol = nrow(Kmatrix))
  diag(Kmatrix) <- diag(Kmatrix) + 0.00001
  Kinv <- chol2inv(chol(Kmatrix))
  lamda <- sigma2e/sigma2a
  X <- Xmatrix
  Mabs <- diag(nrow(X)) - Matprod(Matprod(X, chol2inv(chol(Matprod(t(X), X)))), t(X))
  CuuK = chol2inv(chol(Matprod(Matprod(t(Z), Mabs), Z) + Kinv * lamda))
  PEVK = CuuK * sigma2e
  Management_Unit <- unique(MUScenario)
  switch(statistic,
         PEVD_IdAve = pevd.idAve(PEVK = PEVK, Kmatrix = Kmatrix, Management_Unit = Management_Unit, 
                                 MUScenario = MUScenario, NumofMU = NumofMU, sigma2a = sigma2a, scale = scale),
         PEVD_GrpAve = pevd.grpAve(PEVK = PEVK, Kmatrix = Kmatrix, Management_Unit = Management_Unit, 
                                   MUScenario = MUScenario, NumofMU = NumofMU, sigma2a = sigma2a, scale = scale, diag = diag),
         PEVD_contrast = gc.contrast(statistic = statistic, Management_Unit = Management_Unit, MUScenario = MUScenario, 
                                     Kmatrix = Kmatrix, sigma2a = sigma2a, sigma2e = sigma2e, CuuK = CuuK, NumofMU = NumofMU, scale = scale),
         VED0 = ved.cor(statistic = statistic, NumofMU, Z = Z, Kinv = Kinv, lamda = lamda, X = X, Uidx = Uidx, sigma2a = sigma2a, 
                        sigma2e = sigma2e, Management_Unit = Management_Unit, scale = scale),
         VED1 = ved.cor(statistic = statistic, NumofMU, Z = Z, Kinv = Kinv, lamda = lamda, X = X, Uidx = Uidx, sigma2a = sigma2a, 
                        sigma2e = sigma2e, Management_Unit = Management_Unit, scale = scale),
         VED2 = ved.cor(statistic = statistic, NumofMU, Z = Z, Kinv = Kinv, lamda = lamda, X = X, Uidx = Uidx, sigma2a = sigma2a, 
                        sigma2e = sigma2e, Management_Unit = Management_Unit, scale = scale),
         CD_IdAve = cd.idAve(PEVK = PEVK, Kmatrix = Kmatrix, Management_Unit = Management_Unit,
                             MUScenario = MUScenario, NumofMU = NumofMU, sigma2a = sigma2a),
         CD_GrpAve = cd.grpAve(PEVK = PEVK, Kmatrix = Kmatrix, Management_Unit = Management_Unit,
                               MUScenario = MUScenario, NumofMU = NumofMU, sigma2a = sigma2a, diag = diag),
         CD_contrast = gc.contrast(statistic = statistic, Management_Unit = Management_Unit, MUScenario = MUScenario, Kmatrix = Kmatrix,
                                   sigma2a = sigma2a, sigma2e = sigma2e, CuuK = CuuK, NumofMU = NumofMU),
         CDVED0 = cd.cor(statistic = statistic, NumofMU, Z = Z, Kinv = Kinv, Kmatrix = Kmatrix, lamda = lamda,
                         X = X, Uidx = Uidx, sigma2e = sigma2e, sigma2a =sigma2a, Management_Unit = Management_Unit, MUScenario = MUScenario, diag = diag),
         CDVED1 = cd.cor(statistic = statistic, NumofMU, Z = Z, Kinv = Kinv, Kmatrix = Kmatrix, lamda = lamda,
                         X = X, Uidx = Uidx, sigma2e = sigma2e, sigma2a =sigma2a, Management_Unit = Management_Unit, MUScenario = MUScenario, diag = diag),
         CDVED2 = cd.cor(statistic = statistic, NumofMU, Z = Z, Kinv = Kinv, Kmatrix = Kmatrix, lamda = lamda,
                         X = X, Uidx = Uidx, sigma2e = sigma2e, sigma2a =sigma2a, Management_Unit = Management_Unit, MUScenario = MUScenario, diag = diag),
         r_IdAve = r.idAve(PEVK = PEVK, Management_Unit = Management_Unit, MUScenario = MUScenario, NumofMU = NumofMU),
         r_GrpAve = r.grpAve(PEVK = PEVK, Kmatrix = Kmatrix, Management_Unit = Management_Unit,
                             MUScenario = MUScenario, NumofMU = NumofMU, diag = diag),
         r_contrast = gc.contrast(statistic = statistic, Management_Unit = Management_Unit, MUScenario = MUScenario, Kmatrix = Kmatrix,
                                  sigma2a = sigma2a, sigma2e = sigma2e, CuuK = CuuK, NumofMU = NumofMU),
         CR0 = r.cor (statistic = statistic, NumofMU, Z = Z, Kinv = Kinv, lamda = lamda,
                      X = X, Uidx = Uidx, sigma2e = sigma2e, Management_Unit = Management_Unit),
         CR1 = r.cor (statistic = statistic, NumofMU, Z = Z, Kinv = Kinv, lamda = lamda,
                      X = X, Uidx = Uidx, sigma2e = sigma2e, Management_Unit = Management_Unit),
         CR2 = r.cor (statistic = statistic, NumofMU, Z = Z, Kinv = Kinv, lamda = lamda,
                      X = X, Uidx = Uidx, sigma2e = sigma2e, Management_Unit = Management_Unit)
  )
}


# Contrast for PEVD, CD & r
gc.contrast <- function(statistic = statistic, Management_Unit = Management_Unit, MUScenario = MUScenario, Kmatrix = Kmatrix,
                        sigma2a = sigma2a, sigma2e = sigma2e, CuuK = CuuK, NumofMU = NumofMU, scale = scale) {
  Contrast_Matrix <-  matrix(NA, ncol = length(Management_Unit), nrow = length(Management_Unit))
  for (i in 1 : (length(Management_Unit) - 1)) {
    for(j in (i + 1) : length(Management_Unit)) {
      xcontrast <- matrix(0, ncol = 1, nrow = ncol(Kmatrix))
      indexi <- which(MUScenario == Management_Unit[i])
      indexj <- which(MUScenario == Management_Unit[j])
      xcontrast[indexi] <- 1 / length(which(MUScenario == Management_Unit[i]))
      xcontrast[indexj] <- -1 / length(which(MUScenario == Management_Unit[j]))
      if (statistic == 'PEVD_contrast') {
        PEVK <- CuuK * sigma2e
        if(scale) Contrast_Matrix[i,j] <- Matprod(Matprod(t(xcontrast), PEVK), xcontrast) / sigma2a
        if(!scale) Contrast_Matrix[i,j] <- Matprod(Matprod(t(xcontrast), PEVK), xcontrast)  
      } else if (statistic == 'CD_contrast') {
        lamda <- sigma2e / sigma2a
        Contrast_Matrix[i,j] <- 
          1 - (lamda * Matprod(Matprod(t(xcontrast), CuuK), xcontrast) / Matprod(Matprod(t(xcontrast), Kmatrix), xcontrast))  
      } else if (statistic == 'r_contrast') {
        PEVK <- CuuK * sigma2e
        rijK <- cov2cor(PEVK)
        Contrast_Matrix[i,j] <- Matprod(Matprod(t(xcontrast), rijK), xcontrast)
      }
    }
  }
  if (NumofMU == 'Pairwise') {
    Contrast_Matrix[lower.tri(Contrast_Matrix)] <- t(Contrast_Matrix)[lower.tri(Contrast_Matrix)]
    colnames(Contrast_Matrix) <- rownames(Contrast_Matrix) <- Management_Unit
    return(Contrast_Matrix)
  } else if (NumofMU == 'Overall') {
    Contrast_Matrix.overall <- mean(Contrast_Matrix[upper.tri(Contrast_Matrix)])
    return(Contrast_Matrix.overall)
  }
}


# Individual Ave PEVD
pevd.idAve <- function(PEVK = PEVK, Kmatrix = Kmatrix, Management_Unit = Management_Unit, 
                       MUScenario = MUScenario, NumofMU = NumofMU, sigma2a = sigma2a, scale = scale) {
  PEVD.Pairwise.all <- matrix(NA, ncol = ncol(PEVK), nrow = nrow(PEVK))
  for (i in 1 : (nrow(Kmatrix) - 1)) {
    for (j in (i + 1) : nrow(Kmatrix)) {
      PEVD.Pairwise.all[i,j] <- PEVD.Pairwise.all[j,i] <- (PEVK[i,i] + PEVK[j,j] - 2 * PEVK[i,j])
    }
  }
  if(scale) PEVD.Pairwise.all <- PEVD.Pairwise.all / sigma2a
  if(!scale) PEVD.Pairwise.all <- PEVD.Pairwise.all 
  PEVD.Pairwise.unit <- matrix(NA, ncol = length(Management_Unit), 
                               nrow = length(Management_Unit))
  colnames(PEVD.Pairwise.unit) <- rownames(PEVD.Pairwise.unit) <- Management_Unit
  for (i in 1 : (length(Management_Unit) - 1)) {
    for(j in (i + 1) : length(Management_Unit)) {
      indexi <- which(MUScenario == Management_Unit[i])
      indexj <- which(MUScenario == Management_Unit[j])
      PEVD.Pairwise.unit[i, j] <- PEVD.Pairwise.unit[j, i] <- mean(PEVD.Pairwise.all[indexi, indexj])
    }
  }
  if (NumofMU == 'Pairwise') {
    return(PEVD.Pairwise.unit)
  } else if (NumofMU == 'Overall') {
    PEVD.IdAve.overall <- mean(PEVD.Pairwise.unit[upper.tri(PEVD.Pairwise.unit)])
    return(PEVD.IdAve.overall)
  }
}


# small function to calculate difference with diag argument
diff <- function(x, y, z, diag = TRUE) {
  if(diag) d <- mean(x) + mean(y) - 2 * z
  if(!diag) d <- mean(x[upper.tri(x)]) + mean(y[upper.tri(y)]) - 2 * z
  return(d)
}


# Group Ave PEVD
pevd.grpAve <- function(PEVK = PEVK, Kmatrix = Kmatrix, Management_Unit = Management_Unit, 
                        MUScenario = MUScenario, NumofMU = NumofMU, sigma2a = sigma2a, scale = scale, diag = diag) {
  PEVD.Pairwise.unit <- matrix(NA, ncol = length(Management_Unit), nrow = length(Management_Unit))
  colnames(PEVD.Pairwise.unit) <- rownames(PEVD.Pairwise.unit) <- Management_Unit
  for (i in 1 : (length(Management_Unit) - 1)) {
    for(j in (i + 1) : length(Management_Unit)) {
      indexi <- which(MUScenario == Management_Unit[i])
      indexj <- which(MUScenario == Management_Unit[j])
      PEV.ii <- PEVK[indexi, indexi]
      PEV.jj <- PEVK[indexj, indexj]
      PEC.ij.mean <- mean(PEVK[indexi, indexj])
      PEVD.Pairwise.unit[i, j] <- PEVD.Pairwise.unit[j, i] <- diff(PEV.ii, PEV.jj, PEC.ij.mean, diag = diag)
    }
  }
  if (scale) PEVD.Pairwise.unit <- PEVD.Pairwise.unit / sigma2a
  if (!scale) PEVD.Pairwise.unit <- PEVD.Pairwise.unit
  if (NumofMU == 'Pairwise') {
    return(PEVD.Pairwise.unit)
  } else if (NumofMU == 'Overall') {
    PEVD.GrpAve.overall <- mean(PEVD.Pairwise.unit[upper.tri(PEVD.Pairwise.unit)])
    return(PEVD.GrpAve.overall)
  }
}


# VED0; VED1; VED2
ved.cor <- function(statistic = statistic, NumofMU, Z = Z, Kinv = Kinv, lamda = lamda,
                    X = X, Uidx = Uidx, sigma2a = sigma2a, sigma2e = sigma2e, Management_Unit = Management_Unit, scale = scale) {
  C22_inv <- chol2inv(chol(Matprod(t(Z), Z) + Kinv * lamda))
  var_Bhat <- chol2inv(chol(Matprod(t(X), X) - (Matprod(Matprod(Matprod(t(X), Z), C22_inv), Matprod(t(Z), X))))) * sigma2e
  Summary.Matrix <- matrix(NA, ncol = length(Management_Unit), nrow = length(Management_Unit))
  colnames(Summary.Matrix) <- rownames(Summary.Matrix) <- Management_Unit
  if(statistic == 'VED0') {
    PEV_mean <- var_Bhat
  } else if (statistic == 'VED1') {
    cor_factor <- chol2inv(chol(Matprod(t(X), X))) * sigma2e
    PEV_mean <- var_Bhat - cor_factor
  } else if (statistic == 'VED2') {
    if(is.null(Uidx)) stop("Please specify the length of unit effect in X matrix with argument of Uidx")
    m <- dim(X)[2] # length of total fixed effects
    m1 <- Uidx # lengt of MU
    if(m1 >= m) stop("This statistic requires Uidx to be smaller than the number of columns of X matrix")
    if(m1 <= 1) stop("The number of unit must be larger than 1")
    X1 <- X[, 1 : m1] # design matrix for MU
    X2 <- X[, -c(1 : m1)]
    if(is.vector(X2)) {
      X2 <- as.matrix(X2)
      tX1X1_inv <- chol2inv(chol(Matprod(t(X1), X1)))
      tX1X2 <- Matprod(t(X1), X2)
      tX2X1 <- Matprod(t(X2), X1)
      PEV_mean <- Matprod(Matprod(Matprod(tX1X1_inv, tX1X2) * var_Bhat[(m1 + 1) : m, (m1 + 1) : m], tX2X1), tX1X1_inv) + 
        Matprod(Matprod(tX1X1_inv, tX1X2), matrix(var_Bhat[(m1 + 1) : m, 1 : m1], ncol = m1)) + 
        Matprod(Matprod(matrix(var_Bhat[1 : m1, (m1 + 1) : m], nrow = m1), tX2X1), tX1X1_inv) + 
        var_Bhat[1 : m1, 1 : m1] - sigma2e * tX1X1_inv
    } else {
      tX1X1_inv <- chol2inv(chol(Matprod(t(X1), X1)))
      tX1X2 <- Matprod(t(X1), X2)
      tX2X1 <- Matprod(t(X2), X1)
      PEV_mean <- Matprod(Matprod(Matprod(Matprod(tX1X1_inv, tX1X2), var_Bhat[(m1 + 1) : m, (m1 + 1) : m]), tX2X1), tX1X1_inv) + 
        Matprod(Matprod(tX1X1_inv, tX1X2), var_Bhat[(m1 + 1) : m, 1 : m1]) + 
        Matprod(Matprod(var_Bhat[1 : m1, (m1 + 1) : m], tX2X1), tX1X1_inv) + 
        var_Bhat[1 : m1, 1 : m1] - sigma2e * tX1X1_inv
    }
  }
  for (i in 1 : ncol(Summary.Matrix) - 1) {
    for (j in (i + 1) : ncol(Summary.Matrix)) {
      Summary.Matrix[i, j] <- Summary.Matrix[j, i] <- (PEV_mean[i, i] + PEV_mean[j, j] - 2 * PEV_mean[i, j])
    }
  }
  if(scale) Summary.Matrix <- Summary.Matrix / sigma2a
  if(!scale) Summary.Matrix <- Summary.Matrix
  if (NumofMU == 'Pairwise') {
    return(Summary.Matrix)
  } else if (NumofMU == 'Overall') {
    Summary.Matrix.overall <- mean(Summary.Matrix[upper.tri(Summary.Matrix)])
    return(Summary.Matrix.overall)
  }
}


# Individual Average CD
cd.idAve <- function(PEVK = PEVK, Kmatrix = Kmatrix, Management_Unit = Management_Unit, 
                     MUScenario = MUScenario, NumofMU = NumofMU, sigma2a = sigma2a) {
  PEVD.Pairwise.all <- K.diff.Pairwise.all <- matrix(NA, ncol = ncol(PEVK), nrow = nrow(PEVK))
  for (i in 1 : (nrow(Kmatrix) - 1)) {
    for (j in (i + 1) : nrow(Kmatrix)) {
      PEVD.Pairwise.all[i,j] <- PEVD.Pairwise.all[j,i] <- (PEVK[i,i] + PEVK[j,j] - 2 * PEVK[i,j])
      K.diff.Pairwise.all[i,j] <- K.diff.Pairwise.all[j,i] <- (Kmatrix[i,i] + Kmatrix[j,j] - 2 * Kmatrix[i,j])
    }
  }
  CD.Pairwise.unit <- matrix(NA, ncol = length(Management_Unit), 
                             nrow = length(Management_Unit))
  colnames(CD.Pairwise.unit) <- rownames(CD.Pairwise.unit) <- Management_Unit
  for (i in 1 : (length(Management_Unit) - 1)) {
    for(j in (i + 1) : length(Management_Unit)) {
      indexi <- which(MUScenario == Management_Unit[i])
      indexj <- which(MUScenario == Management_Unit[j])
      CD.Pairwise.unit[i, j] <- CD.Pairwise.unit[j, i] <- 
        1- (mean(PEVD.Pairwise.all[indexi, indexj]) / (mean(K.diff.Pairwise.all[indexi, indexj]) * sigma2a))
    }
  }
  if (NumofMU == 'Pairwise') {
    return(CD.Pairwise.unit)
  } else if (NumofMU == 'Overall') {
    CD.IdAve.overall <- mean(CD.Pairwise.unit[upper.tri(CD.Pairwise.unit)])
    return(CD.IdAve.overall)
  }
}


# Group Average CD
cd.grpAve <- function(PEVK = PEVK, Kmatrix = Kmatrix, Management_Unit = Management_Unit, 
                      MUScenario = MUScenario, NumofMU = NumofMU, sigma2a = sigma2a, diag = diag) {
  CD.Pairwise.unit <- matrix(NA, ncol = length(Management_Unit), nrow = length(Management_Unit))
  colnames(CD.Pairwise.unit) <- rownames(CD.Pairwise.unit) <- Management_Unit
  for (i in 1 : (length(Management_Unit) - 1)) {
    for (j in (i + 1) : length(Management_Unit)) {
      indexi <- which(MUScenario == Management_Unit[i])
      indexj <- which(MUScenario == Management_Unit[j])
      PEV.ii <- PEVK[indexi, indexi]
      K.ii <- Kmatrix[indexi, indexi]
      PEV.jj <- PEVK[indexj, indexj]
      K.jj <- Kmatrix[indexj, indexj]
      PEC.ij.mean <- mean(PEVK[indexi, indexj])
      K.ij.mean <- mean(Kmatrix[indexi, indexj])
      PEVD.ij <- diff(PEV.ii, PEV.jj, PEC.ij.mean, diag = diag)
      K.diff.ij <- diff(K.ii, K.jj, K.ij.mean, diag = diag)
      CD.Pairwise.unit[i, j] <- CD.Pairwise.unit[j, i] <- 1 - (PEVD.ij / (K.diff.ij * sigma2a))
    }
  }
  if (NumofMU == 'Pairwise') {
    return(CD.Pairwise.unit)
  } else if (NumofMU == 'Overall') {
    CD.GrpAve.overall <- mean(CD.Pairwise.unit[upper.tri(CD.Pairwise.unit)])
    return(CD.GrpAve.overall)
  }
}


# CDVED0, CDVED1, CDVED2
cd.cor <- function(statistic = statistic, NumofMU, Z = Z, Kinv = Kinv, Kmatrix = Kmatrix, lamda = lamda,
                   X = X, Uidx = Uidx, sigma2e = sigma2e, sigma2a =sigma2a, Management_Unit = Management_Unit,
                   MUScenario = MUScenario, diag = diag) {
  C22_inv <- chol2inv(chol(Matprod(t(Z), Z) + Kinv * lamda))
  var_Bhat <- chol2inv(chol(Matprod(t(X), X) - (Matprod(Matprod(Matprod(t(X), Z), C22_inv), Matprod(t(Z), X))))) * sigma2e
  Summary.Matrix <- K.Pairwise.unit.diff <- matrix(NA, ncol = length(Management_Unit), nrow = length(Management_Unit))
  colnames(Summary.Matrix) <- rownames(Summary.Matrix) <- Management_Unit
  for (i in 1 : (length(Management_Unit) - 1)) {
    for(j in (i + 1) : length(Management_Unit)) {
      indexi <- which(MUScenario == Management_Unit[i])
      indexj <- which(MUScenario == Management_Unit[j])
      K.ii <- Kmatrix[indexi, indexi]
      K.jj <- Kmatrix[indexj, indexj]
      K.ij <- mean(Kmatrix[indexi, indexj])
      K.Pairwise.unit.diff[i,j] <- K.Pairwise.unit.diff[j,i] <- diff(K.ii, K.jj, K.ij, diag = diag) * sigma2a
    }
  }
  if(statistic == 'CDVED0') {
    PEV_mean <- var_Bhat
  } else if (statistic == 'CDVED1') {
    cor_factor <- chol2inv(chol(Matprod(t(X), X))) * sigma2e
    PEV_mean <- var_Bhat - cor_factor
  } else if (statistic == 'CDVED2') {
    if(is.null(Uidx)) stop("Please specify the length of unit effect in X matrix with argument of Uidx")
    m <- dim(X)[2] # length of total fixed effects
    m1 <- Uidx # lengt of MU
    if(m1 >= m) stop("This statistic requires Uidx to be smaller than the number of columns of X matrix")
    if(m1 <= 1) stop("The number of unit must be larger than 1")
    X1 <- X[, 1 : m1] # design matrix for MU
    X2 <- X[, -c(1 : m1)]
    if(is.vector(X2)) {
      X2 <- as.matrix(X2)
      tX1X1_inv <- chol2inv(chol(Matprod(t(X1), X1)))
      tX1X2 <- Matprod(t(X1), X2)
      tX2X1 <- Matprod(t(X2), X1)
      PEV_mean <- Matprod(Matprod(Matprod(tX1X1_inv, tX1X2) * var_Bhat[(m1 + 1) : m, (m1 + 1) : m], tX2X1), tX1X1_inv) + 
        Matprod(Matprod(tX1X1_inv, tX1X2), matrix(var_Bhat[(m1 + 1) : m, 1 : m1], ncol = m1)) + 
        Matprod(Matprod(matrix(var_Bhat[1 : m1, (m1 + 1) : m], nrow = m1), tX2X1), tX1X1_inv) + 
        var_Bhat[1 : m1, 1 : m1] - sigma2e * tX1X1_inv
    } else {
      tX1X1_inv <- chol2inv(chol(Matprod(t(X1), X1)))
      tX1X2 <- Matprod(t(X1), X2)
      tX2X1 <- Matprod(t(X2), X1)
      PEV_mean <- Matprod(Matprod(Matprod(Matprod(tX1X1_inv, tX1X2), var_Bhat[(m1 + 1) : m, (m1 + 1) : m]), tX2X1), tX1X1_inv) + 
        Matprod(Matprod(tX1X1_inv, tX1X2), var_Bhat[(m1 + 1) : m, 1 : m1]) + 
        Matprod(Matprod(var_Bhat[1 : m1, (m1 + 1) : m], tX2X1), tX1X1_inv) + 
        var_Bhat[1 : m1, 1 : m1] - sigma2e * tX1X1_inv
    }
  }
  for (i in 1 : ncol(Summary.Matrix) - 1) {
    for (j in (i + 1) : ncol(Summary.Matrix)) {
      Summary.Matrix[i, j] <- Summary.Matrix[j, i] <- 1 - (
        (PEV_mean[i, i] + PEV_mean[j, j] - 2 * PEV_mean[i, j]) / K.Pairwise.unit.diff[i, j])
    }
  } 
  if (NumofMU == 'Pairwise') {
    return(Summary.Matrix)
  } else if (NumofMU == 'Overall') {
    Summary.Matrix.overall <- mean(Summary.Matrix[upper.tri(Summary.Matrix)])
    return(Summary.Matrix.overall)
  }
}


# Individual r
r.idAve <- function(PEVK = PEVK, Management_Unit = Management_Unit, MUScenario = MUScenario, NumofMU = NumofMU) {
  rijK <- cov2cor(PEVK)
  rij.Pairwise.unit <- matrix(NA, ncol=length(Management_Unit), nrow=length(Management_Unit))
  colnames(rij.Pairwise.unit) <- rownames(rij.Pairwise.unit) <- Management_Unit
  for (i in 1 : (length(Management_Unit) - 1)) {
    for(j in (i + 1) : length(Management_Unit)) {
      indexi <- which(MUScenario == Management_Unit[i])
      indexj <- which(MUScenario == Management_Unit[j])
      rij.Aacross <- rijK[indexi, indexj]
      rij.Pairwise.unit[i, j] <- rij.Pairwise.unit[j, i] <- mean(rij.Aacross, na.rm = TRUE)
    }
  }
  if (NumofMU == 'Pairwise'){
    return(rij.Pairwise.unit)
  } else if (NumofMU == 'Overall'){
    rij.Pairwise.unit.overall <- mean(rij.Pairwise.unit[upper.tri(rij.Pairwise.unit)])
    return(rij.Pairwise.unit.overall)
  }
}


# small function to calculate cor with diag argument
corr <- function(x, y, z, diag = TRUE){
  if(diag) c <- z / (sqrt(mean(x) * mean(y)))
  if(!diag) c <- z / (sqrt(mean(x[upper.tri(x)]) * mean(y[upper.tri(y)])))
  return(c)
}


# Group Average r
r.grpAve <- function(PEVK = PEVK, Kmatrix = Kmatrix, Management_Unit = Management_Unit, 
                     MUScenario = MUScenario, NumofMU = NumofMU, diag = diag) {
  r.Pairwise.unit <- matrix(NA, ncol = length(Management_Unit), nrow = length(Management_Unit))
  colnames(r.Pairwise.unit) <- rownames(r.Pairwise.unit) <- Management_Unit
  for (i in 1 : (length(Management_Unit) - 1)) {
    for(j in (i + 1) : length(Management_Unit)) {
      indexi <- which(MUScenario == Management_Unit[i])
      indexj <- which(MUScenario == Management_Unit[j])
      PEV.ii <- PEVK[indexi, indexi]
      PEV.jj <- PEVK[indexj, indexj]
      PEC.ij.mean <- mean(PEVK[indexi, indexj])
      r.Pairwise.unit[i, j] <- r.Pairwise.unit[j, i] <- corr(PEV.ii, PEV.jj, PEC.ij.mean, diag = diag)
    }
  }
  if (NumofMU == 'Pairwise') {
    return(r.Pairwise.unit)
  } else if (NumofMU == 'Overall') {
    r.GrpAve.overall <- mean(r.Pairwise.unit[upper.tri(r.Pairwise.unit)])
    return(r.GrpAve.overall)
  }
}


# CR0, CR1, CR2
r.cor <- function(statistic = statistic, NumofMU, Z = Z, Kinv = Kinv, lamda = lamda,
                  X = X, Uidx = Uidx, sigma2e = sigma2e, Management_Unit = Management_Unit) {
  C22_inv <- chol2inv(chol(Matprod(t(Z), Z) + Kinv * lamda))
  var_Bhat <- chol2inv(chol(Matprod(t(X), X) - (Matprod(Matprod(Matprod(t(X), Z), C22_inv), Matprod(t(Z), X))))) * sigma2e
  Summary.Matrix <- matrix(NA, ncol = length(Management_Unit), nrow = length(Management_Unit))
  colnames(Summary.Matrix) <- rownames(Summary.Matrix) <- Management_Unit
  if(statistic == 'CR0') {
    PEV_mean <- var_Bhat
  } else if (statistic == 'CR1') {
    cor_factor <- chol2inv(chol(Matprod(t(X), X))) * sigma2e
    PEV_mean <- var_Bhat - cor_factor
  } else if (statistic == 'CR2') {
    if(is.null(Uidx)) stop("Please specify the length of unit effect in X matrix with argument of Uidx")
    m <- dim(X)[2] # length of total fixed effects
    m1 <- Uidx # lengt of MU
    if(m1 >= m) stop("This statistic requires Uidx to be smaller than the number of columns of X matrix")
    if(m1 <= 1) stop("The number of unit must be larger than 1")
    X1 <- X[, 1 : m1] # design matrix for MU
    X2 <- X[, -c(1 : m1)]
    if(is.vector(X2)) {
      X2 <- as.matrix(X2)
      tX1X1_inv <- chol2inv(chol(Matprod(t(X1), X1)))
      tX1X2 <- Matprod(t(X1), X2)
      tX2X1 <- Matprod(t(X2), X1)
      PEV_mean <- Matprod(Matprod(Matprod(tX1X1_inv, tX1X2) * var_Bhat[(m1 + 1) : m, (m1 + 1) : m], tX2X1), tX1X1_inv) + 
        Matprod(Matprod(tX1X1_inv, tX1X2), matrix(var_Bhat[(m1 + 1) : m, 1 : m1], ncol = m1)) + 
        Matprod(Matprod(matrix(var_Bhat[1 : m1, (m1 + 1) : m], nrow = m1), tX2X1), tX1X1_inv) + 
        var_Bhat[1 : m1, 1 : m1] - sigma2e * tX1X1_inv
    } else {
      tX1X1_inv <- chol2inv(chol(Matprod(t(X1), X1)))
      tX1X2 <- Matprod(t(X1), X2)
      tX2X1 <- Matprod(t(X2), X1)
      PEV_mean <- Matprod(Matprod(Matprod(Matprod(tX1X1_inv, tX1X2), var_Bhat[(m1 + 1) : m, (m1 + 1) : m]), tX2X1), tX1X1_inv) + 
        Matprod(Matprod(tX1X1_inv, tX1X2), var_Bhat[(m1 + 1) : m, 1 : m1]) + 
        Matprod(Matprod(var_Bhat[1 : m1, (m1 + 1) : m], tX2X1), tX1X1_inv) + 
        var_Bhat[1 : m1, 1 : m1] - sigma2e * tX1X1_inv
    }
  }
  for (i in 1 : ncol(Summary.Matrix) - 1) {
    for (j in (i + 1) : ncol(Summary.Matrix)) {
      Summary.Matrix[i, j] <- Summary.Matrix[j, i] <- PEV_mean[i, j] / sqrt(PEV_mean[i, i] * PEV_mean[j, j])
    }
  }
  if (NumofMU == 'Pairwise') {
    return(Summary.Matrix)
  } else if (NumofMU == 'Overall') {
    Summary.Matrix.overall <- mean(Summary.Matrix[upper.tri(Summary.Matrix)])
    return(Summary.Matrix.overall)
  }
}


