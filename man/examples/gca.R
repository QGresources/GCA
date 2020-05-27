# Load cattle data
data(GCcattle)

# Compute genomic relationship matrix
G <- computeG(cattle.W, maf = 0.05, impute = 'mean', method = 'G1')

# The heritability of simulated phenotype was set to 0.6 with additive genetic variace (Vu) = 0.6 and residual variance (Ve) = 0.4
var <- list(Vu = 0.6, Ve = 0.4) 

# Design matrix of fixed effects
## unit effect
X1 <- model.matrix(~ -1 + factor(cattle.pheno$Unit))
## unit effect and sex effect
X2 <- model.matrix(~ -1 + factor(cattle.pheno$Unit) + factor(cattle.pheno$Sex))

# Calculate CD_IdAve
CD_IdAve <- gca(Kmatrix = G, Xmatrix = X1, sigma2a = var$Vu, sigma2e = var$Ve, 
                MUScenario = as.factor(cattle.pheno$Unit), statistic = 'CD_IdAve', 
                NumofMU = 'Overall')

# Calculate CDVED1
CDVED1 <- gca(Kmatrix = G, Xmatrix = X1, sigma2a = var$Vu, sigma2e = var$Ve, 
              MUScenario = as.factor(cattle.pheno$Unit), statistic = 'CDVED1',
              NumofMU = 'Pairwise', diag = TRUE)

# Calculate CDVED2
CDVED2 <- gca(Kmatrix = G, Xmatrix = X2, sigma2a = var$Vu, sigma2e = var$Ve,
              MUScenario = as.factor(cattle.pheno$Unit), statistic = 'CDVED2', 
              NumofMU = 'Pairwise', Uidx = 8, diag = TRUE)



