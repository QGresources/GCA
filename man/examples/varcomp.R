# Load cattle data
data(GCcattle)

# Phenotype, fixed covariates and marker information
str(cattle.pheno)
str(cattle.W)

# Compute genomic relationship matrix
G <- computeG(cattle.W, maf = 0.05, impute = 'rbinom', method = 'G1')

# Eigendecomposition of genomic relationship matrix
EVD <- eigen(G)

# Estimate variance component
var <- varcomp(y = cattle.pheno$Phenotype, Evector = EVD$vectors, Evalue = EVD$values) 

# Retrieve additive genetic variance and residual variance
var$Vu
var$Ve
