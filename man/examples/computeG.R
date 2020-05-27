# Load cattle data
data(GCcattle)

# Marker information
str(cattle.W)

# Compute genomic relationship matrix
G <- computeG(cattle.W, maf = 0.05, impute = 'rbinom', method = 'G1')
