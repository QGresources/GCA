# Load GCA package
library(GCA)

# Import cattle dataset
data(GCcattle)

# Pedigree information
head(cattle.pheno)

# Construct numerator relationship matrix
NRM <- computeA(Progeny = cattle.pheno$Progeny, Sire = cattle.pheno$Sire, 
		Dam = cattle.pheno$Dam)
