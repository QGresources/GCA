## ----global_options, include=FALSE--------------------------------------------
knitr::opts_chunk$set(echo=TRUE, warning=FALSE, message=FALSE)

## ----echo = FALSE, fig.cap = 'Figure 1: An overview of connectedness statistics implemented in the GCA R package.'----
library(DiagrammeR)
grViz("digraph flowchart {
      # node definitions with substituted label text
      # graph [overlap = true, fontsize = 12]
      node [fontname = Helvetica, shape = box, fontsize = 12,
            color = white, fillcolor = Tomato2, style = filled]
      tab1;
      node [fontname = Helvetica, shape = box, fontsize = 11, fixedsize = true, 
            width = 0.6, height = 0.3, color = white, fillcolor = SeaGreen2, style = filled]
      tab4;tab5;tab6;
      node [fontname = Helvetica, shape = box, fontsize = 11, fixedsize = true, 
            width = 0.6, height = 0.3, color = white, fillcolor = deepskyblue2, style = filled]
      tab7;tab8;tab9;
      node [fontname = Helvetica, shape = circle, fixedsize = true, fontsize = 8, 
            width = 0.5, height = 0.3, color = white, fillcolor = SeaGreen1, style = filled]
      tab10;tab11;tab12;
      node [fontname = Helvetica, shape = circle, fixedsize = true, fontsize = 8, 
            width = 0.5, height = 0.3, color = white, fillcolor = deepskyblue1, style = filled]
      tab13;tab14;tab15;
      node [fontname = Helvetica, shape = oval, fixedsize = false, 
            fontsize = 12, color = white, fillcolor = SeaGreen3, style = filled]
      tab2;
      node [fontname = Helvetica, shape = oval, fixedsize = false, 
            fontsize = 12, color = white, fillcolor = deepskyblue3, style = filled]
      tab3;
      tab1 [label = '@@1']
      tab2 [label = '@@2']
      tab3 [label = '@@3']
      tab4 [label = '@@4']
      tab5 [label = '@@5']
      tab6 [label = '@@6']
      tab7 [label = '@@7']
      tab8 [label = '@@8']
      tab9 [label = '@@9']
      tab10 [label = '@@10']
      tab11 [label = '@@11']
      tab12 [label = '@@12']
      tab13 [label = '@@13']
      tab14 [label = '@@14']
      tab15 [label = '@@15']
      
      # edge definitions with the node IDs
      tab1 -> {tab2 tab3}[arrowsize = 0.3];
      tab2 -> {tab4 tab5 tab6}[arrowsize = 0.3];
      tab3 -> {tab7 tab8 tab9}[arrowsize = 0.3];
      tab4 -> {tab10 tab11 tab12}[arrowsize = 0.2];
      tab5 -> {tab10 tab11 tab12}[arrowsize = 0.2];
      tab6 -> {tab10 tab11 tab12}[arrowsize = 0.2];
      tab7 -> {tab13 tab14 tab15}[arrowsize = 0.2];
      tab8 -> {tab13 tab14 tab15}[arrowsize = 0.2];
      tab9 -> {tab13 tab14 tab15}[arrowsize = 0.2];
      }
      
      [1]: 'Connectedness statistics'
      [2]: 'PEV core function'
      [3]: 'VE core function'
      [4]: 'PEVD'
      [5]: 'CD'
      [6]: 'r'
      [7]: 'VED'
      [8]: 'CDVED'
      [9]: 'CR'
      [10]: 'IdAve'
      [11]: 'GrpAve'
      [12]: 'Contrast'
      [13]: '0'
      [14]: '1'
      [15]: '2'
      ")

## ---- results='hide', eval = FALSE--------------------------------------------
#  install.packages("devtools")
#  library(devtools)
#  install_github('HaipengU/GCA')

## ----message = FALSE----------------------------------------------------------
library(BGLR)
library(corrplot)
library(GCA)
library(ggplot2)
library(gridExtra)

## ----message = FALSE----------------------------------------------------------
data(GCcattle)
dim(cattle.pheno) # phenotypes + fixed effects
dim(cattle.W) # marker matrix

## -----------------------------------------------------------------------------
sigma2a <- 0.6 # additive genetic variance
sigma2e <- 0.4 # residual variance

## -----------------------------------------------------------------------------
X2 <- model.matrix(~ -1 + factor(cattle.pheno$Unit) + factor(cattle.pheno$Sex)) # incidence matrix of unit effect and sex
G <- computeG(cattle.W, maf = 0.05) # genomic relationship matrix; markers with minor allele frequency (maf) 
                                    # less than 0.05 are removed


## ---- eval=FALSE--------------------------------------------------------------
#  gca(Kmatrix, Xmatrix, sigma2a, sigma2e, MUScenario, statistic, NumofMU,
#      Uidx = NULL, scale = TRUE, diag = TRUE)

## -----------------------------------------------------------------------------
PEVD_GrpAve <- gca(Kmatrix = G, Xmatrix = X2, sigma2a = sigma2a, sigma2e = sigma2e,
                  MUScenario = as.factor(cattle.pheno$Unit), statistic = 'PEVD_GrpAve',
                  NumofMU = 'Pairwise')
# remove NAs in diagnol to make plot
diag(PEVD_GrpAve) <- 0 
corrplot(PEVD_GrpAve, is.corr = FALSE, method ="circle", type = "upper",
         diag = F, number.cex = 7 / ncol(PEVD_GrpAve), col = cm.colors(10), cl.lim = c(0, 0.08),
         number.digits = 4, mar = c(0,1,1,1), addCoef.col = "black",
         tl.col = "black", tl.srt = 90, tl.cex = 1.2)

## -----------------------------------------------------------------------------
PEVD_GrpAve_Overall <- gca(Kmatrix = G, Xmatrix = X2, sigma2a = sigma2a,
                           sigma2e = sigma2e, MUScenario = as.factor(cattle.pheno$Unit),
                           statistic = 'PEVD_GrpAve', NumofMU = 'Overall')
PEVD_GrpAve_Overall

## -----------------------------------------------------------------------------
CD_GrpAve <- gca(Kmatrix = G, Xmatrix = X2, sigma2a = sigma2a, sigma2e = sigma2e,
                MUScenario = as.factor(cattle.pheno$Unit), statistic = 'CD_GrpAve',
                NumofMU = 'Pairwise')
# replace NAs in diagnol to make plot
diag(CD_GrpAve) <- min(CD_GrpAve, na.rm = T) 
corrplot(CD_GrpAve, is.corr = FALSE, method ="circle", type = "upper",
         diag = F, number.cex = 7 / ncol(CD_GrpAve), col = cm.colors(10), 
         number.digits = 4, mar = c(0,1,1,1), addCoef.col = "black",
         tl.col = "black", tl.srt = 90, tl.cex = 1.2, cl.lim = c(0.5, 0.8))


## -----------------------------------------------------------------------------
VED0  <- gca(Kmatrix = G, Xmatrix = X2, sigma2a = sigma2a, sigma2e = sigma2e,
             MUScenario = as.factor(cattle.pheno$Unit), statistic = 'VED0',
             NumofMU = 'Pairwise')
VED1  <- gca(Kmatrix = G, Xmatrix = X2, sigma2a = sigma2a, sigma2e = sigma2e,
             MUScenario = as.factor(cattle.pheno$Unit), statistic = 'VED1',
             NumofMU = 'Pairwise')
VED2  <- gca(Kmatrix = G, Xmatrix = X2, sigma2a = sigma2a, sigma2e = sigma2e,
               MUScenario = as.factor(cattle.pheno$Unit), statistic = 'VED2',
               NumofMU = 'Pairwise', Uidx = 8)
# replace NAs in diagnol to make plot
diag(VED0) <- diag(VED1) <- diag(VED2) <- floor(min(VED0, na.rm = T)) 
corrplot(VED0, is.corr =FALSE, method="circle", type = "upper",
         diag = F, number.cex = 7 / ncol(VED0), col = cm.colors(10),
         number.digits = 4, addCoef.col = "black", tl.col = "black", 
         tl.srt = 90, tl.cex = 1.2, cl.lim = c(0, 0.08))
corrplot(VED1, is.corr =FALSE, method="circle", type = "upper",
         diag = F, number.cex = 7 / ncol(VED1), col = cm.colors(10),
         number.digits = 4, addCoef.col = "black", tl.col = "black", 
         tl.srt = 90, tl.cex = 1.2, cl.lim = c(0, 0.08))
corrplot(VED2, is.corr =FALSE, method="circle", type = "upper",
         diag = F, number.cex= 7 / ncol(VED2), col = cm.colors(10),
         number.digits = 4, addCoef.col = "black", tl.col = "black", 
         tl.srt = 90, tl.cex = 1.2, cl.lim = c(0, 0.08))


## -----------------------------------------------------------------------------
CDVED0  <- gca(Kmatrix = G, Xmatrix = X2, sigma2a = sigma2a, sigma2e = sigma2e,
             MUScenario = as.factor(cattle.pheno$Unit), statistic = 'CDVED0',
             NumofMU = 'Pairwise')
CDVED1  <- gca(Kmatrix = G, Xmatrix = X2, sigma2a = sigma2a, sigma2e = sigma2e,
             MUScenario = as.factor(cattle.pheno$Unit), statistic = 'CDVED1',
             NumofMU = 'Pairwise')
CDVED2  <- gca(Kmatrix = G, Xmatrix = X2, sigma2a = sigma2a, sigma2e = sigma2e,
             MUScenario = as.factor(cattle.pheno$Unit), statistic = 'CDVED2',
             NumofMU = 'Pairwise', Uidx = 8)
# replace NAs in diagnol to make plot
diag(CDVED0) <- diag(CDVED1) <- diag(CDVED2) <- 0.5 
corrplot(CDVED0, is.corr = FALSE, method ="circle", type = "upper",
         diag = F, number.cex = 7 / ncol(CDVED0), col = cm.colors(10),
         number.digits = 4, addCoef.col = "black", tl.col = "black", 
         tl.srt = 90, tl.cex = 1.2, cl.lim = c(0.5, 0.8))
corrplot(CDVED1, is.corr = FALSE, method ="circle", type = "upper",
         diag = F, number.cex = 7 / ncol(CDVED1), col = cm.colors(10),
         number.digits = 4, addCoef.col = "black", tl.col = "black", 
         tl.srt = 90, tl.cex = 1.2, cl.lim = c(0.5, 0.8))
corrplot(CDVED2, is.corr = FALSE, method ="circle", type = "upper",
         diag = F, number.cex = 7 / ncol(CDVED2), col = cm.colors(10),
         number.digits = 4, addCoef.col = "black", tl.col = "black", 
         tl.srt = 90, tl.cex = 1.2, cl.lim = c(0.5, 0.8))

## ----fig.cap = 'Figure 2: Correlation between PEVD_GrpAve and VED0, VED1, and VED2.'----
df_PEVD_VED <- data.frame(PEVD_GrpAve = PEVD_GrpAve[upper.tri(PEVD_GrpAve)],
                          VED0 = VED0[upper.tri(VED0)], 
                          VED1 = VED1[upper.tri(VED1)],
                          VED2 = VED2[upper.tri(VED2)])

corrplot(cor(df_PEVD_VED), method="circle", type = "upper",
         diag = F, number.cex= 4 / ncol(df_PEVD_VED), col = cm.colors(10), 
         number.digits = 6, addCoef.col = "black", tl.col = "black", 
         tl.srt = 90, tl.cex = 1.2, cl.lim = c(0.8, 1))

## ----fig.cap = 'Figure 3: Correlation between CD_GrpAve and CDVED0, CDVED1, and CDVED2.'----
df_CD_CDVED <- data.frame(CD_GrpAve = CD_GrpAve[upper.tri(CD_GrpAve)],
                          CDVED0 = CDVED0[upper.tri(CDVED0)], 
                          CDVED1 = CDVED1[upper.tri(CDVED1)],
                          CDVED2 = CDVED2[upper.tri(CDVED2)])

corrplot(cor(df_CD_CDVED), method="circle", type = "upper",
         diag = F, number.cex= 4 / ncol(df_CD_CDVED), col = cm.colors(10), 
         number.digits = 6, addCoef.col = "black", tl.col = "black", 
         tl.srt = 90, tl.cex = 1.2, cl.lim = c(0.8, 1))

## ----fig.cap = 'Figure 4: Correlation between r_GrpAve and CR0, CR1, and CR2.'----
r_GrpAve <- gca(Kmatrix = G, Xmatrix = X2, sigma2a = sigma2a, sigma2e = sigma2e,
                  MUScenario = as.factor(cattle.pheno$Unit), 
                statistic = 'r_GrpAve', NumofMU = 'Pairwise')
CR0  <- gca(Kmatrix = G, Xmatrix = X2, sigma2a = sigma2a, sigma2e = sigma2e,
             MUScenario = as.factor(cattle.pheno$Unit), statistic = 'CR0',
             NumofMU = 'Pairwise')
CR1  <- gca(Kmatrix = G, Xmatrix = X2, sigma2a = sigma2a, sigma2e = sigma2e,
             MUScenario = as.factor(cattle.pheno$Unit), statistic = 'CR1',
             NumofMU = 'Pairwise')
CR2  <- gca(Kmatrix = G, Xmatrix = X2, sigma2a = sigma2a, sigma2e = sigma2e,
               MUScenario = as.factor(cattle.pheno$Unit), statistic = 'CR2',
               NumofMU = 'Pairwise', Uidx = 8)
df_r_CR <- data.frame(r_GrpAve = r_GrpAve[upper.tri(r_GrpAve)], 
                      CR0 = CR0[upper.tri(CR0)], CR1 = CR1[upper.tri(CR1)],
                      CR2 = CR2[upper.tri(CR2)])
corrplot(cor(df_r_CR), method = "circle", type = "upper", 
         diag = F, number.cex= 4 / ncol(df_r_CR), col = cm.colors(10),
         number.digits = 6, addCoef.col = "black", tl.col = "black", 
         tl.srt = 90, tl.cex = 1.2, cl.lim = c(0.8, 1))

## -----------------------------------------------------------------------------
# Define a function using gca function in GCA
GC <- function(MU, Kmatrix, sigma2a, sigma2e, stat){
  X <- model.matrix(~ -1 + factor(MU)) 
  GCstat <- gca(Kmatrix, Xmatrix = X, sigma2a, sigma2e,
                MUScenario = as.factor(MU), 
                statistic = stat, NumofMU = 'Overall')
  return(GCstat)
}
# 5 simulated unit scenarios
head(CattleSIM[, 7: 11])
# PEVD_GrpAve
PEVD_GrpAve <- apply(CattleSIM[, 7: 11], 2, function(x) GC(x,
                     Kmatrix = G, sigma2a, sigma2e, stat = 'PEVD_GrpAve'))
print(PEVD_GrpAve)
# CD_GrpAve
CD_GrpAve <- apply(CattleSIM[, 7: 11], 2, function(x) GC(x,
                   Kmatrix = G, sigma2a, sigma2e, stat = 'CD_GrpAve'))
print(CD_GrpAve)

## -----------------------------------------------------------------------------
# Define a function to estimate PA using BGLR package
PA <- function(df, MU, Iter, BurnIn){
  folds <- length(unique(MU))
  Cor_CV <- matrix(NA, nrow = folds, ncol = 1)
  for(fold in 1:folds){ 
    y <- df
    yNa <- y
    whichNa <- which(MU == fold)
    yNa[whichNa] <- NA
    fit.GBLUP <- BGLR(y = yNa, ETA = ETA_GBLUP, nIter = Iter, burnIn = BurnIn, verbose = F)
    Cor_CV[fold, 1] <- cor(fit.GBLUP$yHat[whichNa], y[whichNa])
  }
  return(mean(Cor_CV))
}
# Prediction accuracies
EVD <- eigen(G)
ETA_GBLUP <- list(list(V = EVD$vectors, d = EVD$values, model = 'RKHS'))
PASum <- apply(CattleSIM[, 7: 11], 2, function(x) PA(df = CattleSIM$Phen, MU = x, 6000, 2000))



## -----------------------------------------------------------------------------
dfPEVD <- data.frame("PEVD" = PEVD_GrpAve, 'PA' = PASum, 
"Scenario" = paste('S', 1:5, sep = ''))
dfCD <- data.frame("CD" = CD_GrpAve, 'PA' = PASum, 
"Scenario" = paste('S', 1:5, sep = ''))
p1 <- ggplot(dfPEVD,  aes(PEVD, PA)) + geom_point(aes(shape = factor(Scenario)), color = "blue", size = 3) + 
  xlab("PEVD") + labs(shape="Scenario")
p2 <- ggplot(dfCD,  aes(CD, PA)) + geom_point(aes(shape = factor(Scenario)), color="red", size = 3) + 
  labs(shape = "Scenario")
grid.arrange(p1, p2, nrow = 1, widths = c(10, 10), heights = 4)


