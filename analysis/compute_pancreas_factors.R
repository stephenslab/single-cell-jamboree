# sinteractive 
library(Matrix)
library(sparseMatrixStats)
library(NNLM)
library(flashier)
load("../data/pancreas.RData")

set.seed(1)

# Compute the shifted log counts.
a <- 1
s <- rowSums(counts)
s <- s/mean(s)
Y <- MatrixExtra::mapSparse(counts/(a*s),log1p)

# Remove genes with very low variance in expression.
x <- colSds(Y)
j <- which(x > 0.01)
Y <- Y[,j]

# Set lower bound on variances.
n  <- nrow(counts)
x  <- rpois(1e7,1/n)
s1 <- sd(log(x + 1))

# Fit an NMF using NNLM.
Y_dense <- as.matrix(Y)
t0  <- proc.time()
nmf <- nnmf(Y_dense,k = 25,loss = "mse",method = "scd",max.iter = 200,
            verbose = 2,n.threads = 8)
t1  <- proc.time()
print(t1 - t0)

# Fit an NMF using flashier.
t0 <- proc.time()
fl0 <- flash(Y,ebnm_fn = ebnm_point_exponential,var_type = 0,
             greedy_Kmax = 40,nullcheck = FALSE,backfit = FALSE,
             verbose = 3)
fl_nmf <- flash_init(Y,var_type = 2,S = s1)
fl_nmf <- flash_factors_init(fl_nmf,fl0,ebnm_point_exponential)
fl_nmf <- flash_backfit(fl_nmf,extrapolate = FALSE,maxiter = 100,verbose = 3)
fl_nmf <- flash_backfit(fl_nmf,extrapolate = TRUE,maxiter = 100,verbose = 3)
t1 <- proc.time()
print(t1 - t0)
