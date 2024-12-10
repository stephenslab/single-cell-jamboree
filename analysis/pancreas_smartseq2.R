# This is my first-pass analysis of the pancreas smartseq2 data using
# various matrix factorization methods (flashier, fastTopics, NNLM,
# etc), as well as variants of these methods.
library(tools)
library(Matrix)
library(flashier)
library(fastTopics)
load("../data/pancreas.RData")
set.seed(1)

# Select the Smart-seq2 data (Segerstolpe et al, 2016).
# This should select 2,394 cells.
i           <- which(sample_info$tech == "smartseq2")
sample_info <- sample_info[i,]
counts      <- counts[i,]

# Remove genes that are expressed in fewer than 10 cells.
x      <- colSums(counts > 0)
j      <- which(x > 9)
counts <- counts[,j]

# Compute the shifted log counts.
a <- 1
s <- rowSums(counts)
s <- s/mean(s)
Y <- MatrixExtra::mapSparse(counts/(a*s),log1p)

# Set a lower bound on the variances.
n  <- nrow(counts)
x  <- rpois(1e7,1/n)
s1 <- sd(log(x + 1))

# Fit an NMF using flashier.
# I initially set var_type = 0 to increase the number of
# factors discovered.
fl0 <- flash(Y,ebnm_fn = ebnm_point_exponential,var_type = 0,
             greedy_Kmax = 40,nullcheck = FALSE,backfit = FALSE,
             verbose = 3)
fl_nmf <- flash_init(Y,var_type = 2,S = s1)
fl_nmf <- flash_factors_init(fl_nmf,fl0,ebnm_point_exponential)
fl_nmf <- flash_backfit(fl_nmf,extrapolate = FALSE,maxiter = 100,verbose = 3)
fl_nmf <- flash_backfit(fl_nmf,extrapolate = TRUE,maxiter = 100,verbose = 3)
