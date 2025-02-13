# Analyze the "20 Newsgroups" data using various matrix factorization
# methods (flashier, fastTopics, NNLM).
#
# sinteractive --mem=24G -c 8 --time=24:00:00
# module load R/4.2.0
# .libPaths()[1]
# /home/pcarbo/R_libs_4_20
library(tools)
library(Matrix)
library(NNLM)
library(fastTopics)
library(flashier)
load("../data/newsgroups.RData")
set.seed(1)

# Remove words that appear in fewer than 10 documents.
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

# Set up the "timings" data structure.
timings <- list(nmf        = 0,
                fl_nmf     = 0,
                fl_snmf    = 0,
                fasttopics = 0)

# (1) Fit an NMF using flashier, with k = 10 factors.
t0 <- proc.time()
fl_nmf <- flash(Y,ebnm_fn = ebnm_point_exponential,var_type = 2,
                greedy_Kmax = 10,S = s1,nullcheck = FALSE,
                backfit = FALSE,verbose = 3)
fl_nmf <- flash_backfit(fl_nmf,extrapolate = FALSE,maxiter = 100,verbose = 3)
t1 <- proc.time()
timings$fl_nmf <- t1 - t0
print(timings$fl_nmf)

# (2) Fit an NMF using NNLM, with k = 10 factors.
Y_dense <- as.matrix(Y)
t0  <- proc.time()
nmf <- nnmf(Y_dense,k = 10,loss = "mse",method = "scd",
            max.iter = 200,verbose = 2,n.threads = 8)
t1  <- proc.time()
timings$nmf <- t1 - t0
print(timings$nmf)
rm(Y_dense)

# (3) Fit a semi-NMF using flashier, with k = 10 factors.
t0 <- proc.time()
fl_snmf <- flash(Y,ebnm_fn = c(ebnm_point_exponential,ebnm_point_normal),
                 var_type = 2,greedy_Kmax = 10,S = s1,nullcheck = FALSE,
                 backfit = FALSE,verbose = 3)
fl_snmf <- flash_backfit(fl_snmf,extrapolate = FALSE,maxiter = 100,verbose = 3)
t1 <- proc.time()
timings$fl_snmf <- t1 - t0
print(timings$fl_snmf)

# (4) Fit a Poisson NMF using fastTopics, with k = 10 factors/topics.
t0 <- proc.time()
pnmf <- fit_poisson_nmf(counts,k = 10,numiter = 100,method = "em",
                        control = list(numiter = 4,nc = 8,extrapolate = FALSE),
                        init.method = "random",verbose = "detailed")
pnmf <- fit_poisson_nmf(counts,fit0 = pnmf,numiter = 100,method = "scd",
                        control = list(numiter = 4,nc = 8,extrapolate = TRUE),
                        verbose = "detailed")
t1 <- proc.time()
timings$fasttopics <- t1 - t0
print(timings$fasttopics)

# Perform the "grade of membership" differential expression analysis
# using the fitted Poisson NMF model.
de_le <- de_analysis(pnmf,counts,shrink.method = "ash",
                     lfc.stat = "le",pseudocount = 0.1,
                     control = list(ns = 1e4,nc = 8,nsplit = 1000))
de_vsnull <- de_analysis(pnmf,counts,shrink.method = "ash",
                         lfc.stat = "vsnull",pseudocount = 0.1,
                         control = list(ns = 1e4,nc = 8,nsplit = 1000))

# Save the model fits and other outputs to an .Rdata file.
fl_nmf_ldf   <- ldf(fl_nmf,type = "i")
fl_snmf_ldf  <- ldf(fl_snmf,type = "i")
session_info <- sessionInfo()
save(list = c("nmf","fl_nmf_ldf","fl_snmf_ldf","pnmf","de_le",
              "de_vsnull","timings","session_info"),
     file = "newsgroups_factors.RData")
resaveRdaFiles("newsgroups_factors.RData")
