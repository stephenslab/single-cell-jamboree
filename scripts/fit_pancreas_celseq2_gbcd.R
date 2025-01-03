# remotes::install_github("stephenslab/gbcd@form-YYT-option")
# packageVersion("gbcd")
# 0.2.4
# 
library(Matrix)
library(ebnm)
library(flashier)
library(gbcd)
load("../data/pancreas.RData")
set.seed(1)

# All the matrix factorizations will have this many topics or factors.
k <- 9

# Select the CEL-seq2 data (Muraro et al, 2016).
# This should select 2,285 cells.
i           <- which(sample_info$tech == "celseq2")
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

# Fit a semi-NMF covariance decomposition.
# Note that Kmax = 5 results in a factorization with at most 9 factors.
fl_cd <- fit_gbcd(Y,Kmax = 5,form_YYT = TRUE,
                  prior = flash_ebnm(prior_family = "point_exponential"),
                  maxiter1 = 100,maxiter2 = 100,maxiter3 = 100,
                  verbose = 3)

# Save the model fits to an .Rdata file.
# fl_cd_ldf <- ldf(fl_cd$fit.cov,type = "i")
