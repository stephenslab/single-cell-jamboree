# Analyze the budding yeast data using flashier and fastTopics.
#
# sinteractive --mem=16G -c 8 --time=12:00:00 -p mstephens \
#   --account=pi-mstephens
# module load R/4.2.0
# .libPaths()[1]
# /home/pcarbo/R_libs_4_20
#
library(tools)
library(Matrix)
library(fastTopics)          # 0.7.25
library(flashier)            # 1.0.56
library(singlecelljamboreeR) # 0.1.41
load("../data/yeast.RData")
set.seed(1)

# Set up the "timings" data structure.
timings <- list(tm = 0,fl_nmf = 0)

# Compute the shifted log counts.
a <- 1
s <- rowSums(counts)
s <- s/mean(s)
shifted_log_counts <- MatrixExtra::mapSparse(counts/(a*s),log1p)

# Set a lower bound on the variances.
n  <- nrow(counts)
x  <- rpois(1e7,1/n)
s1 <- sd(log(x + 1))

# Fit a Poisson NMF using fastTopics.
# This step is expected to take about 2 h.
t0 <- proc.time()
tm0 <- fit_poisson_nmf(counts,k = 15,numiter = 200,method = "em",
                      control = list(numiter = 4,nc = 8,extrapolate = FALSE),
                      init.method = "random",verbose = "detailed")
tm <- fit_poisson_nmf(counts,fit0 = tm0,numiter = 200,method = "scd",
                      control = list(numiter = 4,nc = 8,extrapolate = TRUE),
                      verbose = "detailed")
t1 <- proc.time()
timings$tm <- t1 - t0

# Fit an NMF using flashier.
# This step is expected to take about 2 h.
t0 <- proc.time()
fl_nmf <- flashier_nmf(shifted_log_counts,k = 15,greedy_init = TRUE,
                       var_type = 2,S = s1,verbose = 2,maxiter = 200)
t1 <- proc.time()
timings$fl_nmf <- t1 - t0

# Fit an NMF with "cross-cutting" factors.
# TO DO.

# Save the model fits to an .Rdata file.
session_info <- sessionInfo()
fl_nmf_ldf <- ldf(fl_nmf,type = "i")
save(list = c("tm","fl_nmf_ldf","timings","session_info"),
     file = "yeast_factors.RData")
resaveRdaFiles("yeast_factors.RData")
