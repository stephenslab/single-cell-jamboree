# Analyze the pancreas LSA data set using flashier and fastTopics.
# See prepare_pancreas_lsa_data.R for background on these data.
#
# sinteractive --mem=20G -c 8 --time=12:00:00 -p mstephens \
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
load("../data/pancreas_cytokine_lsa_v2.RData")
colnames(counts) <- genes$ensembl
set.seed(1)
# k <- 15
k <- 13

# Filter out genes that are expressed in fewer than 10 cells.
j      <- which(colSums(counts > 0) > 9)
genes  <- genes[j,]
counts <- counts[,j]

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
tm0 <- fit_poisson_nmf(counts,k = k,numiter = 200,method = "em",
                      control = list(numiter = 4,nc = 8,extrapolate = FALSE),
                      init.method = "random",verbose = "detailed")
tm <- fit_poisson_nmf(counts,fit0 = tm0,numiter = 200,method = "scd",
                      control = list(numiter = 4,nc = 8,extrapolate = TRUE),
                      verbose = "detailed")
t1 <- proc.time()
print(t1 - t0)
timings$tm <- t1 - t0

# Fit an NMF using flashier.
# This step is expected to take about 3 h.
t0 <- proc.time()
fl_nmf <- flashier_nmf(shifted_log_counts,k = k,greedy_init = TRUE,
                       var_type = 2,S = s1,verbose = 2,maxiter = 200)
t1 <- proc.time()
print(t1 - t0)
timings$fl_nmf <- t1 - t0

# Fit an NMF with "cross-context" factors.
cc_factors <- c(5,12)
fl_nmf_cc <- convert_factors_nn_to_pn(fl_nmf,cc_factors,shifted_log_counts,
                                      S = s1,var_type = 2)
fl_nmf_cc <- flash_backfit(fl_nmf_cc,extrapolate = FALSE,maxiter = 100,
                           verbose = 2)
fl_nmf_cc <- flash_backfit(fl_nmf_cc,extrapolate = TRUE,maxiter = 100,
                           verbose = 2)

# Save the model fits to an .Rdata file.
session_info <- sessionInfo()
fl_nmf_ldf <- ldf(fl_nmf,type = "i")
fl_nmf_cc_ldf <- ldf(fl_nmf_cc,type = "i")
save(list = c("tm","fl_nmf_ldf","fl_nmf_cc_ldf","timings","session_info"),
     file = "pancreas_lsa_factors.RData")
resaveRdaFiles("pancreas_lsa_factors.RData")
