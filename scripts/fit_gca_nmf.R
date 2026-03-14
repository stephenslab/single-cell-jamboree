# Fit an NMF to the combined GCA data using flashier.
#
# sinteractive --mem=40G -c 8 --time=72:00:00 -p mstephens \
#   --account=pi-mstephens
# module load R/4.2.0
# .libPaths()[1]
# /home/pcarbo/R_libs_4_20
#
library(tools)
library(Matrix)
library(ebnm)                # 1.1.42
library(flashier)            # 1.0.58
library(singlecelljamboreeR) # 0.1.41
load("../data/gca.RData")
set.seed(1)
k <- 30
outfile <- sprintf("gca_nmf_k=%d.RData",k)
cat("k =",k,"\n")
cat("outfile =",outfile,"\n")

# Filter out genes that are expressed in fewer than 10 cells.
j <- which(colSums(counts > 0) > 9)
counts <- counts[,j]

# Compute the shifted log counts.
a <- 1
s <- rowSums(counts)
s <- s/mean(s)
shifted_log_counts <- MatrixExtra::mapSparse(counts/(a*s),log1p)

# Set a lower bound on the variances.
n  <- nrow(counts)
x  <- rpois(1e7,1/n)
s1 <- sd(log(x + 1))

# Fit an NMF using flashier.
# This step is expected to take about 30 h with k = 20.
t0 <- proc.time()
fl_nmf <- flashier_nmf(shifted_log_counts,k = k,greedy_init = TRUE,
                       var_type = 2,S = s1,verbose = 2,maxiter = 200)
t1 <- proc.time()
print(t1 - t0)

# Save the model fit to an .Rdata file.
session_info <- sessionInfo()
fl_nmf_ldf <- ldf(fl_nmf,type = "i")
save(list = c("fl_nmf_ldf","session_info"),
     file = outfile)
resaveRdaFiles(outfile)
