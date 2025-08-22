# Analyze the pancreas LSA data using flashier and fastTopics.
# See prepare_pancreas_lsa_data.R for background on these data.
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
load("../data/pancreas_cytokine_lsa_v2.RData")
set.seed(1)

# Filter out genes that are expressed in fewer than 10 cells.
j <- which(colSums(counts > 0) > 9)
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
#
# TO DO: Use k = 15.
#

# Fit an NMF using flashier.
#
# TO DO: Use k = 15.
#

# Fit an NMF with "cross-cutting" factors.
# TO DO.

# Save the model fits to an .Rdata file.
# TO DO.
