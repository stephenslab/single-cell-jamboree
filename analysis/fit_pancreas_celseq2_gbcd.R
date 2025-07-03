#fits gbcd to the pancreas celseq2 data and saves results

library("gbcd")
library("Matrix")

load("../data/pancreas.RData")
set.seed(1)

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

#randomly divide rows of Y into 2
subset = sample(1:nrow(Y), nrow(Y)/2)

fit.gbcd = fit_gbcd(Y, Kmax = 20)
fit.gbcd.1 = fit_gbcd(Y[subset,],Kmax = 20)
fit.gbcd.2 = fit_gbcd(Y[-subset,],Kmax=20)
session_info <- sessionInfo()

save(list = c("fit.gbcd","fit.gbcd.1","fit.gbcd.2","session_info"),
     file = "../output/pancreas_celseq2_gbcd.RData")

