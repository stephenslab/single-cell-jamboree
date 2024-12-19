library(tools)
library(Matrix)
library(flashier)
library(fastTopics)
library(ggplot2)
library(cowplot)
load("../data/pancreas.RData")
set.seed(1)

# Number of factors.
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

# Fit an NMF using flashier.
# I initially set var_type = 0 to increase the number of
# factors discovered.
fl0 <- flash(Y,ebnm_fn = ebnm_point_exponential,var_type = 0,
             greedy_Kmax = k,nullcheck = FALSE,backfit = FALSE,
             verbose = 3)
fl <- flash_init(Y,var_type = 2,S = s1)
fl <- flash_factors_init(fl,fl0,ebnm_point_exponential)
fl <- flash_backfit(fl,extrapolate = FALSE,maxiter = 100,verbose = 3)
fl <- flash_backfit(fl,extrapolate = TRUE,maxiter = 100,verbose = 3)

# Visualize the factors with a Structure plot.
# I omit the first factor from the Structure plot because it is a
# "baseline" factor.
celltype <- sample_info$celltype
celltype <-
 factor(celltype,
        c("acinar","ductal","activated_stellate","quiescent_stellate",
		  "endothelial","macrophage","mast","schwann","alpha","beta",
		  "delta","gamma","epsilon"))
fl_ldf <- ldf(fl,type = "i")
L <- fl_ldf$L
colnames(L) <- paste0("k",1:k)
p1 <- structure_plot(L[,-1],grouping = celltype,gap = 20,
                     perplexity = 70,n = Inf) +
  labs(y = "membership",fill = "factor",color = "factor")
print(p1)
