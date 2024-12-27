library(Matrix)
library(flashier)
library(fastTopics)
library(reshape2)
library(ggplot2)
library(cowplot)
source("../code/annotation_plots.R")
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

# Fit a semi-NMF using flashier.
fl0 <- flash(Y,ebnm_fn = c(ebnm_point_exponential,ebnm_point_normal),
             var_type = 0,greedy_Kmax = 9,nullcheck = FALSE,
             backfit = FALSE,verbose = 3)
fl <- flash_init(Y,var_type = 2,S = s1)
fl <- flash_factors_init(fl,fl0,
                         c(ebnm_point_exponential,
                           ebnm_point_normal))
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
p1 <- structure_plot(L[,-c(1,5)],grouping = celltype,gap = 20,
                     perplexity = 70,n = Inf) +
  labs(y = "membership",fill = "factor",color = "factor")
print(p1)

# TO DO: Explain what this function does, and how to use it.
driving_genes_heatmap <- function (effects_matrix,
                                   factors = 1:ncol(effects_matrix),
                                   n = 3,
                                   genes = NULL) {
}

# Create heatmaps to summarize the "driving genes" for each factor.
F <- with(fl_ldf,F %*% diag(D))
colnames(F) <- paste0("k",1:k)
genes <- NULL
ks <- c(2:4,6:9)
for (i in ks) {
  genes <- c(genes,
             head(order(F[,i],decreasing = TRUE),n = 3),
             head(order(F[,i],decreasing = FALSE),n = 3))
}
genes <- rownames(F)[genes]
genes <- unique(genes)
p2 <- effect_plot(F[genes,])
