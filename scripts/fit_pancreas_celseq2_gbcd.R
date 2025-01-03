# remotes::install_github("stephenslab/gbcd@form-YYT-option")
# packageVersion("gbcd")
# 0.2.4
# 
library(Matrix)
library(ebnm)
library(flashier)
library(gbcd)
library(fastTopics)
library(ggplot2)
library(cowplot)
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
fl_cd <- fit_gbcd(Y,Kmax = 6,form_YYT = TRUE,
                  prior = flash_ebnm(prior_family = "point_exponential"),
                  maxiter1 = 100,maxiter2 = 100,maxiter3 = 100,
                  verbose = 3)

stop()

# Save the model fits to an .Rdata file.
fl_cd_ldf <- list(L = fl_cd$L,F = fl_cd$F$lfc)

# Structure plot.
celltype <- sample_info$celltype
celltype <-
 factor(celltype,
        c("acinar","ductal","activated_stellate","quiescent_stellate",
		  "endothelial","macrophage","mast","schwann","alpha","beta",
		  "delta","gamma","epsilon"))
L <- fl_cd_ldf$L
k <- ncol(L)
colnames(L) <- paste0("k",1:k)
celltype_factors  <- c(2,3,4,5,6,7,8)
other_factors <- c(1,9)
p1 <- structure_plot(L[,celltype_factors],grouping = celltype,
                     gap = 20,perplexity = 70,n = Inf) +
  labs(y = "membership",fill = "factor",color = "factor",
       title = "cell-type factors")
p2 <- structure_plot(L[,other_factors],grouping = celltype,
                     gap = 20,perplexity = 70,n = Inf) +
  scale_color_manual(values = other_colors) +
  scale_fill_manual(values = other_colors) +
  labs(y = "membership",fill = "factor",color = "factor",
       title = "other factors")
plot_grid(p1,p2,nrow = 2,ncol = 1)
