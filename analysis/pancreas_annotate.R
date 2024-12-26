library(tools)
library(Matrix)
library(ebnm)
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

# Create a scatterplot to explore the "driving genes" for the beta cells.

# Create a heatmap to explore the "driving genes" for the beta cells
# factor (k = 4).
#
# TO DO: Work this into a function that we can use to annotate the
# topics.
# 
k <- 9
F <- with(fl_ldf,F %*% diag(D))
genes       <- colnames(Y)
rownames(F) <- genes
colnames(F) <- paste0("k",1:k)
k    <- 4
pdat <- cbind(data.frame(gene = genes,stringsAsFactors = FALSE),
              F)
x    <- F[,k] # - apply(F[,c(2,3,5:9)],1,max)
rows <- order(x,decreasing = TRUE)
pdat <- pdat[rows,]
pdat$gene <- factor(pdat$gene,rev(pdat$gene))
i    <- 1:20
pdat2 <- melt(pdat[i,],id.vars = "gene",variable.name = "factor",
              value.name = "lfc")
p2 <- ggplot(pdat2,aes(x = factor,y = gene,size = lfc)) +
  geom_point(color = "white",fill = "royalblue",shape = 21) +
  scale_size(range = c(0.5,6)) +
  labs(x = "",y = "") +
  theme_cowplot(font_size = 9)

# TO DO: Produce a similar function for annotating the fastTopics and
# flashier semi-NMF results.

#
# *** OLD STUFF BELOW ***
#

# We can view the entries of the F matrix in the non-negative matrix
# factorization as LFCs (or, more correctly, log-fold *increases* in
# expression). What I'd like to do now is re-estimate the LFCs without
# the variational approximation, using lm() followed by ebnm().
m <- ncol(Y)
F <- matrix(0,m,k)
factors     <- colnames(L)
rownames(F) <- colnames(Y)
colnames(F) <- factors
fl_stats <- list(F_lse  = F,
                 F_se   = F,
                 F_pm   = F,
                 F_psd  = F,
                 F_lfsr = F)

# TO DO:
# + Explain why we rescale the memberships to be between 0 and 1.
# + Extract the 
fl_ldf <- ldf(fl,type = "i")
L <- fl_ldf$L
colnames(L) <- paste0("k",1:k)
for (i in 1:m) {
  y     <- Y[,i]
  dat   <- as.data.frame(cbind(y,L))
  # TO DO:
  # + Check these calculations.
  # + Remove the intercept.
  fit   <- lm(y ~ .,dat) 
  coefs <- summary(fit)$coefficients
  fl_stats$F_est[i,] <- coefs[factors,"Estimate"]
  fl_stats$F_se[i,]  <- coefs[factors,"Std. Error"]
}
for (j in 1:k) {
  x  <- fl_stats$F_est[,j]
  se <- fl_stats$F_se[,j]
  fit <- ebnm(x,se,prior_family = "point_exponential")
  fl_stats$F_pm[,j]  <- fit$posterior$mean
  fl_stats$F_psd[,j] <- fit$posterior$sd
}
k <- 2
plot(fl$F_pm[,k],fl_stats$F_pm[,k],pch = 20)
plot(fl$F_psd[,k],fl_stats$F_psd[,k],pch = 20,log = "xy")
# TO DO:
# + Fix the scaling for these comparisons.
# + Compute and store the lfsr's, too.
# + Compare the z-scores.
