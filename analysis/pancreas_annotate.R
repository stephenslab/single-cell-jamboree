library(tools)
library(Matrix)
library(ebnm)
library(flashier)
library(fastTopics)
library(ggplot2)
library(ggrepel)
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

# Interpret the factors with "distinctive genes" scatterplots.
F <- with(fl_ldf,F %*% diag(D))
colnames(F) <- paste0("k",1:k)
compare_dims <- c(2:3,5:9)
p1 <- distinctive_genes_scatterplot(F,k = 2,compare_dims = compare_dims,
                                    label_gene = function(x,y) x > 4 | y > 4,
                                    font_size = 9,label_size = 2) +
  ggtitle("factor 2 (acinar cells)") +
  theme(plot.title = element_text(size = 9,face = "plain"))
p2 <- distinctive_genes_scatterplot(F,k = 9,compare_dims = compare_dims,
                                    label_gene = function(x,y) x > 3.5 | y > 3,
                                    font_size = 9,label_size = 2) +
  ggtitle("factor 9 (ductal cells)") +
  theme(plot.title = element_text(size = 9,face = "plain"))
p3 <- distinctive_genes_scatterplot(F,k = 3,compare_dims = compare_dims,
                                    label_gene = function(x,y) x > 3.5 | y > 3,
                                    font_size = 9,label_size = 2) +
  ggtitle("factor 3 (stellate cells)") +
  theme(plot.title = element_text(size = 9,face = "plain"))
p4 <- distinctive_genes_scatterplot(F,k = 5,compare_dims = compare_dims,
                                    label_gene = function(x,y) x > 3.5 | y > 3,
                                    font_size = 9,label_size = 2) +
  ggtitle("factor 5 (alpha cells)") +
  theme(plot.title = element_text(size = 9,face = "plain"))
p5 <- distinctive_genes_scatterplot(F,k = 4,compare_dims = compare_dims,
                                    label_gene = function(x,y) x > 3 | y > 2,
                                    font_size = 9,label_size = 2) +
  ggtitle("factor 4 (beta cells)") +
  theme(plot.title = element_text(size = 9,face = "plain"))
p6 <- distinctive_genes_scatterplot(F,k = 6,compare_dims = compare_dims,
                                    label_gene = function(x,y) x > 3 | y > 2,
                                    font_size = 9,label_size = 2) +
  ggtitle("factor 6 (delta + gamma cells)") +
  theme(plot.title = element_text(size = 9,face = "plain"))
plot_grid(p1,p2,p3,p4,p5,p6,nrow = 2,ncol = 3)

le_lfc <- compute_le_diff(F,compare_dims = c(2:3,5:9))

# Create a scatterplot to explore the "driving genes" for the beta cells.

# Create a heatmap to explore the "driving genes" for the beta cells
# factor (k = 4).
#
# TO DO: Work this into a function that we can use to annotate the
# topics.
# 
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

