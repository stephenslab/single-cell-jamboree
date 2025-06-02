library(Matrix)
library(rsvd)
library(fastTopics)
library(fastglmpca)
library(NNLM)
library(ggplot2)
library(cowplot)

data(pbmc_facs,package = "fastTopics")
counts  <- pbmc_facs$counts
samples <- pbmc_facs$samples
genes   <- pbmc_facs$genes

# Remove genes that are expressed in fewer than 5 cells.
j <- which(colSums(counts > 0) > 4)
counts <- counts[,j]
genes  <- genes[j,]

# Compute the shifted log counts.
a <- 1
s <- rowSums(counts)
s <- s/mean(s)
Y <- log1p(counts/(a*s))

# Get the ribosomal protein genes.
x       <- substr(genes$symbol,1,3)
rpgenes <- which(x == "RPS" | x == "RPL")

# Get proportion of expression due to RP Genes.
s   <- rowSums(counts)
srp <- rowSums(counts[,rpgenes])

# # PCA.
# subpop_colors <- c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e")
# pca <- rpca(Y,k = 10,center = TRUE,scale = FALSE)
# colnames(pca$x) <- paste0("PC",1:10)
# pdat <- data.frame(samples,pca$x)
# p1 <- ggplot(pdat,aes(x = PC1,y = PC2,color = subpop)) +
#   geom_point() +
#   scale_color_manual(values = subpop_colors) +
#   theme_cowplot(font_size = 12)

# # fastglmpca
# fit_glmpca <- init_glmpca_pois(t(counts),K = 10)
# fit_glmpca <- fit_glmpca_pois(t(counts),fit0 = fit_glmpca,
#                               verbose = TRUE,control = list(maxiter = 20))

# fastTopics
set.seed(1)
fit_tm <- fit_topic_model(counts,k = 6,
                          control.main = list(nc = 8,numiter = 4),
                          control.refine = list(nc = 8,numiter = 4,
                                                extrapolate = TRUE))
p1 <- structure_plot(fit_tm,gap = 25,grouping = samples$subpop)
plot(fit_tm$L[,"k4"],srp/s,pch = 20)

# NMF
nmf <- nnmf(as.matrix(Y),k = 6,loss = "mse",method = "scd",
            max.iter = 100,verbose = 2,n.threads = 8)
scale_cols <- function (A, b)
  t(t(A) * b)
W <- nmf$W
k <- ncol(W)
d <- apply(W,2,max)
W <- scale_cols(W,1/d)
p2 <- structure_plot(W,gap = 25,grouping = samples$subpop)
plot(W[,1],srp/s,pch = 20)
round(cor(fit_tm$L,W[,c(6,5,2,3,1,4)]),digits = 3)
