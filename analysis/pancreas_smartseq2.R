# This is just a first draft of some analysis code that will be
# incorporated into pancreas_closer_look.Rmd.
set.seed(1)
load("../data/pancreas.RData")
load("../output/pancreas_smartseq2_factors.RData")
i           <- which(sample_info$tech == "smartseq2")
sample_info <- sample_info[i,]
counts      <- counts[i,]
sample_info <- transform(sample_info,celltype = factor(celltype))

celltype <- sample_info$celltype
celltype <-
 factor(celltype,
        c("acinar","ductal","activated_stellate","quiescent_stellate",
		  "endothelial","macrophage","mast","schwann","t_cell","alpha",
		  "beta","delta","gamma","epsilon"))

# (1) Topic model (fastTopics).
L <- poisson2multinom(pnmf)$L
p1 <- structure_plot(L,grouping = celltype,
                     gap = 20,perplexity = 70,n = Inf)

# (2) NMF (flashier)
L <- fl_nmf_ldf$L
k <- ncol(L)
colnames(L) <- paste0("k",1:k)
celltype_topics  <- 2:7
other_topics <- c(1,8,9)
p2 <- structure_plot(L[,celltype_topics],grouping = celltype,gap = 20,
                     perplexity = 70,n = Inf)
p3 <- structure_plot(L[,other_topics],grouping = celltype,gap = 20,
                     perplexity = 70,n = Inf)
plot_grid(p2,p3,nrow = 2,ncol = 1)

# (3) NMF (NNLM)
scale_cols <- function (A, b)
  t(t(A) * b)
W <- nmf$W
k <- ncol(W)
d <- apply(W,2,max)
W <- scale_cols(W,1/d)
colnames(W) <- paste0("k",1:k)
celltype_topics  <- c(3:6,8,9)
other_topics <- c(1,2,7)
p4 <- structure_plot(W[,celltype_topics],grouping = celltype,
                     gap = 20,perplexity = 70,n = Inf)
p5 <- structure_plot(W[,other_topics],grouping = celltype,
                     gap = 20,perplexity = 70,n = Inf)
plot_grid(p4,p5,nrow = 2,ncol = 1)

# (4) semi-NMF (flashier)
L <- fl_snmf_ldf$L
k <- ncol(L)
colnames(L) <- paste0("k",1:k)
celltype_topics  <- c(3,4,6,7,9)
other_topics <- c(1,5,8)
p6 <- structure_plot(L[,celltype_topics],grouping = celltype,
                     gap = 20,perplexity = 70,n = Inf)
p7 <- structure_plot(L[,other_topics],grouping = celltype,
                     gap = 20,perplexity = 70,n = Inf)
plot_grid(p6,p7,nrow = 2,ncol = 1)
