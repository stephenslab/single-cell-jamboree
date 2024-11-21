library(Matrix)
library(fastTopics)
library(ggplot2)
library(cowplot)
set.seed(1)
load("../data/pancreas.RData")
load("../output/pancreas_factors.RData")
timings0 <- timings
load("../output/pancreas_factors2.RData")
rm(session_info)
timings <- c(timings0,timings)

# scale.cols(A,b) scales each column A[,i] by b[i].
scale.cols <- function (A, b)
  t(t(A) * b)

# Function to subsample the cell types.
subsample <- function (x, n = 1000) {
  rows <- NULL
  groups <- levels(x)
  for (g in groups) {
    i  <-  which(x == g)
    n0 <- min(n,length(i))
    i  <- sample(i,n0)
    rows <- c(rows,i)
  }
  return(sort(rows))
}

# fastTopics.
cells <- subsample(sample_info$celltype,n = 500)
L <- poisson2multinom(pnmf)$L
batch_topics     <- c(2,5,6,11,12)
celltype_topics  <- c(4,9,8,15,16,17,18,19,20)
celltype_topics2 <- c(1,3,7,10,13,14,21,22,23)
p1 <- structure_plot(L[,batch_topics],grouping = sample_info[,"tech"],
                     gap = 10,perplexity = 70)
p2 <- structure_plot(L[cells,celltype_topics],
                     grouping = sample_info[cells,"celltype"],
                     gap = 20,perplexity = 70,n = Inf)
p3 <- structure_plot(L[cells,celltype_topics2],
                     grouping = sample_info[cells,"celltype"],
                     gap = 20,perplexity = 70,n = Inf)
plot_grid(p1,p2,p3,nrow = 3,ncol = 1)

stop()

W <- nmf$W
k <- ncol(W)
d <- apply(W,2,max)
W <- scale.cols(W,1/d)
colnames(W) <- paste0("k",1:k)
batch_topics <- c(1,6,7,10,13,14,16,17)
p1 <- structure_plot(W[,batch_topics],grouping = sample_info$tech,gap = 10,
                     perplexity = 70)
celltype_topics <- c(2,5,8,9,11,12,15,20,21,22,23)
p2 <- structure_plot(W[cells,celltype_topics],
                     grouping = sample_info[cells,"celltype"],
                     gap = 20,n = Inf,perplexity = 70)
plot_grid(p1,p2,nrow = 2,ncol = 1)

# NMF.
W <- nmf$W
k <- ncol(W)
d <- apply(W,2,max)
W <- scale.cols(W,1/d)
colnames(W) <- paste0("k",1:k)
batch_topics <- c(1,6,7,10,13,14,16,17)
p1 <- structure_plot(W[,batch_topics],grouping = sample_info$tech,gap = 10,
                     perplexity = 70)
celltype_topics <- c(2,5,8,9,11,12,15,20,21,22,23)
p2 <- structure_plot(W[cells,celltype_topics],
                     grouping = sample_info[cells,"celltype"],
                     gap = 20,n = Inf,perplexity = 70)
plot_grid(p1,p2,nrow = 2,ncol = 1)

# flashier, NMF.
set.seed(1)
L <- fl_nmf_ldf$L
k <- ncol(L)
colnames(L) <- paste0("k",1:k)
print(colSums(L > 0.1))
batch_topics    <- paste0("k",c(2,3,4,5,7,8,20))
celltype_topics <- paste0("k",c(6,11:19,21))
other_topics    <- paste0("k",c(1,9:10,22:23))
p1 <- structure_plot(L,topics = batch_topics,grouping = sample_info$tech,
                     gap = 10,perplexity = 70)
p2 <- structure_plot(L[cells,],topics = celltype_topics,
                     grouping = sample_info[cells,"celltype"],
                     gap = 20,n = Inf,perplexity = 70)
p3 <- structure_plot(L[cells,],topics = other_topics,
                     grouping = sample_info[cells,"celltype"],
                     gap = 20,n = Inf,perplexity = 70)
plot_grid(p1,p2,p3,nrow = 3,ncol = 1)

# flashier, semi-NMF.
L <- fl_snmf_ldf$L
colnames(L) <- paste0("k",1:k)
batch_topics <- c(4,5,6,11,22)
p1 <- structure_plot(L[,batch_topics],grouping = sample_info$tech,gap = 10,
                     perplexity = 70)
celltype_topics <- c(3,7,8,9,10,12,13,14,15,16,17,18,19,20,21,23)
p2 <- structure_plot(L[cells,celltype_topics],
                     grouping = sample_info[cells,"celltype"],
                     gap = 40,n = Inf,perplexity = 70)

