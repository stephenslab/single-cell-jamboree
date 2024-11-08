library(Matrix)
library(fastTopics)
set.seed(1)
load("../data/pancreas.RData")
load("../output/pancreas_factors.RData")

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

# NMF.
cells <- subsample(sample_info$celltype,n = 500)
W <- nmf$W
k <- ncol(W)
d <- apply(W,2,max)
W <- scale.cols(W,1/d)
colnames(W) <- paste0("k",1:k)
# celseq = k17
# celseq2 = k13
# fluidigmc1 = k14
# inDrop1:4 = k7
# smarter = k6
# smartseq2 = k1
batch_topics <- c(1,6,7,10,13,14,16,17)
p1 <- structure_plot(W[,batch_topics],grouping = sample_info$tech,gap = 10)
celltype_topics <- c(2,5,8,9,11,12,15,20,21,22,23)
p2 <- structure_plot(W[cells,celltype_topics],
                     grouping = sample_info[cells,"celltype"],
                     gap = 20,n = Inf)

# flashier, NMF
# TO DO.

# flashier, semi-NMF
L <- fl_snmf_ldf$L
colnames(L) <- paste0("k",1:k)
batch_topics <- c(4,5,6,11,22)
p1 <- structure_plot([,batch_topics],grouping = sample_info$tech,gap = 10)
celltype_topics <- c(3,7,8,9,10,12,13,14,15,16,17,18,19,20,21,23)
p2 <- structure_plot(L[cells,celltype_topics],
                     grouping = sample_info[cells,"celltype"],
                     gap = 40,n = Inf)

