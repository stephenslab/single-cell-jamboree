# This is just a first draft of some analysis code that will go into
# pancreas_closer_look.Rmd.

# Topic model (fastTopics).
L <- poisson2multinom(pnmf)$L
p1 <- structure_plot(L,grouping = sample_info$celltype,
                     gap = 10,perplexity = 70,n = Inf)

# NMF (NNLM)
scale_cols <- function (A, b)
  t(t(A) * b)
W <- nmf$W
k <- ncol(W)
d <- apply(W,2,max)
W <- scale_cols(W,1/d)
colnames(W) <- paste0("k",1:k)
p2 <- structure_plot(W,grouping = sample_info$celltype,
                     gap = 10,perplexity = 70,n = Inf)

# NMF (flashier)
L <- fl_nmf_ldf$L
k <- ncol(L)
colnames(L) <- paste0("k",1:k)
p3 <- structure_plot(L[,-1],grouping = sample_info$celltype,
                     gap = 10,perplexity = 70,n = Inf)

# semi-NMF (flashier)
L <- fl_snmf_ldf$L
k <- ncol(L)
colnames(L) <- paste0("k",1:k)
p4 <- structure_plot(L[,-1],grouping = sample_info$celltype,
                     gap = 10,perplexity = 70,n = Inf)

plot_grid(p1,p2,p3,p4,nrow = 4,ncol = 1)
