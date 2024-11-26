# flashier, NMF with cross-cutting factors (NMF-CC).
set.seed(1)
L <- fl_nmf_cc_ldf$L
k <- ncol(L)
colnames(L) <- paste0("k",1:k)
print(colSums(L > 0.1))
batch_topics    <- paste0("k",c(1,2,3,4,5,6,7,8,13))
celltype_topics <- paste0("k",c(9,12,14:21))
other_topics    <- paste0("k",c(10,11,22:23))
p1 <- structure_plot(L,topics = batch_topics,grouping = sample_info$tech,
                     gap = 10,perplexity = 70)
p2 <- structure_plot(L[cells,],topics = celltype_topics,
                     grouping = sample_info[cells,"celltype"],
                     gap = 20,n = Inf,perplexity = 70)
p3 <- structure_plot(L[cells,],topics = other_topics,
                     grouping = sample_info[cells,"celltype"],
                     gap = 20,n = Inf,perplexity = 70)
plot_grid(p1,p2,p3,nrow = 3,ncol = 1)
