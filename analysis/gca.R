# First look at the results of fitting NMF models to the combined GCA
# data. (If promising, I will migrate this analysis to an .Rmd file.)
load("../output/gca_nmf_k=20.RData")
L <- fl_nmf_ldf$L
k <- ncol(L)
colnames(L) <- paste0("k",1:k)
pdat <- cbind(sample_info,L)
p1 <- ggplot(subset(pdat,k14 > 0.01),aes(x = k14)) +
  geom_histogram(fill = "black",color = "white",bins = 24) +
  facet_grid(rows = vars(origin),scales = "free_y") +
  theme_cowplot(font_size = 12)
p2 <- ggplot(subset(pdat,k14 > 0.01),aes(x = k14)) +
  geom_histogram(fill = "black",color = "white",bins = 24) +
  facet_wrap(facets = vars(celltype),scales = "free_y") +
  theme_cowplot(font_size = 8)
with(sample_info,table(celltype,origin))

