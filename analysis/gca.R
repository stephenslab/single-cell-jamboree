# First look at the results of fitting NMF models to the combined GCA
# data. (If promising, I will migrate this analysis to an .Rmd file.)
library(ggplot2)
library(cowplot)
library(fastTopics)
load("../output/gca_nmf_k=20.RData")
L <- fl_nmf_ldf$L
k <- ncol(L)
colnames(L) <- paste0("k",1:k)
pdat <- cbind(sample_info,L)
# Top candidates for being shared by tissues and organoids:
# k = 2, *5*, *7*, 9, 10, 11, 16, 19, 
p1 <- ggplot(subset(pdat,k9 > 0.01),aes(x = k9)) +
  geom_histogram(fill = "black",color = "white",bins = 24) +
  facet_grid(rows = vars(origin),scales = "free_y") +
  theme_cowplot(font_size = 12)
