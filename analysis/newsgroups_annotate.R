library(Matrix)
library(flashier)
library(fastTopics)
library(ggplot2)
library(cowplot)
load("../data/newsgroups.RData")
load("../output/newsgroups_factors.RData")
set.seed(1)
L <- poisson2multinom(pnmf)$L
newsgroups_topics <- paste0("k",c(1:3,5:10))
other_topics      <- paste0("k",c(2,4))
p1 <- structure_plot(L,grouping = topics,gap = 10,topics = newsgroups_topics) +
  labs(title = "newsgroups topics")
p2 <- structure_plot(L,grouping = topics,gap = 10,topics = other_topics) +
  labs(title = "other topics")
print(plot_grid(p1,p2,nrow = 2,ncol = 1))
