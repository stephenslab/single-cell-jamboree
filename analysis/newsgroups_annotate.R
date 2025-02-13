# TO DO NEXT, PART 2: Move this code to newsgroups_annotate.Rmd.
library(Matrix)
library(fastTopics)
library(ggplot2)
library(cowplot)
library(ggrepel)
source("../code/annotation_plots.R")
load("../data/newsgroups.RData")
load("../output/newsgroups_factors.RData")
stop()
set.seed(1)

# Remove words that appear in fewer than 10 documents.
x <- colSums(counts > 0)
j <- which(x > 9)
counts <- counts[,j]

# Topic model.
L <- poisson2multinom(pnmf)$L
newsgroups_topics <- paste0("k",c(1,3,5:10))
other_topics      <- paste0("k",c(2,4))
p1 <- structure_plot(L,grouping = topics,gap = 10,topics = newsgroups_topics) +
  labs(title = "newsgroups topics")
p2 <- structure_plot(L,grouping = topics,gap = 10,topics = other_topics) +
  labs(title = "other topics")
plot_grid(p1,p2,nrow = 2,ncol = 1)

# Volcano plots using "least extreme" LFCs.
p1 <- volcano_plot(de_le,k = "k1",ymax = 30)
p2 <- volcano_plot(de_le,k = "k2",ymax = 100)
p3 <- volcano_plot(de_le,k = "k3")
p5 <- volcano_plot(de_le,k = "k5",ymax = 80)
p6 <- volcano_plot(de_le,k = "k6")
p7 <- volcano_plot(de_le,k = "k7",ymax = 30)
p8 <- volcano_plot(de_le,k = "k8",ymax = 60)
p9 <- volcano_plot(de_le,k = "k9",ymax = 40)
p10 <- volcano_plot(de_le,k = "k10",ymax = 40)

# TO DO NEXT, PART 1a: Create plotting functions for visualizing and/or
# annotating the semi-NMF factors (fl_snmf_ldf). This includes "effect
# plots".
L <- fl_snmf_ldf$L
k <- ncol(L)
colnames(L) <- paste0("k",1:k)
newsgroups_topics <- c(1,2,4,5,6,7,8,9,10)
p1 <- structure_plot(L,grouping = topics,gap = 10,topics = newsgroups_topics) +
  labs(title = "newsgroups topics")
p2 <- structure_plot(L,grouping = topics,gap = 10,topics = other_topics) +
  labs(title = "other topics")
plot_grid(p1,p2,nrow = 2,ncol = 1)

#
# PART 1b: Create plotting functions for visualizing and/or annotating
# the flashier NMF factors (fl_nmf_ldf). This includes "driving gene
# scatterplots" and "effect plots". (Although I should replace "gene"
# with something else more general, e.g., feature.)
normalize.rows <- function (A)
  A / rowSums(A)
L <- fl_nmf_ldf$L
k <- ncol(L)
colnames(L) <- paste0("k",1:k)
newsgroups_topics <- paste0("k",c(1,2,3,5,6,9,10))
L <- L[,newsgroups_topics]
L <- normalize.rows(L)
p1 <- structure_plot(L,grouping = topics,gap = 10)
p2 <- structure_plot(L,grouping = topics,gap = 10,topics = other_topics) +
  labs(title = "other topics")
plot_grid(p1,p2,nrow = 2,ncol = 1)

# NNLM
L <- nmf$W
k <- ncol(L)
colnames(L) <- paste0("k",1:k)
newsgroups_topics <- c(1,2,3,4,5,6,7,8,9,10)
p1 <- structure_plot(L,grouping = topics,gap = 10,topics = newsgroups_topics) +
  labs(title = "newsgroups topics")
p2 <- structure_plot(L,grouping = topics,gap = 10,topics = other_topics) +
  labs(title = "other topics")
plot_grid(p1,p2,nrow = 2,ncol = 1)

