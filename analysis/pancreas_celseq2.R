# This is just a first draft of some analysis code that will go into
# pancreas_closer_look.Rmd.
library(Matrix)
library(fastTopics)
library(ggplot2)
library(cowplot)
set.seed(1)
load("../data/pancreas.RData")
load("../output/pancreas_celseq2_factors.RData")
i           <- which(sample_info$tech == "celseq2")
sample_info <- sample_info[i,]
counts      <- counts[i,]
sample_info <- transform(sample_info,celltype = factor(celltype))
