# Structure plots.
L <- poisson2multinom(pnmf)$L
newsgroups_topics <- paste0("k",c(1,3,5:10))
other_topics      <- paste0("k",c(2,4))
p1 <- structure_plot(L,grouping = topics,gap = 10,topics = newsgroups_topics) +
  labs(title = "newsgroups topics")
p2 <- structure_plot(L,grouping = topics,gap = 10,topics = other_topics) +
  labs(title = "other topics")
plot_grid(p1,p2,nrow = 2,ncol = 1)

#
# TO DO NEXT:
# - Move this code to newsgroups_annotate.Rmd.
# - Generate the "effects matrices" (F_le, F_vsnull), that we will analyze.
# - Create "distinctive features" scatterplots.
# - Design the "effects plot" on paper first (no R code).
#

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
