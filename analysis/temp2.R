# Code to be incorporated into pancreas_cytokine_S1_factors.Rmd.
celltype_topics <- c("k2+k11","k3","k5","k7","k8","k9","k13")
lfc <- de$postmean
rownames(lfc) <- genes$symbol
lfc <- lfc[,celltype_topics]
out <- rank_effects(lfc)
rows <- NULL
for (i in colnames(lfc))
  rows <- c(rows,which(out[,i] <= 6))
rows <- unique(rows)
lfc <- lfc[rows,]
p1 <- effects_heatmap(lfc,font_size = 9)
print(p1)

celltype_factors <- c("k1","k2","k3","k4","k5","k9","k10","k11","k13")
F <- ldf(fl_nmf,type = "i")$F
rownames(F) <- genes$symbol
colnames(F) <- paste0("k",1:13)
F <- F[,celltype_factors]
out <- rank_effects(F)
rows <- NULL
for (i in colnames(F))
  rows <- c(rows,which(out[,i] <= 6))
rows <- unique(rows)
p2 <- effects_heatmap(F[rows,],font_size = 9)

out <- rank_effects(compute_le_effects(F))
rows <- NULL
for (i in colnames(F))
  rows <- c(rows,which(out[,i] <= 6))
rows <- unique(rows)
p3 <- effects_heatmap(F[rows,],font_size = 9)
print(plot_grid(p2,p3,nrow = 1,ncol = 2))
