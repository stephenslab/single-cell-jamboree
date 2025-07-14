# Code to be incorporated into pancreas_cytokine_S1_factors.Rmd.
celltype_topics <- c("k2+k11","k3","k5","k7","k8","k9","k13")
lfc <- de$postmean
rownames(lfc) <- genes$symbol
lfc <- lfc[,celltype_topics]
out <- rank_effects(lfc)
rows <- NULL
for (i in colnames(lfc))
  rows <- c(rows,which(out[,i] <= 5))
lfc <- lfc[rows,]
p <- effects_heatmap(lfc)
print(p)
