library(ggrepel)

# Compute the "least extreme" (l.e.) effect differences for a
# non-negative effects matrix.
compute_le_diff <- function (effects_matrix,
                             compare_dims = seq(1,ncol(effects_matrix))) {
  m <- ncol(effects_matrix)
  out <- effects_matrix
  for (i in 1:m) {
    dims <- setdiff(compare_dims,i)
    out[,i] <- effects_matrix[,i] - apply(effects_matrix[,dims],1,max)
  }
  return(out)
}

# This creates a "distinctive genes plot"; this is a plot in which the
# effect estimate is shown on the x axis and the "least extreme"
# difference between the estimated effects is shown on the y axis. The
# idea is that these scatterplots should better highlight the
# "interesting" genes for a given dimension/factor. The "label_gene"
# input argument is a function that returns TRUE when the gene should
# be labeled in the plot; the default is that it always returns FALSE
# (so that no genes are labeled in the plot).
distinctive_genes_scatterplot <- function (effects_matrix, k,
                                           effect_quantile_prob = 0.999,
                                           lediff_quantile_prob = 0.999) {
  lediff <- compute_le_diff(effects_matrix)
  genes  <- rownames(effects_matrix)
  pdat   <- data.frame(gene    = genes,
                       effect  = effects_matrix[,k],
                       lediff = lediff[,k])
  effect_quantile <- quantile(pdat$effect,effect_quantile_prob)
  lediff_quantile <- quantile(pdat$lediff,lediff_quantile_prob)
  i <- which(pdat$effect < effect_quantile & pdat$lediff < lediff_quantile)
  pdat[i,"gene"] <- NA
  return(ggplot(pdat,aes(x = effect,y = lediff,label = gene)) +
         geom_point(color = "dodgerblue") +
         geom_hline(yintercept = 0,color = "magenta",linetype = "dotted",
                    linewidth = 0.5) +
         geom_text_repel(color = "black",size = 2.25,
                         fontface = "italic",segment.color = "black",
                         segment.size = 0.25,min.segment.length = 0,
                         max.overlaps = Inf,na.rm = TRUE) +
         labs(x = "log-fold change",y = "l.e. difference") +
         theme_cowplot(font_size = 10))
}

F <- fl_nmf_ldf$F
colnames(F) <- paste0("k",1:9)
kset <- paste0("k",4:6)
p1 <- distinctive_genes_scatterplot(F[,kset],"k4") + ggtitle("factor k4")
p2 <- distinctive_genes_scatterplot(F[,kset],"k5") + ggtitle("factor k5")
p3 <- distinctive_genes_scatterplot(F[,kset],"k6") + ggtitle("factor k6")
print(plot_grid(p1,p2,p3,nrow = 1,ncol = 3))
