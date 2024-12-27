# This function creates an "effect plot". The effects_matrix input
# should be a matrix in which rows are features and columns are
# dimensions (e.g., factors in a matrix factorization). Effects
# smaller than zero_value in magnitude are not included in the plot.
effect_plot <- function (effects_matrix, font_size = 9,
                         zero_value = 0.01) {
  features <- rownames(effects_matrix)
  pdat <- data.frame(feature_name = features,
                     stringsAsFactors = FALSE)
  pdat <- cbind(pdat,effects_matrix)
  pdat <- melt(pdat,id.vars = "feature_name",variable.name = "dim",
               value.name = "value")
  pdat <- transform(pdat,
                    effect_size = abs(value),
                    effect_sign = factor(sign(value),c(-1,0,1)),
                    feature_name = factor(feature_name,rev(features)))
  pdat$effect_size[pdat$effect_size < zero_value] <- NA
  return(ggplot(pdat,aes(x = dim,y = feature_name,size = effect_size,
                         fill = effect_sign)) +
         geom_point(color = "white",shape = 21,na.rm = TRUE) +
         scale_size(range = c(1,6),
                    breaks = range(pdat$effect_size,na.rm = TRUE),
                    labels = round(range(pdat$effect_size,na.rm = TRUE),
                                   digits = 2)) +
         scale_fill_manual(values = c("navy","gray","orangered"),
                           drop = FALSE) +
         guides(size = guide_legend(override.aes = list(fill = "black")),
                fill = guide_legend(override.aes = list(size = 3))) +
         labs(x = "dimension",y = "feature name",
              fill = "effect sign",size = "effect size") +
         theme_cowplot(font_size = font_size))
}

# This function selects the top "driving" genes for the selected
# dimensions ("dims"). The effects_matrix input should be a matrix in
# which rows are genes and columns are dimensions. It is assumed the
# names of the rows give the names or ids of the genes. If
# select_down_effects = TRUE, the genes with the largest negative
# effects are also included. Note that a gene is never selected more
# than once.
select_driving_genes <- function (effect_matrix, dims, n,
                                  select_down_effects = TRUE) {
  genes <- NULL
  for (i in dims) {
    genes <- c(genes,head(order(effect_matrix[,i],decreasing = TRUE),n))
    if (select_down_effects)
      genes <- c(genes,head(order(effect_matrix[,i],decreasing = FALSE),n))
  }
  genes <- unique(genes)
  return(rownames(effect_matrix)[genes])
}
