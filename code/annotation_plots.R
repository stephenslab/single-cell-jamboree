# This function creates an "effect plot". The effects_matrix input
# should be a matrix in which rows are to features and columns are
# dimension (e.g., factors in a matrix factorization). Effects
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
         scale_fill_manual(values = c("darkblue","gray","tomato"),
                           drop = FALSE) +
         guides(size = guide_legend(override.aes = list(fill = "black")),
                fill = guide_legend(override.aes = list(size = 3))) +
         labs(x = "dimension",y = "feature name",
              fill = "effect sign",size = "effect size") +
         theme_cowplot(font_size = font_size))
}
