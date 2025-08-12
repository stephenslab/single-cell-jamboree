# Some code to be added to the lps.Rmd workflowr analysis.

# --- structure-plot-flashier-nmf-k9 ---
rows <- order(samples$timepoint)
topic_colors <- c("powderblue","dodgerblue","olivedrab","limegreen",
                  "forestgreen","red","darkmagenta","gray","darkorange",
                  "cyan","royalblue","darkblue","lightskyblue",
                  "gold","sienna")
L <- ldf(fl_nmf_k9,type = "i")$L
structure_plot(L,grouping = samples$tissue,gap = 4,
               topics = 1:9,colors = topic_colors,
               loadings_order = rows) +
  labs(fill = "",y = "membership")

# --- structure-plot-flashier-nmf ---
L <- ldf(fl_nmf,type = "i")$L
p <- structure_plot(L,grouping = samples$tissue,gap = 4,
               topics = 1:15,colors = topic_colors,
			   loadings_order = rows) +
  labs(fill = "",y = "membership") +
  theme(axis.text.x = element_text(angle = 0,hjust = 0.5),
        legend.key.height = unit(0.01,"cm"),
        legend.key.width = unit(0.25,"cm"),
        legend.text = element_text(size = 7))
