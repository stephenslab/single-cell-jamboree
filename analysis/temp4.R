# Some code to be added to the lps.Rmd workflowr analysis.

# --- flashier-nmf-k9 ---
set.seed(1)
n  <- nrow(shifted_log_counts)
x  <- rpois(1e7,1/n)
s1 <- sd(log(x + 1))
set.seed(1)
fl_nmf_k9 <- flashier_nmf(shifted_log_counts,k = 9,greedy_init = TRUE,
                          var_type = 2,S = s1)

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

# --- flashier-nmf-k15 ---
set.seed(1)
fl_nmf <- flashier_nmf(shifted_log_counts,k = 15,greedy_init = TRUE,
                       var_type = 2,S = s1)

# --- structure-plot-flashier-nmf ---
#
# TO DO
#
