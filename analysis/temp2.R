# Some code to be incorporated into lps.Rmd.
library(rsvd)
library(uwot)
set.seed(1)
U <- rsvd(shifted_log_counts,k = 40)$u
Y <- umap(U,n_neighbors = 20,metric = "cosine",min_dist = 0.3,
          n_threads = 8,verbose = TRUE)
x <- Y[,1]
y <- Y[,2]
samples$umap1 <- x
samples$umap2 <- y

# UMAP plot by (1) tissue and (2) factor 6.
tissue_colors <- c("magenta","darkorange","darkblue","forestgreen",
                   "dodgerblue","red","olivedrab","darkmagenta",
                   "sienna","limegreen","royalblue","lightskyblue",
                   "gold")
L <- ldf(fl_nmf,type = "i")$L
colnames(L) <- paste0("k",1:15)
pdat <- cbind(data.frame(k6 = L[,"k6"]),samples)
p1 <- ggplot(pdat,aes(x = umap1,y = umap2,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  theme_cowplot(font_size = 10)
p2 <- ggplot(pdat,aes(x = umap1,y = umap2,color = k6)) +
  geom_point(size = 1) +
  scale_color_gradient2(low = "deepskyblue",mid = "gold",high = "tomato",
                        midpoint = 0.66) +
  theme_cowplot(font_size = 10)
plot_grid(p1,p2,nrow = 1,ncol = 2)
