# TO DO: Figure out way to scale the F matrix so that they are
# comparable.

F <- ldf(fl_nmf,type = "i")$F
colnames(F) <- paste0("k",1:4)
pdat <- cbind(genes[c("GeneID","Symbol")],F)
p1 <- ggplot(pdat,aes(x = k2,y = k3)) +
  geom_point() +
  geom_abline(intercept = 0,slope = 1,color = "magenta",
              linetype = "dotted") +
  theme_cowplot(font_size = 12)
p2 <- ggplot(pdat,aes(x = k2,y = k4)) +
  geom_point() +
  geom_abline(intercept = 0,slope = 1,color = "magenta",
              linetype = "dotted") +
  theme_cowplot(font_size = 12)
p3 <- ggplot(pdat,aes(x = k3,y = k4)) +
  geom_point() +
  geom_abline(intercept = 0,slope = 1,color = "magenta",
              linetype = "dotted") +
  theme_cowplot(font_size = 12)
print(plot_grid(p1,p2,p3,nrow = 1,ncol = 3))

