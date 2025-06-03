de <- de_analysis(tm,counts,pseudocount = 0.1,lfc.stat = "k3",
                  control = list(ns = 1e4,nc = 4))
library(ggrepel)
pdat <- data.frame(gene = genes$Symbol,
                   lfc  = de$postmean[,"k2"],
				   stringsAsFactors = FALSE)
pdat$rank <- rank(-pdat$lfc)
j <- which(with(pdat,
                lfc > quantile(lfc,0.001) &
                lfc < quantile(lfc,0.999)))
pdat[j,"gene"] <- ""
p1 <- ggplot(pdat,aes(x = rank,y = lfc,label = gene)) +
  geom_point(color = "dodgerblue") +
  geom_text_repel(color = "black",size = 2.25,fontface = "italic",
                  segment.color = "black",segment.size = 0.25,
                  min.segment.length = 0,max.overlaps = Inf) +
  scale_x_continuous(limits = range(pdat$rank)) +
  labs(x = "rank",y = "posterior mean LFC") +
  theme_cowplot(font_size = 10)
print(p1)
