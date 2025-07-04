# Run this script after running the code in lps.Rmd.
library(Matrix)
library(pathways)
library(readr)
library(singlecelljamboreeR)
F <- poisson2multinom(tm)$F[,c("k1","k7")]
F <- scale(F,center = FALSE,scale = TRUE)
set.seed(1)
data(gene_sets_mouse)
gene_sets     <- gene_sets_mouse$gene_sets
gene_info     <- gene_sets_mouse$gene_info
gene_set_info <- gene_sets_mouse$gene_set_info
j <- which(with(gene_sets_mouse$gene_set_info,
                (database == "MSigDB-C2" &
                 grepl("CP",sub_category_code,fixed = TRUE)) |
                (database == "MSigDB-C5") &
                 grepl("GO",sub_category_code,fixed = TRUE)))
genes <- sort(intersect(rownames(F),gene_info$Symbol))
i <- which(is.element(gene_info$Symbol,genes))
gene_info     <- gene_info[i,]
gene_set_info <- gene_set_info[j,]
gene_sets     <- gene_sets[i,j]
rownames(gene_sets) <- gene_info$Symbol
rownames(gene_set_info) <- gene_set_info$id
gene_set_info <- gene_set_info[,-2]
gsea_tm <- singlecelljamboreeR::perform_gsea(F,gene_sets,gene_set_info,
                                             L = 36)
F <- ldf(fl_nmf,type = "i")$F
F <- scale(F,center = FALSE,scale = TRUE)
gsea_fl_nmf <- singlecelljamboreeR::perform_gsea(F[,6],gene_sets,
                                                 gene_set_info,L = 36)
out <- gsea$selected_gene_sets
out$top_genes <- sapply(out$top_genes,function (x) paste(x,collapse = " "))
write_csv(out,"lps_gsea.csv",quote = "none")

F <- ldf(fl_nmf,type = "i")$F
colnames(F) <- paste0("k",1:15)
pdat <- data.frame(tm   = poisson2multinom(tm)$F[,"k7"],
                   nmf  = F[,"k6"],
                   gene = rownames(F),
                   pathway = FALSE)
pathway_genes <- names(which(gene_sets[,"M16779"] > 0))
pdat[pathway_genes,"pathway"] <- TRUE
rows <- which(with(pdat,tm < 0.005 & nmf < 0.8))
pdat[rows,"gene"] <- ""
ggplot(pdat,aes(x = tm,y = nmf,label = gene,color = pathway)) +
  geom_point() +
  geom_text_repel(color = "gray",size = 2.5,max.overlaps = Inf) +
  scale_x_continuous(trans = "sqrt") +
  scale_color_manual(values = c("darkblue","darkorange")) +
  labs(x = "topic 7",y = "factor 6") +
  theme_cowplot(font_size = 12)
