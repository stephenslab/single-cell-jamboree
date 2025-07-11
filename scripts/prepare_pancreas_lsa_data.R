## ----knitr-opts, include=FALSE------------------------------------------------
knitr::opts_chunk$set(comment = "#",collapse = TRUE,results = "hold",
                      fig.align = "center",dpi = 120)


## ----load-pkgs, message=FALSE-------------------------------------------------
library(Matrix)
library(tools)
library(readr)
library(fastTopics)
library(dplyr)
library(Seurat)


## ----set-seed-----------------------------------------------------------------
set.seed(1)


## ----import-barcodes, message=FALSE-------------------------------------------
barcodes1<-read_tsv("../data/GSE156175_RAW/GSM5842388_Rep1_S1_barcodes.tsv.gz",
                    col_names = "barcode")
barcodes2 <- read_tsv("../data/GSE156175_RAW/GSM4726017_S2_barcodes.tsv.gz",
                      col_names = "barcode")
barcodes3 <- read_tsv("../data/GSE156175_RAW/GSM4726018_S3_barcodes.tsv.gz",
                      col_names = "barcode")
barcodes4 <- read_tsv("../data/GSE156175_RAW/GSM4726019_S4_barcodes.tsv.gz",
                      col_names = "barcode")
barcodes1$barcode <- substr(barcodes1$barcode,1,nchar(barcodes1$barcode) - 2)
barcodes2$barcode <- substr(barcodes2$barcode,1,nchar(barcodes2$barcode) - 2)
barcodes3$barcode <- substr(barcodes3$barcode,1,nchar(barcodes3$barcode) - 2)
barcodes4$barcode <- substr(barcodes4$barcode,1,nchar(barcodes4$barcode) - 2)
barcodes1$condition <- "Untreated"
barcodes2$condition <- "IL-1B"
barcodes3$condition <- "IFNg"
barcodes4$condition <- "IL-1B_IFNg"
barcodes <- rbind(barcodes1, barcodes2, barcodes3, barcodes4)
barcodes$cell_bc <- paste0(barcodes$barcode, "_", barcodes$condition)
barcodes$condition <- factor(barcodes$condition)


## ----import-gene-info, message=FALSE------------------------------------------
genes1 <- read_tsv("../data/GSE156175_RAW/GSM5842388_Rep1_S1_features.tsv.gz",
                   col_names = c("ensembl", "symbol", "type"))
genes2 <- read_tsv("../data/GSE156175_RAW/GSM4726017_S2_features.tsv.gz",
                   col_names = c("ensembl", "symbol", "type"))
genes3 <- read_tsv("../data/GSE156175_RAW/GSM4726018_S3_features.tsv.gz",
                   col_names = c("ensembl", "symbol", "type"))
genes4 <- read_tsv("../data/GSE156175_RAW/GSM4726019_S4_features.tsv.gz",
                   col_names = c("ensembl", "symbol", "type"))
genes <- genes1
genes$type <- factor(genes$type)


## ----import-counts, message=FALSE---------------------------------------------
counts1 <- t(readMM("../data/GSE156175_RAW/GSM5842388_Rep1_S1_matrix.mtx.gz"))
counts2 <- t(readMM("../data/GSE156175_RAW/GSM4726017_S2_matrix.mtx.gz"))
counts3 <- t(readMM("../data/GSE156175_RAW/GSM4726018_S3_matrix.mtx.gz"))
counts4 <- t(readMM("../data/GSE156175_RAW/GSM4726019_S4_matrix.mtx.gz"))
rownames(counts1) <- paste0(barcodes1$barcode, "_", barcodes1$condition)
rownames(counts2) <- paste0(barcodes2$barcode, "_", barcodes2$condition)
rownames(counts3) <- paste0(barcodes3$barcode, "_", barcodes3$condition)
rownames(counts4) <- paste0(barcodes4$barcode, "_", barcodes4$condition)
colnames(counts1) <- genes1$symbol
colnames(counts2) <- genes2$symbol
colnames(counts3) <- genes3$symbol
colnames(counts4) <- genes4$symbol
counts <- rbind(
  as.matrix(counts1),
  as.matrix(counts2),
  as.matrix(counts3),
  as.matrix(counts4)
)
counts <- as(counts, "CsparseMatrix")


## ----filter-genes-1-----------------------------------------------------------
x      <- colSums(counts)
j      <- which(x > 0)
genes  <- genes[j,]
counts <- counts[,j]


## ----filter-cells-------------------------------------------------------------
x        <- rowSums(counts > 0)
i        <- which(x > 2000)
barcodes <- barcodes[i,]
counts   <- counts[i,]
x        <- rowSums(counts)
i        <- which(x <= 60000)
barcodes <- barcodes[i,]
counts   <- counts[i,]


## ----filter-cells-2-----------------------------------------------------------
mito_genes <- which(substr(genes$symbol,1,2) == "mt")
s          <- rowSums(counts)
s_mito     <- rowSums(counts[,mito_genes])
prop_mito  <- s_mito/s
i          <- which(prop_mito < 0.1)
barcodes   <- barcodes[i,]
counts     <- counts[i,]


## ----filter-genes-2-----------------------------------------------------------
j <- which(!(grepl("^mt-",genes$symbol) | grepl("^Rp[sl]",genes$symbol)))
genes  <- genes[j,]
counts <- counts[,j]
j      <- which(genes$symbol != "Malat1")
genes  <- genes[j,]
counts <- counts[,j]


## ----filter-genes-3-----------------------------------------------------------
x      <- colSums(counts > 0)
j      <- which(x > 2)
genes  <- genes[j,]
counts <- counts[,j]


## ----write-outputs, eval=FALSE------------------------------------------------
# save(list = c("barcodes","genes","counts"),
#      file = "pancreas_cytokine_lsa_v2.RData")
# resaveRdaFiles("pancreas_cytokine_lsa_v2.RData")


## ----fit-topic-model, cache=TRUE----------------------------------------------
set.seed(1)
tm_fit <- fit_poisson_nmf(X = counts,k = 12,control = list(nc = 7),
                          verbose = "none")


## ----structure-plot, fig.height=1.75, fig.width=7, results="hide"-------------
set.seed(1)
structure_plot(tm_fit,n = Inf)


## ----annotation-heatmap, fig.height=2.5, fig.width=4--------------------------
scale_rows <- function (A)
  A / apply(A,1,max)
marker_genes <- c("Ins1","Ins2","Gcg","Sst","Ppy","Krt17","Ccr5",
                  "Col1a1","Cd34","Cpa1")
F <- poisson2multinom(tm_fit)$F
F <- F[marker_genes,]
F <- scale_rows(F)
rownames(F) <- marker_genes
annotation_heatmap(F,select_features = "all",verbose = FALSE)


## ----clustering, fig.height=2.5, fig.width=7, results="hide"------------------
L <- poisson2multinom(tm_fit)$L
L_df <- as.data.frame(L)
L_df <- L_df %>%
  mutate(
    cluster = case_when(
      k2 > (1/3) ~ "Acinar",
      k3 > (1/3) ~ "Endothelial/Mesnchymal",
      k5 > (1/3) ~ "Ductal",
      k8 > (1/3) ~ "Delta",
      k9 > (1/3) ~ "Alpha",
      k10 > (1/3) ~ "Gamma",
      k11 > (1/3) ~ "Macrophage",
      TRUE ~ "Beta"
    )
  )
structure_plot(tm_fit, grouping = L_df$cluster, gap = 25, n = Inf)


## ----elbo-plot, fig.height=3, fig.width=3-------------------------------------
barcodes$celltype <- L_df$cluster
counts <- counts[,!duplicated(colnames(counts))]
so <- CreateSeuratObject(counts = Matrix::t(counts), meta.data = barcodes)
so <- SCTransform(so, verbose = FALSE)
so <- RunPCA(so, verbose = FALSE)
ElbowPlot(so)


## ----umap-plot, fig.height=3, fig.width=5-------------------------------------
so <- RunUMAP(so, dims = 1:10, verbose = FALSE)
DimPlot(so, group.by = "celltype") +
  theme_cowplot(font_size = 10)

