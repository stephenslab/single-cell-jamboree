# Script to prepare the pancreas LSA data (Stancill et al, 2021) for
# analysis with NMF methods (including the topic model). See
# pancreas_cytokine_lsa_clustering.Rmd for background on these data.
#
# These data preparation steps are slightly different than what is in
# pancreas_cytokine_lsa_clustering.Rmd, and therefore I save the final
# output as "pancreas_cytokine_lsa_v2.RData".
#
library(Matrix)
library(tools)
library(readr)
library(dplyr)

# Set the seed for reproducibility.
set.seed(1)

# Import the sample information.
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
barcodes$condition <- factor(barcodes$condition)

# Import the gene information.
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

# Now import in the read counts.
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

# Clean up the working environment.
rm(counts1,counts2,counts3,counts4)
rm(genes1,genes2,genes3,genes4)
rm(barcodes1,barcodes2,barcodes3,barcodes4)

# Remove genes that are not expressed in any cell.
x      <- colSums(counts)
j      <- which(x > 0)
genes  <- genes[j,]
counts <- counts[,j]

# Remove cells in terms (i) total UMI count and (ii) very few
# expressed genes.
x        <- rowSums(counts > 0)
i        <- which(x > 1250)
barcodes <- barcodes[i,]
counts   <- counts[i,]
x        <- rowSums(counts)
i        <- which(x <= 60000)
barcodes <- barcodes[i,]
counts   <- counts[i,]

# Remove cells with high mitochondrial counts.
mito_genes <- which(substr(genes$symbol,1,2) == "mt")
s          <- rowSums(counts)
s_mito     <- rowSums(counts[,mito_genes])
prop_mito  <- s_mito/s
i          <- which(prop_mito < 0.15)
barcodes   <- barcodes[i,]
counts     <- counts[i,]

# Remove MALAT1, ribosomal genes, and mitochondrial genes.
j <- which(!(grepl("^mt-",genes$symbol) | grepl("^Rp[sl]",genes$symbol)))
genes  <- genes[j,]
counts <- counts[,j]
j      <- which(genes$symbol != "Malat1")
genes  <- genes[j,]
counts <- counts[,j]

# Remove again genes that not expressed in any cells.
x      <- colSums(counts > 0)
j      <- which(x > 0)
genes  <- genes[j,]
counts <- counts[,j]

# Save the processed data to an .RData file.
barcodes <- as.data.frame(barcodes)
genes    <- as.data.frame(genes)
save(list = c("barcodes","genes","counts"),
     file = "pancreas_cytokine_lsa_v2.RData")
resaveRdaFiles("pancreas_cytokine_lsa_v2.RData")
