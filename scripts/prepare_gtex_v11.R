# Before running this script, download the following files from the
# GTEx Portal, gtexportal.org (under "Downloads"), and store them in
# the "data" subdirectory:
#
#   GTEx_Analysis_v11_Annotations_SampleAttributesDD.xlsx
#   GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt
#   GTEx_Analysis_2026-05-19_v11_RNASeQCv2.4.3_gene_reads.gct.gz
#
library(tools)
library(data.table)
library(pathways)
library(fastTopics)

# Prepare the sample annotations.
sample_info <- 
    fread("../data/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt.gz",
          sep = "\t",header = TRUE)
class(sample_info) <- "data.frame"
cols <- c("SAMPID","SMCENTER","SMPTHNTS","SMTS","SMTSD","SMUBRID",
          "SMNABTCH","SMNABTCHT",
          "SMGEBTCHT","ANALYTE_TYPE","SMAFRZE")
sample_info <- sample_info[cols]
sample_info <- transform(sample_info,
                         SMCENTER     = factor(SMCENTER),
                         SMTS         = factor(SMTS),
                         SMTSD        = factor(SMTSD),
                         SMUBRID      = factor(SMUBRID),
                         SMNABTCH     = factor(SMNABTCH),
                         SMNABTCHT    = factor(SMNABTCHT),
                         SMGEBTCHT    = factor(SMGEBTCHT),
                         ANALYTE_TYPE = factor(ANALYTE_TYPE),
                         SMAFRZE      = factor(SMAFRZE))
                         
# Prepare the RNA-seq count data.
counts_file <- 
  "../data/GTEx_Analysis_2026-05-19_v11_RNASeQCv2.4.3_gene_reads.gct.gz"
counts <- fread(counts_file,header = TRUE,sep = "\t",skip = 2)
class(counts) <- "data.frame"
gene_info <- counts[c("Name","Description")]
counts <- counts[-(1:2)]
counts <- as.matrix(counts)
counts <- t(counts)
colnames(counts) <- gene_info$Name

# Align the RNA-seq count data to the sample annotations.
rownames(sample_info) <- sample_info$SAMPID
ids <- rownames(counts)
sample_info <- sample_info[ids,]

# Retain the samples used in the RNA-seq analysis.
i <- which(sample_info$SMAFRZE == "RNASEQ")
sample_info <- sample_info[i,]
counts      <- counts[i,]

# Keep only the genes in the "gene_info" database.
data(gene_sets_human)
ids <- sapply(strsplit(gene_info$Name,".",fixed = TRUE),"[[",1)
colnames(counts) <- ids
j <- is.element(ids,gene_sets_human$gene_info$Ensembl)
counts <- counts[,j]
gc(verbose = TRUE)
storage.mode(counts) <- "double"
gene_info <- gene_sets_human$gene_info
ids       <- colnames(counts)
j         <- match(ids,gene_info$Ensembl)
gene_info <- gene_info[j,]
colnames(counts) <-  gene_info$Symbol

# Remove unexpressed genes.
j         <- which(colSums(counts) > 0)
gene_info <- gene_info[j,]
counts    <- counts[,j]

# Remove all ribosomal protein genes.
j <- which(!grepl("^RP[SL]",gene_info$Symbol))
gene_info <- gene_info[j,]
counts    <- counts[,j]

# Remove MALAT1, which may show variation mainly for technical reasons:
j         <- which(gene_info$Symbol != "MALAT1")
gene_info <- gene_info[j,]
counts    <- counts[,j]

# This is the dimension of the counts matrix after completing these
# data preparation steps.
print(dim(counts))

# Sanity checks.
rownames(sample_info) <- NULL
rownames(gene_info) <- NULL
print(all(sample_info$SAMPID == rownames(counts)))
print(all(gene_info$Symbol == colnames(counts)))

# Save the prepared data to an .RData file.
save(list = c("sample_info","gene_info"),
     file = "gtex_meta_data_v11.RData")
save(list = c("sample_info","gene_info","counts"),
     file = "gtex_v11.RData")
resaveRdaFiles("gtex_meta_data_v11.RData")
resaveRdaFiles("gtex_v11.RData")
