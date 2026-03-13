# Script to prepare the combined GCA (tissue + organoid) data.
#
# sinteractive -c 2 --mem=36G --time=24:00:00 -p mstephens \
#   --account=pi-mstephens
# module load R/4.2.0
# > .libPaths()[1]
# "/home/pcarbo/R_libs_4_20"
#
library(tools)
library(dplyr)
library(Matrix)
library(Seurat)
dat <- readRDS("tissue_RNA_integrated.rds")
tissue <- dat[["RNA"]]$counts
meta_tissue <- dat@meta.data
dat <- readRDS("org_RNA_integrated.rds")
organoid <- dat[["RNA"]]$counts
meta_organoid <- dat@meta.data

# Merge the meta data.
names(meta_organoid) <- tolower(names(meta_organoid))
cols <- c("adjusted_lineage","batch","cellannotation","celltype",
          "condition","group","main_label","origin","patient","reannotation",
          "region")
meta_tissue <- meta_tissue[cols]
names(meta_tissue)   <- tolower(names(meta_tissue))
cols <- c("batch","celltype","condition","group","orig.ident","patient",
          "region")
meta_organoid <- meta_organoid[cols]
names(meta_organoid)[5] <- "origin"
meta_organoid <- transform(meta_organoid,celltype = as.character(celltype))
sample_info <- bind_rows(meta_tissue,meta_organoid)
n <- length(sample_info)
for (i in 1:n)
  sample_info[[i]] <- factor(sample_info[[i]])
                         
# Merge the count data.
genes1   <- rownames(tissue)
genes2   <- rownames(organoid)
genes    <- intersect(genes1,genes2)
tissue   <- tissue[genes,]
organoid <- organoid[genes,]
counts   <- cbind(tissue,organoid)
rm(tissue,organoid)
counts <- t(counts)

# Save the combined data to an .RData file.
save(list = c("sample_info","counts"),file = "gca.RData")
resaveRdaFiles("gca.RData")
