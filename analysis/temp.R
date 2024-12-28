# NOTES:
#
# Downloaded GSE132188_adata.h5ad.h5 from
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132188
#
# Ran prepare_pancreas_endocrine_data.py
# to generate pancreas_endocrine_alldays.h5ad.
# 
# $ conda activate base
# $ conda list | grep anndata
# anndata 0.11.1 pyhd8ed1ab_1 conda-forge
#
library(Matrix)
library(anndata)
library(reticulate)
library(tools)
use_python("/Users/pcarbo/miniforge3/bin/python")

# Retrieve the count data prepared using the Python script.
dat1 <- read_h5ad("../data/pancreas_endocrine_alldays.h5ad")
counts <- dat1$X
counts <- as(counts,"CsparseMatrix")

# Get the meta-data downloaded from GEO.
dat2 <- read_h5ad("../data/GSE132188_adata.h5ad.h5")
ids1 <- rownames(dat1$obs)
ids2 <- rownames(dat2$obs)
ids2 <- paste0("e",10*as.numeric(as.character(dat2$obs$day)),"-",ids2)
ids2 <- substr(ids2,1,23)

# Merge the count data with the meta-data downloaded from GEO.
rows   <- which(is.element(ids1,ids2))
ids1   <- ids1[rows]
counts <- counts[rows,]
obs1   <- dat1$obs[rows,]

# Check that the sample ids and genes are the same.
print(all(ids1 == ids2))
print(all(rownames(dat1$var) == rownames(dat2$var)))

# Extract the gene info.
gene_info <- dat2$var
gene_info <- cbind(gene = rownames(gene_info),gene_info)
rownames(gene_info) <- NULL

# Extract the sample info.
sample_info <- dat2$obs
umap <- dat2$obsm$X_umap
colnames(umap) <- c("umap1","umap2")
sample_info <- cbind(umap,sample_info)

# Write the data to an .RData file.
save(list = c("gene_info","sample_info","counts"),
     file = "pancreas_endocrine.RData")
resaveRdaFiles("pancreas_endocrine.RData")
