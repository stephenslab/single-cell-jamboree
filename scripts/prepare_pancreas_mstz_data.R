# Before running this script, download and extract GSE128565_RAW.tar
# from GEO, accession GSE128565. After extracting the tar archives and
# moving files around, you should have the following directory
# structure in ../data/GSE128565:
#
# ├── GSM3680207
# │   ├── barcodes.tsv
# │   ├── genes.tsv
# │   └── matrix.mtx
# ├── GSM3680208
# │   ├── barcodes.tsv
# │   ├── genes.tsv
# │   └── matrix.mtx
# ├── GSM3680209
# │   ├── barcodes.tsv
# │   ├── genes.tsv
# │   └── matrix.mtx
# ├── GSM3680210
# │   ├── barcodes.tsv
# │   ├── genes.tsv
# │   └── matrix.mtx
# ├── GSM3680211
# │   ├── barcodes.tsv
# │   ├── genes.tsv
# │   └── matrix.mtx
# ├── GSM3680212
# │   ├── barcodes.tsv
# │   ├── genes.tsv
# │   └── matrix.mtx
# └── GSM3680213
#     ├── barcodes.tsv
#     ├── genes.tsv
#     └── matrix.mtx
# 
# This is the correspondence between the GEO accessions and the
# treatments/conditions:
#
# GSM3680207  No STZ
# GSM3680208  STZ
# GSM3680209  STZ insulin
# GSM3680210  STZ GLP-1
# GSM3680211  STZ estrogen
# GSM3680212  STZ GLP-1/estrogen
# GSM3680213  STZ GLP-1/estrogen + insulin
# 
library(Matrix)
library(tools)
library(readr)

# Set the seed for reproducibility.
set.seed(1)

# Import the barcodes.
barcodes1 <- read_tsv("../data/GSE128565/GSM3680207/barcodes.tsv",
                      col_names = "barcode")
barcodes2 <- read_tsv("../data/GSE128565/GSM3680208/barcodes.tsv",
                      col_names = "barcode")
barcodes3 <- read_tsv("../data/GSE128565/GSM3680209/barcodes.tsv",
                      col_names = "barcode")
barcodes4 <- read_tsv("../data/GSE128565/GSM3680210/barcodes.tsv",
                      col_names = "barcode")
barcodes5 <- read_tsv("../data/GSE128565/GSM3680211/barcodes.tsv",
                      col_names = "barcode")
barcodes6 <- read_tsv("../data/GSE128565/GSM3680212/barcodes.tsv",
                      col_names = "barcode")
barcodes7 <- read_tsv("../data/GSE128565/GSM3680213/barcodes.tsv",
                      col_names = "barcode")
barcodes1 <- cbind(barcodes1,data.frame(condition = "no_STZ"))
barcodes2 <- cbind(barcodes2,data.frame(condition = "STZ"))
barcodes3 <- cbind(barcodes3,data.frame(condition = "insulin"))
barcodes4 <- cbind(barcodes4,data.frame(condition = "GLP-1"))
barcodes5 <- cbind(barcodes5,data.frame(condition = "estrogen"))
barcodes6 <- cbind(barcodes6,data.frame(condition = "GLP-1-estrogen"))
barcodes7 <- cbind(barcodes7,data.frame(condition = "GLP-1-estrogen+insulin"))
barcodes <- rbind(barcodes1,barcodes2,barcodes3,barcodes4,
                  barcodes5,barcodes6,barcodes7)
barcodes <- 
  transform(barcodes,
            condition = factor(condition,
                               c("no_STZ","STZ","estrogen","GLP-1",
                                 "GLP-1-estrogen","insulin",
                                 "GLP-1-estrogen+insulin")))
