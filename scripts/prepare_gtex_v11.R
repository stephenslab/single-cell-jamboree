# Before running this script, download the following files from the
# GTEx Portal, gtexportal.org (under "Downloads"), and store them in
# the "data" subdirectory:
#
#   GTEx_Analysis_v11_Annotations_SampleAttributesDD.xlsx
#   GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt
#   GTEx_Analysis_2026-05-19_v11_RNASeQCv2.4.3_gene_reads.gct.gz
#
library(data.table)
library(fastTopics)
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
                         
# Save the prepared data to an .RData file.
# TO DO.

# Fit a topic model to the data.
# TO DO.
