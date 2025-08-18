# Analyze the pancreas data using flashier and fastTopics.
#
# sinteractive --mem=16G -c 8 --time=24:00:00
# module load R/4.2.0
# .libPaths()[1]
# /home/pcarbo/R_libs_4_20
library(tools)
library(Matrix)
library(fastTopics)
library(flashier)
load("../data/yeast.RData")
set.seed(1)

# (1) Fit a Poisson NMF using fastTopics.
# TO DO.

# (2) Fit an NMF using flashier.
# TO DO.

# Save the model fits to an .Rdata file.
# TO DO.
