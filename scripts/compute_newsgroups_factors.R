# TO DO: Explain here what this script is for, and how to use it.
#
# sinteractive --mem=16G -c 8 --time=24:00:00
# module load R/4.2.0
# .libPaths()[1]
# /home/pcarbo/R_libs_4_20
library(tools)
library(Matrix)
library(NNLM)
library(fastTopics)
library(flashier)
load("../data/pancreas.RData")
set.seed(1)
