# gtex_v11.RData can be downloaded from this shared Box folder: 
# https://uchicago.box.com/s/l0emrjcqpw5yat1zaygtr1v3akq7de1t
#
# sinteractive -p mstephens --mem=24G -c 20 --nodelist=midway2-0440 \
#   --time=48:00:00
# module load R/4.2.0
# .libPaths()[1]
# /home/pcarbo/R_libs_4_20
# 
library(tools)
library(fastTopics) # 0.7.46
k <- 31
outfile <- sprintf("gtex_tm_k=%d.RData",k)
cat("k =",k,"\n")
cat("outfile =",outfile,"\n")
load("../data/gtex_v11.RData")
set.seed(1)

# Fit a Poisson NMF model using fastTopics.
# This step is expected to take X h.
t0 <- proc.time()
tm <- fit_poisson_nmf(counts,k = k,numiter = 100,method = "em",
                      control = list(numiter = 4,nc = 20,extrapolate = FALSE),
                      init.method = "random",verbose = "detailed")
tm <- fit_poisson_nmf(counts,fit0 = tm,numiter = 100,method = "scd",
                      control = list(numiter = 4,nc = 20,extrapolate = TRUE),
                      verbose = "detailed")
t1 <- proc.time()
print(t1 - t0)

# Save the model fits to an .Rdata file.
session_info <- sessionInfo()
save(list = c("tm","session_info"),file = outfile)
resaveRdaFiles(outfile)
