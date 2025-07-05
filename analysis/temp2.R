F1 <- ldf(fl_nmf,type = "i")$F
rownames(F1) <- genes$symbol
colnames(F1) <- paste0("k",1:13)
F2 <- poisson2multinom(tm)$F
plot(F1[,"k1"],log10(F2[,"k2"] + 1e-5),pch = 20)
plot(F1[,"k5"],log10(F2[,"k3"] + 1e-5),pch = 20)
plot(F1[,"k13"],log10(F2[,"k8"] + 1e-5),pch = 20)
plot(F1[,"k3"],log10(F2[,"k5"] + 1e-5),pch = 20)
plot(F1[,"k4"],log10(F2[,"k7"] + 1e-5),pch = 20)
plot(F1[,"k6"],log10(F2[,"k6"] + 1e-5),pch = 20)
