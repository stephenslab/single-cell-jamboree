F1 <- ldf(fl_nmf,type = "i")$F
F2 <- poisson2multinom(tm)$F
rownames(F1) <- genes$symbol
colnames(F1) <- paste0("k",1:13)
rownames(F2) <- genes$symbol
plot(F1[,"k1"],log10(F2[,"k2"] + 1e-5),pch = 20)
plot(F1[,"k5"],log10(F2[,"k3"] + 1e-5),pch = 20)
plot(F1[,"k13"],log10(F2[,"k8"] + 1e-5),pch = 20)
plot(F1[,"k3"],log10(F2[,"k5"] + 1e-5),pch = 20)
plot(F1[,"k4"],log10(F2[,"k7"] + 1e-5),pch = 20)
plot(F1[,"k6"],log10(F2[,"k6"] + 1e-5),pch = 20)
plot(F1[,"k12"],F[,"k1"],pch = 20)
points(F1[cell_cycle,"k12"],F[cell_cycle,"k1"],pch = 20,col = "orangered")
plot(F1[,"k7"],F[,"k1"],pch = 20)
plot(F1[,"k12"],F[,"k13"],pch = 20)
plot(F1[,"k1"],F[,"k4"],pch = 20)
cor(F1[,"k12"],F[,"k1"])
cor(F1[,"k12"],F[,"k13"])
cor(F1[,"k1"],F[,"k13"])
points(F1[cell_cycle,"k7"],F[cell_cycle,"k1"],pch = 20,col = "cyan")
plot(F1[,"k7"],F[,"k13"],pch = 20)
plot(F1[,"k1"],F[,"k4"],pch = 20)


