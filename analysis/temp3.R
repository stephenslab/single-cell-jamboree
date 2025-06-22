# NOTE: Relative sizes of the F estimates can really matter.

# plot(de_vsnull$postmean[,"k8"],de_vsnull$postmean[,"k9"],pch = 20)
# plot(de_vsnull$postmean[,"k9"],de_le$postmean[,"k9"],pch = 20)
# p1 <- volcano_plot(de_vsnull,k = "k8",ymax = 50,labels = genes$symbol)
# p2 <- volcano_plot(de_le,k = "k8",ymax = 50,labels = genes$symbol)

F <- ldf(fl_nmf,type = "i")$F
rownames(F) <- genes$symbol
colnames(F) <- paste0("k",1:9)
F_distinct <- rank_transform_effects_matrix(F,compare_cols = TRUE,
                                            compare_dims = celltype_factors)
plot(F[,6],F_distinct[,6],pch = 20)

#
# TO DO NEXT:
#
# (1) Implement function compute_distinctive_changes in
#     singlecelljamboreeR.
#
# (2) Develop a function in the singlecelljamboreeR package to
#    create plots of "changes" vs. "distinctive changes".
#

