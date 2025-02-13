library(NNLM)
load("../data/newsgroups.RData")
set.seed(1)

# Remove words that appear in fewer than 10 documents.
x      <- colSums(counts > 0)
j      <- which(x > 9)
counts <- counts[,j]

# Remove documents with fewer than 40 words.
s <- rowSums(counts)
i <- which(s > 39)
counts <- counts[i,]
topics <- topics[i]

# Compute the shifted log counts.
a <- 1
s <- rowSums(counts)
s <- s/mean(s)
Y <- MatrixExtra::mapSparse(counts/(a*s),log1p)

# Set a lower bound on the variances.
n  <- nrow(counts)
x  <- rpois(1e7,1/n)
s1 <- sd(log(x + 1))

# Try running flashier NMF again.
fl0 <- flash(Y,ebnm_fn = ebnm_point_exponential,var_type = 0,
             greedy_Kmax = 10,nullcheck = FALSE,backfit = FALSE,
             verbose = 3)
fl_nmf <- flash_init(Y,var_type = 2,S = s1)
fl_nmf <- flash_factors_init(fl_nmf,fl0,ebnm_point_exponential)
fl_nmf <- flash_backfit(fl_nmf,extrapolate = FALSE,maxiter = 100,verbose = 3)
fl_nmf_ldf <- ldf(fl_nmf,type = "i")

stop()

normalize.rows <- function (A)
  A / rowSums(A)
L <- fl_nmf_ldf$L
k <- ncol(L)
colnames(L) <- paste0("k",1:k)
k <- paste0("k",c(1,2,3,4,5,6,7,8,9,10))
L <- L[,k]
# L <- normalize.rows(L)
p1 <- structure_plot(L,grouping = topics,gap = 10)
