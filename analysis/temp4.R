set.seed(1)
n  <- nrow(shifted_log_counts)
x  <- rpois(1e7,1/n)
s1 <- sd(log(x + 1))
set.seed(1)
fl_nmf_new <- 
  singlecelljamboreeR::flashier_nmf(shifted_log_counts,k = 13,n.threads = 8,
                                    var_type = 2,S = s1)
