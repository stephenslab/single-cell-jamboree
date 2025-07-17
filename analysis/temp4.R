set.seed(1)
x  <- rpois(1e7,1/n)
s1 <- sd(log(x + 1))
