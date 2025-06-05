k <- 14
tm <- fit_poisson_nmf(counts,k = k,init.method = "random",method = "em",
                      numiter = 20,verbose = "none",
                      control = list(numiter = 4,nc = 8,extrapolate = FALSE))
tm <- fit_poisson_nmf(counts,fit0 = tm,method = "scd",numiter = 40,
                      control = list(numiter = 4,nc = 8,extrapolate = TRUE),
					  verbose = "none")
