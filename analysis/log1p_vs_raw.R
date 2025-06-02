# TO DO: Convert this script to R Markdown.

# First set the seed so that the results are reproducible.
set.seed(1)

# This is the number of RNA-seq samples to simulate.
n <- 200

# Simulate the size factors.
s <- pmax(10,1000 * rnorm(n,10,4))

# Simulate the membership matrix.
k <- 9
L <- matrix(0,n,k)
for (i in 1:n) {
  # TO DO.
}

# Simulate the gene matrix.
# TO DO.

# Simulate the counts.
counts <- rpois(n*m,tcrossprod(s*L,F))
counts <- as.double(counts)
counts <- matrix(counts,n,m)
counts <- generate_poisson_nmf_counts(F,s*L)
