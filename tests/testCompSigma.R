#test compSigma function that runs in C++
library("IceCast")
library("fields")

###########################
#Check covariance functions
###########################
n_lines <- 10
sigma <- c(rep(1, n_lines/2), rep(2, n_lines/2))
kappa <- 5
dists <- dist_mat(n_lines)

Sigma_fn <- compSigma(sigma = sigma, kappa = kappa, dists = dists)
image.plot(Sigma_fn)

#compute by hand
Sigma_hand <- matrix(nrow = n_lines, ncol = n_lines)
for (i in 1:n_lines) {
  for (j in 1:n_lines) {
    Sigma_hand[i, j] <- sigma[i]*sigma[j]*exp(-dists[i,j]/kappa)
  }
}
image.plot(Sigma_hand)

sum(Sigma_hand == Sigma_fn)
sum(Sigma_hand != Sigma_fn)
