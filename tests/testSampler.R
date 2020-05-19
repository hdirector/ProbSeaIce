#Test script shows that MCMC works with large number of samples and drawing
#from the truth

#Assess coverage under particular conditions
rm(list = ls())
set.seed(103)
library("IceCast")

#true  settings
n_lines <- 20
mu_true <- rep(logit(.5), n_lines)
angs <- seq(0, 2*pi, length = n_lines + 1)
angs <- angs[1:n_lines]
dists <- theta_dist_mat(thetas = angs)
kappa_true <- 3
sigma_true <- rep(1, n_lines)
Sigma_true <- compSigma(sigma = sigma_true, kappa = kappa_true, dists = dists)

#priors
mu0 <- rep(logit(.48), n_lines)
Lambda0 <- 3.5^2*diag(n_lines)
alpha_kappa0 <- .1
beta_kappa0 <- 8
beta_sigma0 <- rep(2.5, n_lines)

#MCMC parameters
burn_in <- 5000
n_iter <- 10000
mu_prop_sd <- .15
sigma_prop_sd <- .1
kappa_prop_sd <-.05

#generate data from truth
n_obs <- 200
prop_tilde <- t(MASS::mvrnorm(n = n_obs, mu = mu_true, Sigma = Sigma_true))

#initial values
kappa_ini <-  2
sigma_min = .01
mu_ini <- apply(prop_tilde, 1, mean)
emp_cov_y <- cov(t(prop_tilde))
sigma_ini <- sqrt(diag(emp_cov_y))
sigma_high <- which(sigma_ini >= beta_sigma0)
sigma_ini[sigma_high] <- beta_sigma0[sigma_high] - .2
sigma_ini[which(sqrt(diag(emp_cov_y)) == 0)] <- sigma_min


#test MCMC
fits <- RunMCMC(nIter = n_iter, prop_tilde = prop_tilde,
               mu = mu_ini, mu0 = mu0, Lambda0 = Lambda0, muPropSD = mu_prop_sd,
               kappa = kappa_ini, alphaKappa0 = alpha_kappa0, betaKappa0 = beta_kappa0,
               kappaPropSD = kappa_prop_sd, sigma = sigma_ini,
               betaSigma0  = beta_sigma0, sigmaPropSD = sigma_prop_sd,
               dists = dists)


#parameter estimates
mu_est <- apply(fits$mu[,(burn_in + 1):n_iter], 1, mean)
kappa_est <- mean(fits$kappa[(burn_in + 1):n_iter])
sigma_est <- apply(fits$sigma[,(burn_in + 1):n_iter],1, mean)

par(mfrow = c(2, 2))
Sigma_obs <- cov(t(prop_tilde))
Sigma_ini <- compSigma(sigma_ini, kappa_ini, dists)
Sigma_est <- compSigma(sigma_est, kappa_est, dists)
library("fields")
image.plot(Sigma_est, zlim  = c(0, 1.6), main =  "est")
image.plot(Sigma_obs, zlim = c(0, 1.6), main = "obs")
image.plot(Sigma_true, zlim = c(0, 1.6), main = "true")
image.plot(Sigma_est - Sigma_true, main = "diff")

par(mfrow = c(1, 1))
plot(sigma_est, type= "l", ylim = c(.8, 1.2))
points(sigma_ini, type= "l", col = "blue")
points(sigma_true, type= "l", col = 'red')
points(beta_sigma0, type = 'l', col = 'purple')
legend("bottom", fill = c("black", "blue", "red", "purple"),
       legend = c("est", "ini", "true", "priot"),
       ncol = 2, cex = .75)

plot(mu_est, type= "l", ylim = c(-.5, .15))
points(mu_ini, col= 'blue', type= "l")
points(mu_true, type= 'l', col = 'red')
points(mu0, type= 'l', col = 'purple')
legend("bottom", fill = c("black", "blue", "red", "purple"),
       legend = c("est", "ini", "true", "priot"),
       ncol = 2, cex = .75)

plot(fits$kappa, type= "l")
abline(h = kappa_true, col = 'red')
kappa_est; kappa_true


