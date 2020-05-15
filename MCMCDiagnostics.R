#Script for MCMC Diagnostics tables and figures in supplement
#Focuses on forecast for September 2005 at 2-month lead

#load MCMC results
load("Results/cont_fits/fit_Month9_Train1995_2004_init7.rda")

#libraries
library("coda") #for raftery.diag()
library("xtable")

#general info
regs_to_fit <- which(sapply(res, function(x){!is.null(x)}))
reg_names <- c("Central Arctic", "Baffin Bay", "Greenland Sea") #need to update manually!
n_reg <- length(regs_to_fit)
n_iter <- length(res[[1]]$kappa)
quants_to_assess <- c(.5, .85, 1) #what quantiles to assess when there is more than on parameter value

#--------------------------
#Assess sigma
#--------------------------
ub_prop <- c(.99, .99, .99, .99, .7263) #need to update manually!
lb_prop <- c(.15, .01, .01, .01, .01) #need to update manually!

N_sigma_rd <- matrix(ncol = length(quants_to_assess) + 1, nrow = n_reg)
M_sigma_rd <- matrix(ncol = length(quants_to_assess) + 1, nrow = n_reg)
colnames(N_sigma_rd) <- c("Region", paste("n_", quants_to_assess, sep = ""))
#png(filename = sprintf("/users/hdirector/Dropbox/SeaIce_InProgress/probContours_ECMWF/Paper/Figures/traceplot_sigma.png",
#                       r), res = 100)
par(mfrow = c(3, 1), oma = rep(1, 4), mar = c(2, 4, 2, 2))
for (j in 1:n_reg) {
  r <- regs_to_fit[j]
  n_lines <- dim(res[[r]]$sigma)[1]
  beta_sigma0 <- ((logit(ub_prop[r]) - logit(lb_prop[r]))/2)/qnorm(.995)
  N_sigma <- rep(NA, n_lines)
  M_sigma <- rep(NA, n_lines)
  for (i in 1:n_lines) {
    #raftery.diag doesn't make sense if chain is right near boundary,
    #if 95% of samples are right by boundary don't assess
    #WANT sigma's at boundary here
    at_bd <- FALSE
    if ((sum(abs(res[[r]]$sigma[i,] - rep(0, n_iter)) <= .05)/n_iter >= .95) ||
          (sum(abs(res[[r]]$sigma[i,] - rep(beta_sigma0, n_iter)) <= .05)/n_iter >= .95)) {
        at_bd <- TRUE
        N_sigma[i] <- NA
        M_sigma[i] <- NA
      } else  {
        #For values away from bound, compute raftery.diag()
        lb <- raftery.diag(res[[r]]$sigma[i,], q = .025, r = .0125)
        ub <- raftery.diag(res[[r]]$sigma[i,], q = .975, r = .0125)
        N_sigma[i] <- max(lb$resmatrix[,"N"], ub$resmatrix[,"N"])
        M_sigma[i] <- max(lb$resmatrix[,"M"], ub$resmatrix[,"M"])
    }
  }

  #pull out quantiles
  N_quants <- quantile(N_sigma, quants_to_assess, na.rm = T)
  N_sigma_rd[j, ] <- c(reg_names[r],  N_quants)
  M_quants <- quantile(M_sigma, quants_to_assess, na.rm = T)
  M_sigma_rd[j,] <- c(reg_names[r], M_quants)
  ind_85 <- which.min(abs(N_sigma - N_quants[2]))

  #plot results
  plot(seq(1, n_iter), res[[r]]$sigma[ind_85,], type= "l",
          xlab = "Iteration", ylab = expression(sigma),
          main = sprintf("%s", reg_names[j]))
}
#dev.off()
xtable(N_sigma_rd)


#--------------------------
#Assess mu
#--------------------------
N_mu_rd <- M_mu_rd <- matrix(ncol = length(quants_to_assess) + 1, nrow = n_reg)
colnames(N_mu_rd) <- colnames(M_mu_rd) <- c("Region",  
                                            paste("n_", quants_to_assess, sep = ""))
#png(filename = sprintf("/users/hdirector/Dropbox/SeaIce_InProgress/probContours_ECMWF/Paper/Figures/traceplot_mu.png",
#                       r), res = 100)
par(mfrow = c(3, 1), oma = rep(1, 4), mar = c(2, 4, 2, 2))
for (j in 1:n_reg) {
  r <- regs_to_fit[j]
  n_lines <- dim(res[[r]]$sigma)[1]
  N_mu <- rep(NA, n_lines)
  M_mu <- rep(NA, n_lines)
  for (i in 1:n_lines) {
    #compute raftery.diag()
    lb <- raftery.diag(res[[r]]$mu[i,], q = .025, r = .0125)
    ub <- raftery.diag(res[[r]]$mu[i,], q = .975, r = .0125)
    N_mu[i] <- max(lb$resmatrix[,"N"], ub$resmatrix[,"N"])
    M_mu[i] <- max(lb$resmatrix[,"M"], ub$resmatrix[,"M"])
  }
  
  #pull out quantiles
  N_quants <- quantile(N_mu, quants_to_assess, na.rm = T)
  N_mu_rd[j, ] <- c(reg_names[j],  N_quants)
  M_mu <- quantile(M_sigma, quants_to_assess, na.rm = T)
  M_mu_rd[j,] <- c(reg_names[j], M_quants)
  ind_85 <- which.min(abs(N_mu - N_quants[2]))
  
  #plot results
  plot(seq(1, n_iter), res[[r]]$mu[ind_85,], type= "l",
       xlab = "Iteration", ylab = expression(sigma),
       main = sprintf("%s", reg_names[j]))
}
#dev.off()
xtable(N_mu_rd)

#--------------------------
#Assess kappa
#--------------------------
n_kappa_rd <- matrix(ncol = 2, nrow = n_reg)
colnames(n_kappa_rd) <- c("Region",  "n")
#png(sprintf("/users/hdirector/Dropbox/SeaIce_InProgress/probContours_ECMWF/Paper/Figures/traceplot_rho.png"),
#    res = 100)
par(mfrow = c(3, 1), oma = rep(1, 4), mar = c(2, 4, 2, 2))
for (j in 1:n_reg) {
  r <- regs_to_fit[j]
  lb <- raftery.diag(res[[r]]$kappa, q = .025, r = .0125)
  ub <- raftery.diag(res[[r]]$kappa, q = .975, r = .0125)
  n_kappa_rd[j,] <- max(lb$resmatrix[,"N"], ub$resmatrix[,"N"])
  plot(1:n_iter, res[[r]]$kappa, type= 'l', xlab = "Iteration",
       ylab = expression(kappa), main = sprintf("Traceplot, %s", reg_names[j]))
}
#dev.off()
xtable(n_kappa_rd)


