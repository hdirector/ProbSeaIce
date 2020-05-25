rm(list = ls())
library("IceCast")

#Set up
n_wght_years <- 3
years <- 2005:2016
n_years <- length(years)
nX <- 304; nY <- 448
months <- 1:12; lags <- 0:6
n_months <- length(months); n_lags <- length(lags)
stat_train <- 10
sip_filepath <- "Data/ecmwfsipn/forecast/ecmwfsipn_sip"

#region and ocean masks
all_regions_mask <- conv_to_grid(all_regions)
non_reg_ocean <- matrix(nrow = nX, ncol = nY, data = 0)
non_reg_ocean[all_regions_mask == 0] <- 1

#Load and process observations
obs_file_path <-  "Data/bootstrapV3_1/"
obs_all_temp <- read_monthly_BS(start_year = years[1], end_year = years[length(years)],
                                version = 3.1, file_folder = obs_file_path)
obs_all <- obs_all_temp
obs_all[obs_all_temp == 120] <- NA #land index
obs_all[obs_all == 110] <- 100 #satellite hole is ice
obs_all <- obs_all/100
obs_all[obs_all >= .15] <- 1
obs_all[obs_all < .15] <- 0

###load one observation and one prediction to identify differences in NA patterns
#between dynamic model and post-processed
load(sprintf("%s/initMonth1.rda", sip_filepath))
dyn_prob <- sip[1,1,,]
obs_temp <- obs_all[1, 1,,]
obs_temp[obs_temp == 120] <- NA #land index
obs_temp[obs_temp == 110] <- 100 #satellite hole is ice
obs_curr <- matrix(nrow = nX, ncol = nY)
obs_curr[obs_temp >= .15] <- 1
obs_curr[obs_temp < .15] <- 0
obs_curr[land_mat == 1] <- NA
NA_in_dyn <- which(is.na(dyn_prob) & !is.na(obs_curr))
rm(sip)

#load areas of each grid box and compute area weighting
library("R.matlab")
temp <- readMat("Data/grids/psn25area_v3.mat")
grid_area <- temp$area
grid_area[is.na(dyn_prob)] <- NA
grid_area[non_reg_ocean == 1] <- NA
tot_area <- sum(grid_area, na.rm = T)
prop_area <- grid_area/tot_area
stopifnot(sum(prop_area, na.rm = T) == 1)

#compute  weights
# weights <- array(dim = c(n_years - n_wght_years, n_months, n_lags))
# for (y in (n_wght_years + 1):n_years) {
#   wght_years <- (years[y] - n_wght_years):(years[y] - 1)
#   for (m in months) {
#     for (l in 1:n_lags) {
#       #pull out obs in current month during training period
#       obs <- obs_all[wght_years - years[1] + 1, months[m],,]
#       for (i in 1:dim(obs)[1]) {
#         obs[i,,][NA_in_dyn] <- NA
#         obs[i,,][land_mat == 1] <- NA
#         obs[i,,][non_reg_ocean == 1] <- NA
#       }
# 
#       #download both forecasts and organize
#       clim <- conts <- array(dim = dim(obs))
# 
#       #climatology probabilistic ("clim_prob")
#       for (k in 1:n_wght_years) {
#         train_start_year <- wght_years[k] - stat_train
#         train_end_year <- wght_years[k] - 1
# 
#         #climatology probabilistic ("clim_prob")
#         load(sprintf("Results/clim_prob/clim_prob_month%i_year%i.rda", months[m], wght_years[k]))
#         clim_prob[NA_in_dyn] <- NA
#         clim_prob[non_reg_ocean == 1] <- NA
#         clim[k,,] <- clim_prob
# 
#         #contour model ("cont_prob")
#         init_month <- get_init_month(months[m], lags[l])
#         f <- Sys.glob(file.path('Results/cont_prob', 
#                       sprintf("cont_prob_Task*_Month%i_Year%i_Train%i_%i_Init%i.rda",
#                       months[m], wght_years[k], train_start_year, train_end_year, 
#                       init_month)))
#         load(f)
#         cont_prob[NA_in_dyn] <- NA
#         cont_prob[non_reg_ocean == 1] <- NA
#         conts[k,,] <- cont_prob
#       }
# 
#       #Fit and store the weights
#       weights[y - n_wght_years, m, l] <- fit_weights(mod1 = conts, mod2 = clim,
#                                                      obs = obs, prop_area = prop_area)
#       print(c(y, m, l))
#     }
#   }
# }
#save(weights, file = "Results/summaries/weights.rda")
load("Results/summaries/weights.rda")


#plot weights for understanding
library("fields")
library("viridis")
pdf("Paper/Figures/av_weights.pdf", height = 4, width = 8)
weights_av <- apply(weights, 2:3, mean)
image.plot(1:12, 1:7, weights_av, col = viridis(10), zlim = c(0.18, 0.82),
           main = "Weight on Post-Proccessed Ensemble", xaxt = "n",
           yaxt = "n", ylab = "Lead Time", xlab = "Month")
axis(1, labels = c("", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",
                   "Aug", "Sep", "Oct", "Nov", "Dec", ""),
     at = 0:13)
axis(2, labels = c("", 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, ""),
     at = 0:8, las = 1)
points(which(weights_av >= 0.4, arr.ind = T), pch = 20)
dev.off()

#compute weighted forecasts conditional on the weights
for (y in (n_wght_years + 1):n_years) {
  for (m in 1:n_months) {
    for (l in 1:n_lags) {
      init_month <- get_init_month(months[m], lags[l])
      train_start_year <- years[y] - stat_train
      train_end_year <- years[y] - 1

      #climatology probabilistic ("clim_prob")
      load(sprintf("Results/clim_prob/clim_prob_month%i_year%i.rda", months[m], 
                   years[y]))
      clim_prob[NA_in_dyn] <- NA
      clim_prob[non_reg_ocean == 1] <- NA

      #cont model probabilistic ("cont_prob")
      f <- Sys.glob(file.path('Results/cont_prob', 
                              sprintf("cont_prob_Task*_Month%i_Year%i_Train%i_%i_Init%i.rda",
                                      months[m], years[y], train_start_year,
                                      train_end_year, init_month)))
      load(f)
      cont_prob[NA_in_dyn] <- NA
      cont_prob[non_reg_ocean == 1] <- NA
      mcf_prob <- wght_mod(weights[y - n_wght_years, m, l], cont_prob,
                            clim_prob)
      save(mcf_prob, file = sprintf("Results/mcf_prob/mcf_prob_month%i_year%i_train%i_%i_init%i.rda",
                                     months[m], years[y], train_start_year, train_end_year, init_month))

      mcf_bin <- matrix(nrow = nX, ncol = nY, data = NA)
      mcf_bin[mcf_prob >= .5] <- 1
      mcf_bin[mcf_prob <.5] <- 0
      save(mcf_bin, file = sprintf("Results/mcf_bin/mcf_bin_month%i_year%i_train%i_%i_init%i.rda",
                                    months[m], years[y], train_start_year, train_end_year, init_month))

      rm(clim_prob); rm(cont_prob); rm(mcf_prob); rm(mcf_bin)
      print(c(y, m, l))
    }
  }
}


