#-------------------------------------------------------------------------------
# Script to evaluate length of training period for weighting in MCF (produces
# Table 4 in Supplement)
#
# Requires Data: ecmwfsipn_sip (arrays named 'sip' with proportion of members
#                               predicting sea ice from the ECWMF ensemble
#                               on Polar Stereographic grid for each 
#                               initialization month. Dimensions are year
#                               x month x longitude x latitude. Forecasts from 
#                               1993 - 2018)
#                bootstrapV3_1 (Bootstrap sea ice observations, version 3.1
#                               downloaded from the National Snow and Ice Data
#                               Center, in original binary form)
#                psn25area_v3.mat (Matlab file with dimension 304 x 448 that 
#                                  gives area of which grid box in Polar 
#                                  Stereographic grid)
#
# Requires Results: clim_prob (computed climatology reference forecasts obtained
#                              with clim_ref.R script
#                   cont_prob (forecast from contour model, obtained with 
#                              fitAndGen.R script)
#-------------------------------------------------------------------------------

library("IceCast")
library("tidyverse")

#' Function to compute Brier Score
#' @param predField matrix of predicted values (either binary or probabilities)
#' @param obsField matrix of observed values (binary)
#' @param prop_area matrix of area weight in each grid box
field_brier <- function(pred_field, obs_field, prop_area) {
  pred_field[non_reg_ocean == 1] <- NA
  obs_field[non_reg_ocean == 1] <- NA
  #data checks
  stopifnot(dim(pred_field) == dim(obs_field))
  stopifnot(is.na(pred_field) == is.na(obs_field))
  stopifnot(sum(is.na(pred_field)) == sum(is.na(obs_field)))
  
  #compute area-weighted brier score
  brier <- sum(prop_area*((pred_field - obs_field)^2), na.rm = T)
  return(brier)
}


#Set up
years <- 2005:2016
train_lengths <- 1:7
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

#load one observation and one prediction to identify differences in NA patterns
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
ligrid_area <- temp$area
grid_area[is.na(dyn_prob)] <- NA
grid_area[non_reg_ocean == 1] <- NA
tot_area <- sum(grid_area, na.rm = T)
prop_area <- grid_area/tot_area
stopifnot(sum(prop_area, na.rm = T) == 1)

# compute weights for all the different training periods
weights_all <- list()
n_test <- n_years - max(train_lengths) 
for (n_wght_years in train_lengths) {
  weights <- array(dim = c(n_test, n_months, n_lags))
  for (y in (n_years - n_test + 1):n_years) {
    wght_years <- (years[y] - n_wght_years):(years[y] - 1)
    for (m in months) {
      for (l in 1:n_lags) {
        #pull out obs in current month during training period
        obs <- obs_all[wght_years - years[1] + 1, months[m],,]
        if (n_wght_years == 1) {
          temp <- array(dim = c(1, nX, nY))
          temp[1,,] <- obs
          obs <- temp
        }
        for (i in 1:dim(obs)[1]) {
            obs[i,,][NA_in_dyn] <- NA
            obs[i,,][land_mat == 1] <- NA
            obs[i,,][non_reg_ocean == 1] <- NA
        }

        #download both forecasts and organize
        clim <- cont <- array(dim = dim(obs))

        #climatology probabilistic ("clim_prob")
        for (k in 1:n_wght_years) {
          train_start_year <- wght_years[k] - stat_train
          train_end_year <- wght_years[k] - 1

          #climatology probabilistic ("clim_prob")
          load(sprintf("Results/clim_prob/clim_prob_month%i_year%i.rda", months[m], wght_years[k]))
          clim_prob[NA_in_dyn] <- NA
          clim_prob[non_reg_ocean == 1] <- NA
          clim[k,,] <- clim_prob

          #contour probabilistic ("cont_prob")
          init_month <- get_init_month(months[m], lags[l])
          f <- Sys.glob(file.path('Results/cont_prob', 
                                  sprintf("cont_prob_Task*_Month%i_Year%i_Train%i_%i_Init%i.rda",
                                          months[m], wght_years[k], 
                                          train_start_year, train_end_year, 
                                          init_month)))
          load(f)
          cont_prob[NA_in_dyn] <- NA
          cont_prob[non_reg_ocean == 1] <- NA
          cont[k,,] <- cont_prob
        }

        #Fit and store the weights
        weights[y - (n_years - n_test), m, l] <- fit_weights(mod1 = cont, mod2 = clim,
                                                       obs = obs, prop_area = prop_area)
        print(c(y, m, l))
      }
    }
  }
  weights_all[[n_wght_years]] <- weights
}

#compute weighted forecasts conditional on the weights
test_years <- (n_years - n_test + 1):n_years
wghts_exper <- expand.grid(years[test_years], months, lags, train_lengths)
wghts_exper <- cbind(wghts_exper, rep(NA, nrow(wghts_exper)))
colnames(wghts_exper) <- c("year", "month", "lag", "train_lengths", "brier")
for (y in test_years) {
  for (m in 1:n_months) {
    for (l in 1:n_lags) {

      #observation for comparison
      obs_temp <- obs_all[years[y] - years[1] + 1, months[m],,]
      obs_temp[obs_temp == 120] <- NA #land index
      obs_temp[obs_temp == 110] <- 100 #satellite hole is ice
      obs_curr <- matrix(nrow = nX, ncol = nY)
      obs_curr[obs_temp >= .15] <- 1
      obs_curr[obs_temp < .15] <- 0
      obs_curr[land_mat == 1] <- NA
      obs_curr[NA_in_dyn] <- NA

      #general info
      init_month <- get_init_month(months[m], lags[l])
      train_start_year <- years[y] - stat_train
      train_end_year <- years[y] - 1

      #climatology probabilistic ("clim_prob")
      load(sprintf("Results/clim_prob/clim_prob_month%i_year%i.rda", months[m], years[y]))
      clim_prob[NA_in_dyn] <- NA

      #contour probabilistic ("cont_prob")
      f <- Sys.glob(file.path('Results/cont_prob', 
                              sprintf("cont_prob_Task*_Month%i_Year%i_Train%i_%i_Init%i.rda",
                                      months[m], years[y], train_start_year,
                                      train_end_year, init_month)))
      load(f)
      cont_prob[NA_in_dyn] <- NA

      for (t in train_lengths) {
        wght_prob <- wght_mod(weights_all[[t]][y - test_years[1] + 1, m, l],
                              cont_prob, clim_prob)
        ind <- which(wghts_exper$year == years[y] & wghts_exper$month == months[m]
                     & wghts_exper$lag == lags[l] & wghts_exper$train_lengths == t)
        stopifnot(length(ind) == 1)
        wghts_exper[ind,]$brier <- field_brier(pred_field = wght_prob, 
                                               obs_field = obs_curr,
                                               prop_area = prop_area)
        rm(wght_prob)

      }
      rm(clim_prob); rm(cont_prob);
      print(c(y, m, l))
    }
  }
}

#results table
wghts_exper_sum <- wghts_exper %>%
  group_by(train_lengths) %>%
  dplyr::summarize(mean = mean(brier))
wghts_exper_sum

library("xtable")
xtable(wghts_exper_sum[,1:2], digits = 5)

