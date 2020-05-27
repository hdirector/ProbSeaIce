#Script to produce all Brier score plots in paper and supplement

#Note that abbreviated names of forecasts are used here compared to the paper
# dPersis = Damped persistence, clim_bin = binary climatology,
# clim_prob = probabilistic climatology, dyn_bin = binary ensemble
# dyn_prob = probabilistic ensemble, 
# cont_prob = post-processed ensemble probabilistic
# cont_bin = post-processed ensemble binary (contour-shifting), 
# taqm = trend adjusted quantile mapping (Dirkson et al 2019),
# mcf_prob = Mixture Contour Forecasts probabilistic
# mcf_bin = Mixture Contour Forecasts binary


#set up
library("IceCast")
library('RcppCNPy')
months <- 1:12
years <- 2008:2016
lags <- 0:6
nX <- 304; nY <- 448
stat_train <- 10
sip_filepath <- "/Users/hdirector/Dropbox/SeaIce_InProgress/probContours_ECMWF/Data/ecmwfsipn/forecast/ecmwfsipn_sip"
sip_start_year <- 1993

#identify non-regional ocean
all_regions_mask <- conv_to_grid(all_regions)
non_reg_ocean <- matrix(nrow = nX, ncol = nY, data = 0)
non_reg_ocean[all_regions_mask == 0] <- 1

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


#basic info
n_months <- length(months)
n_years <- length(years)
n_lags <- length(lags)

#Load observations,
obs_file_path <-  "Data/bootstrapV3_1/"
obs <- read_monthly_BS(start_year = years[1], 
                       end_year = years[length(years)], 
                       version = 3.1,
                       file_folder = obs_file_path)

#make data frame of all test cases and models
cases <- expand.grid(years, months, lags)
colnames(cases) <- c("year", "month", "lag")
fcasts <- matrix(nrow = nrow(cases), ncol = 10)
colnames(fcasts) <- c("dPersis", "clim_bin", "clim_prob",  "dyn_bin", "dyn_prob",
                      "cont_bin", "cont_prob", "taqm", "mcf_prob", "mcf_bin")
brier <- data.frame(cbind(cases, fcasts))


###load one observation and one prediction to identify differences in NA patterns
#between dynamic model and post-processed
load(sprintf("%s/initMonth1.rda", sip_filepath))
dyn_prob <- sip[1,1,,]
#observation for comparison
obs_temp <- obs[1, 1,,]
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

# #Compute brier scores for all cases and forecast types
# for (y in 1:n_years) {
#   for (m in 1:n_months) {
#     #observation for comparison
#     obs_temp <- obs[length(years[1]:years[y]), months[m],,]
#     obs_temp[obs_temp == 120] <- NA #land index
#     obs_temp[obs_temp == 110] <- 100 #satellite hole is assumed to be ice
#     obs_curr <- matrix(nrow = nX, ncol = nY)
#     obs_curr[obs_temp >= .15] <- 1
#     obs_curr[obs_temp < .15] <- 0
#     obs_curr[land_mat == 1] <- NA
#     obs_curr[NA_in_dyn] <- NA
# 
#     #identify years where training begins for stat model
#     train_start_year <- years[y] - stat_train
#     train_end_year <- years[y] - 1
# 
#     ###compute brier scores for forecasts where the lag does not matter
#     inds <- which(brier$year == years[y] & brier$month == months[m])
#     stopifnot(length(inds) == n_lags)
# 
#     #climatology binary ("clim_bin")
#     load(sprintf("Results/clim_bin/clim_bin_month%i_year%i.rda", months[m], years[y]))
#     clim_bin[NA_in_dyn] <- NA
#     brier[inds,]$clim_bin <- field_brier(pred_field = clim_bin, obs_field = obs_curr,
#                                          prop_area)
# 
#     #climatology probabilistic ("clim_prob")
#     load(sprintf("Results/clim_prob/clim_prob_month%i_year%i.rda", months[m], years[y]))
#     clim_prob[NA_in_dyn] <- NA
#     brier[inds,]$clim_prob <- field_brier(pred_field = clim_prob, obs_field = obs_curr,
#                                           prop_area)
# 
#     #compute brier scores for forecasts where the lag matters
#     for (j in 1:n_lags) {
#       init_month <- get_init_month(months[m], lags[j])
#       ind <- which(brier$year == years[y] & brier$month == months[m]
#                    & brier$lag == lags[j])
#       stopifnot(length(ind) == 1)
# 
#       #damped persistence ("dPersis")
#       load(sprintf("Results/dPersis/dPersis_month%i_year%i_lag%i.rda", months[m], years[y],
#                    lags[j]))
#       dPersis[NA_in_dyn] <- NA
#       brier[ind,]$dPersis <- field_brier(pred_field = dPersis, obs_field = obs_curr,
#                                           prop_area)
#       rm(dPersis)
# 
#       #ensemble probabilistic forecast ("dyn_prob")
#       load(sprintf("%s/initMonth%i.rda", sip_filepath, init_month))
#       dyn_prob <- sip[length(sip_start_year:years[y]), months[m],,]
#       brier[ind,]$dyn_prob <- field_brier(pred_field = dyn_prob, obs_field = obs_curr,
#                                           prop_area)
# 
#       #ensemble binary forecast ("dyn_bin")
#       dyn_bin <- matrix(nrow = nX, ncol = nY, data = 0)
#       dyn_bin[land_mat == 1] <- NA
#       dyn_bin[NA_in_dyn] <- NA
#       dyn_bin[dyn_prob >= 0.5] <- 1
#       brier[ind,]$dyn_bin <- field_brier(pred_field = dyn_bin, obs_field = obs_curr,
#                                          prop_area)
# 
#       #trend adjusted quantile mapping ("TAQM")
#       taqm <- npyLoad(sprintf("Results/taqm/taqm_month%i_year%i_init%i.npy", months[m],
#                               years[y], init_month))
#       brier[ind,]$taqm <- field_brier(pred_field = taqm, obs_field = obs_curr,
#                                       prop_area)
#       rm(taqm)
# 
#       #post-processed ensemble binary ("cont_bin")
#       f <- Sys.glob(file.path('Results/cont_bin', 
#                               sprintf("cont_bin_Task*_Month%i_Year%i_Train%i_%i_Init%i.rda",
#                                       months[m], years[y], train_start_year, 
#                                       train_end_year, init_month)))
#       load(f)
#       cont_bin[NA_in_dyn] <- NA
#       brier[ind,]$cont_bin <- field_brier(pred_field = cont_bin,
#                                             obs_field = obs_curr, prop_area)
#       rm(cont_bin)
# 
#       #contour probabilistic ("cont_prob")
#       f <- Sys.glob(file.path('Results/cont_prob', 
#                               sprintf("cont_prob_Task*_Month%i_Year%i_Train%i_%i_Init%i.rda",
#                                       months[m], years[y], train_start_year,
#                                       train_end_year, init_month)))
#       load(f)
#       cont_prob[NA_in_dyn] <- NA
#       brier[ind,]$cont_prob <- field_brier(pred_field = cont_prob,
#                                              obs_field = obs_curr, prop_area)
#       rm(cont_prob); rm(dyn_prob); rm(dyn_bin)
# 
#       #mixture contour forecast probabilistic ("mcf_prob")
#       load(sprintf("Results/mcf_prob/mcf_prob_month%i_year%i_train%i_%i_init%i.rda",
#                    months[m], years[y], train_start_year, train_end_year, init_month))
#       brier[ind,]$mcf_prob <- field_brier(pred_field = mcf_prob,
#                                            obs_field = obs_curr, prop_area)
#       rm(mcf_prob)
# 
#       #mixture contour forecast binary ("mcf_bin")
#       load(sprintf("Results/mcf_bin/mcf_bin_month%i_year%i_train%i_%i_init%i.rda",
#                    months[m], years[y], train_start_year, train_end_year, init_month))
#       brier[ind,]$mcf_bin <- field_brier(pred_field = mcf_bin,
#                                           obs_field = obs_curr, prop_area)
#       rm(mcf_bin)
#     }
#     print(c(m, y))
# 
#   }
# }
# save(brier, file = "Results/summaries/brier.rda")
load("Results/summaries/brier.rda")


####Brier Score by  month, lag, model
library("reshape2")
library("tidyverse")
library("dplyr")

#rearrange results
brier_short <- melt(brier, id.vars = c("year", "month", "lag"))
colnames(brier_short)[4:5] <- c("mod", "BSS")

#make lag time numeric and remove month rounding
brier_short$lag <- as.numeric(brier_short$lag) + 0.5

#assign'season' categories
brier_short <- brier_short %>%  
  mutate(season = ifelse(month %in% c(11:12, 1), "Nov/Dec/Jan",
                  ifelse(month %in% 2:4, "Feb/Mar/Apr", 
                  ifelse(month %in% 5:7, "May/Jun/Jul", "Aug/Sep/Oct")))) %>%
  mutate(mod_type = ifelse(mod == 'dyn_prob'|mod == 'cont_prob'|
                           mod == "clim_prob"|mod == 'taqm'| mod == "mcf_prob",
                           "prob", "bin"))
brier_short$season <- factor(brier_short$season, 
                             levels = c("Nov/Dec/Jan", "Feb/Mar/Apr", "May/Jun/Jul",
                                        "Aug/Sep/Oct"))

#brier summary
brier_sum <- brier_short %>%
  group_by(month, lag, season, mod, mod_type) %>%
  dplyr::summarize(mean = mean(BSS, na.rm = T))
brier_sum$month <- factor(month.abb[brier_sum$month], levels = month.abb)

#set colors to ensure consistency among plots
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
six_colors <- gg_color_hue(6)
six_colors <- six_colors[c(2:4, 6, 5, 1)] #order colors
five_colors <- six_colors[c(1:4, 6)]

#formal names of forecasts
prob_formal <- c("Climatology", "Ensemble", "Contour Model",
                 "Mixture Contour Forecast","Trend Adjusted Quantile Mapping",
                 "Damped Persistence")
bin_formal <- c("Climatology", "Ensemble", "Contour-Shifted Ensemble",
                "Mixture Contour Forecast", "Damped Persistence")

#--------------------------------------------------
#Prob. Brier score for months ASO, results section
#--------------------------------------------------
brier_prob_month <- brier_sum %>%
  filter(mod_type == "prob") %>%
  filter(month == "Aug"| month == "Sep"| month == "Oct") %>%
  group_by(lag, mod, month, mod_type) %>%
  dplyr::summarize(mean = mean(mean))

p_brier_prob_ASO <- ggplot(data = brier_prob_month, 
                            aes(x = lag, y = mean, group = mod, col = mod)) +
  geom_point() + 
  geom_line(aes(linetype=mod_type), show.legend = FALSE)+ 
  xlab("Lead Time (in months)") +
  ylab("Mean Brier Score") +
  ylim(0, .083) + 
  theme(legend.position="bottom", plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12)) +
  scale_color_manual(breaks = c("clim_prob", "dyn_prob", "cont_prob", 
                                "mcf_prob", "taqm"),
                     values = six_colors[1:5],
                     labels = prob_formal[1:5]) +
  scale_linetype_manual(breaks = c("prob"),
                        values = c("solid"),
                        labels = c("Probabilistic")) +
  guides(col = guide_legend(nrow = 2, title = "")) +
  ggtitle("Probabilistic  Forecast Performance,
          2008-2016")

pdf("Paper/Figures/brier_prob_ASO.pdf", width = 10, height = 5)
p_brier_prob_ASO + facet_grid(cols = vars(month))
dev.off()

#------------------------------------------------
#Prob. Brier score by season, supplement
#-----------------------------------------------
brier_prob_seas <- brier_sum %>%
  filter(mod_type == "prob"|mod == "dPersis") %>%
  group_by(lag, mod, season, mod_type) %>%
  dplyr::summarize(mean = mean(mean))
p_brier_prob_seas <- ggplot(data = brier_prob_seas, 
                             aes(x = lag, y = mean, group = mod, col = mod)) +
  geom_point() + 
  geom_line(aes(linetype=mod_type), show.legend = FALSE)+ 
  xlab("Lead Time (in months)") +
  ylab("Mean Brier Score") +
  ylim(0, .095) + 
  theme(legend.position="bottom", plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12)) +
  scale_color_manual(breaks = c("clim_prob", "dyn_prob", "cont_prob", 
                                "mcf_prob", "taqm", "dPersis"),
                       values = six_colors,
                       labels = prob_formal) +
  scale_linetype_manual(breaks = c("prob", "bin"),
                        values = c("solid", "dashed"),
                        labels = c("Probabilistic", "Binary")) +
  guides(col = guide_legend(nrow = 2, title = "")) +
  ggtitle("Probabilistic  Forecast Performance,
          2008-2016")

pdf("Paper/Figures/brier_prob_seas.pdf", width = 10, height = 5)
p_brier_prob_seas + facet_grid(cols = vars(season))
dev.off()

#------------------------------------------------
#Binary Brier score by month ASO, results section
#------------------------------------------------
brier_bin_ASO <- brier_sum %>%
  filter(mod_type == "bin") %>%
  filter(month == "Aug"| month == "Sep"| month == "Oct") %>%
  group_by(lag, mod, month, mod_type) %>%
  dplyr::summarize(mean = mean(mean))
p_brier_bin_ASO <- ggplot(data = brier_bin_ASO, 
                           aes(x = lag, y = mean, group = mod, color = mod)) +
  geom_point() +
  geom_line() +
  xlab("Lead Time (in months)") +
  ylab("Mean Brier Score") +
  ylim(0, .103) +
  theme(legend.position ="bottom", plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12)) +
  scale_color_manual(breaks = c("clim_bin", "dyn_bin", "cont_bin",
                                "mcf_bin", "dPersis"),
                     values = five_colors, labels = bin_formal) +
  guides(col = guide_legend(nrow = 2, title = "")) +
  ggtitle("Binary Forecast Performance, 2008-2016")

pdf("Paper/Figures/brier_bin_ASO.pdf", width = 10, height = 5)
p_brier_bin_ASO + facet_grid(cols = vars(month)) 
dev.off()


#------------------------------------------------
#Binary Brier score by season, supplement
#------------------------------------------------
brier_bin_seas <- brier_sum %>%
  filter(mod_type == "bin") %>%
  group_by(lag, mod, season, mod_type) %>%
  dplyr::summarize(mean = mean(mean))
p_brier_bin_seas <- ggplot(data = brier_bin_seas, 
                            aes(x = lag, y = mean, group = mod, color = mod)) +
  geom_point() +
  geom_line() +
  xlab("Lead Time (in months)") +
  ylab("Mean Brier Score") +
  ylim(0, .095) +
  theme(legend.position="bottom", plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12)) +
  scale_color_manual(breaks = c("clim_bin", "dyn_bin", 
                                "cont_bin", "mcf_bin", "dPersis"),
                     values = five_colors, labels = bin_formal) +
  guides(col = guide_legend(nrow = 2, title = "")) +
  ggtitle("Binary Forecast Performance, 2008-2016")

pdf("Paper/Figures/brier_bin_seas.pdf", width = 10, height = 5)
p_brier_bin_seas + facet_grid(cols = vars(season)) 
dev.off()

#-------------------------------
#Overall performance, supplement
#------------------------------
overall_formal <- c("Climatology","Ensemble", "Post-Processed Ensemble",
                    "Mixture Contour Forecast", "Damped Persistence")
#probabilistic forecasts, overall
brier_all <- brier_sum %>%
  filter(mod_type == "prob"|mod == "dPersis") %>%
  group_by(lag, mod, mod_type) %>%
  dplyr::summarize(mean = mean(mean))
p_brier_all <- ggplot(data = brier_all, 
                           aes(x = lag, y = mean, group = mod, col = mod)) +
  geom_point() + 
  geom_line(aes(linetype=mod_type), show.legend = FALSE)+ 
  xlab("Lead Time (in months)") +
  ylab("Mean Brier Score") +
  ylim(0, .065) + 
  theme(legend.position="bottom", plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12)) +
  scale_color_manual(breaks = c("clim_prob", "dyn_prob", "cont_prob", 
                                "mcf_prob", "taqm", "dPersis"),
                     values = six_colors,
                     labels = prob_formal)  +
  scale_linetype_manual(breaks = c("prob", "bin"),
                        values = c("solid", "dashed"),
                        labels = c("Probabilistic", "Binary")) +
  guides(col = guide_legend(nrow = 2, title = "")) +
  ggtitle("Probabilistic Forecasts")

pdf("Paper/Figures/overall_brier.pdf", width = 8, height = 5)
p_brier_all
dev.off()
