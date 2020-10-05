#-------------------------------------------------------------------------------
# Script to make Figure 4 in main paper showing four different forecasts
#
# Code to compute the calibration numbers is currently commented out. Figures 
# can be generated with just the summary data table: calib.rda.
#
# To run sections of code that are currently commented out, the following data 
# and results are required
# Requires Data: ecmwfsipn_sip (arrays named 'sip' with proportion of members
#                               predicting sea ice from the ECWMF ensemble
#                               on Polar Stereographic grid for each 
#                               initialization month. Dimensions are year
#                               x month x longitude x latitude. Forecasts from 
#                               1993 - 2018)
#                bootstrapV3_1 (Bootstrap sea ice observations, version 3.1
#                               downloaded from the National Snow and Ice Data
#                               Center, in original binary form)
#
# Requires Results: clim_prob (computed probabilistic climatology reference 
#                             forecasts obtained with clim_ref.R script)
#                   cont_prob (computed probabilistic forecast from contour 
#                              model, obtained with fitAndGen.R script)
#                   mcf_prob (computed probabilistic MCF forecast, obtained with
#                             wght_EM.R script)
#-------------------------------------------------------------------------------
#library
library("IceCast")
library("fields")
library("viridis") 

#general info
nX <- 304; nY <- 448
xmn = -3850; xmx = 3750; ymn = -5350; ymx = 5850
xBd <- seq(xmn, xmx, 25)
yBd <- seq(ymn, ymx, 25)
n_train <- 10
sip_filepath <- "Data/ecmwfsipn/forecast/ecmwfsipn_sip"
sip_start_year <- 1993

#demo setting of interest
month <- 9
lag <- 1
year <- 2008

#Load observations
obs_file_path <-  "Data/bootstrapV3_1/"
obs <- read_monthly_BS(start_year = year, 
                       end_year = year, 
                       version = 3.1,
                       file_folder = obs_file_path)

###load one observation and one prediction and identify differences in NA patterns
#between dynamic model and post-processed
load(sprintf("%s/initMonth1.rda", sip_filepath))
dyn_prob <- sip[1,1,,]
obs_temp <- obs[1, 1,,]
obs_temp[obs_temp == 120] <- NA #land index
obs_temp[obs_temp == 110] <- 100 #satellite hole is ice
obs_curr <- matrix(nrow = nX, ncol = nY)
obs_curr[obs_temp >= .15] <- 1
obs_curr[obs_temp < .15] <- 0
obs_curr[land_mat == 1] <- NA
NA_in_dyn <- which(is.na(dyn_prob) & !is.na(obs_curr))
rm(sip)

#info about forecasts
train_start_year <- year - n_train
train_end_year <- year - 1
init_month <- get_init_month(month, lag)

#region and ocean masks
all_regions_mask <- conv_to_grid(all_regions)
non_reg_ocean <- matrix(nrow = nX, ncol = nY, data = 0)
non_reg_ocean[all_regions_mask == 0] <- 1

#read in climatology probabilistic (clim_prob/climatology)
load(sprintf("Results/clim_prob/clim_prob_month%i_year%i.rda", month, year))
clim_prob[NA_in_dyn] <- NA
clim_prob[is.na(clim_prob)] <- 101/100
clim_prob[non_reg_ocean == 1] <- 102/100

#read in dynamic ensemble forecast (dyn_prob/ensemble)
load(sprintf("%s/initMonth%i.rda", sip_filepath, init_month))
dyn_prob <- sip[length(sip_start_year:year), month,,]
dyn_prob[is.na(dyn_prob)] <- 101/100
dyn_prob[non_reg_ocean == 1] <- 102/100

#read in contour probabilistic
f <- Sys.glob(file.path('Results/cont_prob', 
                        sprintf("cont_prob_Task*_Month%i_Year%i_Train%i_%i_Init%i.rda",
                                month, year, train_start_year, train_end_year, init_month)))
load(f)

cont_prob[NA_in_dyn] <- NA
cont_prob[is.na(cont_prob)] <- 101/100
cont_prob[non_reg_ocean == 1] <- 102/100

#read MCF prob
load(sprintf("Results/mcf_prob/mcf_prob_month%i_year%i_train%i_%i_init%i.rda",
             month, year, train_start_year, train_end_year, init_month))
mcf_prob[is.na(mcf_prob)] <- 101/100
mcf_prob[non_reg_ocean == 1] <- 102/100

#make contour of observation for comparison
obs_temp <- obs[1, month,,]
obs_temp[obs_temp == 120] <- NA #land index
obs_temp[obs_temp == 110] <- 100 #satellite hole is ice
obs_curr <- matrix(nrow = nX, ncol = nY)
obs_curr[obs_temp >= .15] <- 1
obs_curr[obs_temp < .15] <- 0
obs_curr[land_mat == 1] <- NA
obs_curr[NA_in_dyn] <- NA
obs_poly <- get_region(obs_curr, dat_type = "simple")

#Visualize all forecasts
con_color <- "red"
con_lwd <- 2
cwid = 8
main_size = 1.6
#pdf("Paper/Figures/visSep.pdf", height = 7, width = 7)
par(oma = c(.1, 1, 3, 3), mar = c(.1, 1, 1.5, 0))
layout(matrix(nrow = 14, ncol = 2*cwid + 1, byrow = TRUE,
              data = c(rep(c(rep(1:2, each = cwid), 6), 6),
                       rep(c(rep(3:4, each = cwid), 6), 6),
                       c(rep(5, 2*cwid), 6),
                       c(rep(5, 2*cwid), 6))))
xBdInd <- 60:215; xBdN <- length(xBdInd) 
yBdInd <- 130:310; yBdN <- length(yBdInd)
image(xBd[xBdInd], yBd[yBdInd], dyn_prob[xBdInd[1:(xBdN - 1)], yBdInd[1:(yBdN - 1)]], 
      xaxt = "n",   yaxt = "n", col = c(rep(viridis(10), each = 10), "grey55", "black"),
      main = "Ensemble", zlim = c(0, 102/100), xlab = "", ylab = "",
      cex.main = main_size)
plot(obs_poly, add = T, border = con_color, lwd = con_lwd)
image(xBd[xBdInd], yBd[yBdInd], cont_prob[xBdInd[1:(xBdN - 1)], yBdInd[1:(yBdN - 1)]],
      xaxt = "n",   yaxt = "n", col = c(rep(viridis(10), each = 10), "grey55", "black"),
      zlim = c(0, 102/100), xlab = "", ylab = "", main = "Contour",
      cex.main = main_size)
plot(obs_poly, add = T, border = con_color, lwd = con_lwd)
image(xBd[xBdInd], yBd[yBdInd], clim_prob[xBdInd[1:(xBdN - 1)], yBdInd[1:(yBdN - 1)]], 
      xaxt = "n",  yaxt = "n", col = c(rep(viridis(10), each = 10), "grey55", "black"),
      zlim = c(0, 102/100), xlab = "", ylab = "",
      main = "Climatology", cex.main = main_size)
plot(obs_poly, add = T, border = con_color, lwd = con_lwd)
image(xBd[xBdInd], yBd[yBdInd], mcf_prob[xBdInd[1:(xBdN - 1)], yBdInd[1:(yBdN - 1)]],
      xaxt = "n", yaxt = "n", col = c(rep(viridis(10), each = 10), "grey55", "black"),
      zlim = c(0, 102/100), xlab = "", ylab = "", main = "Mixture Contour Forecast",
      cex.main = main_size)
plot(obs_poly, add = T, border = con_color, lwd = con_lwd)
mtext(sprintf("%s %i, Lead Time %s Months", month.abb[month], year,
              lag +  0.5), outer= TRUE, cex = 1.2, font = 2, line = .5)
plot.new()
legend("center", fill = c("grey55", "black", NA), 
       legend = c("land", "Outside Region", "Observed Contour"),
       lty = c(NA, NA, 1), pch = c(22, 22, NA), ncol = 2, cex = 2,
       col = c(NA, NA, "red"), lwd = c(NA, NA, 8),
       bg = 'white', border = c("black", "black", "white"),
       pt.cex = c(3, 3, 3), text.font = 2, bty = "n")
plot.new()
image.plot(mcf_prob, zlim = c(0, 1), col = viridis(10), legend.only = TRUE,
           legend.width = 8, legend.shrink = .8, axis.args=list(cex.axis=1.8))
#dev.off()


