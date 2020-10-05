#-------------------------------------------------------------------------------
# Script to compute contour forecast for a particular year, month, lead 
# times, and training length (one line in task_table). Designed to be used on
# a cluster
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
#
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#General set up
#-------------------------------------------------------------------------------
set.seed(104)
n_gen <- 100

#load task information
args <- commandArgs(trailingOnly =TRUE)
task_id <- as.numeric(args[1])
task_path <- ".../exper_design/ecmwfExper.rda" 
load(task_path)
task <- task_table[task_id,]
attach(task)
print(sprintf("Executing fitContours taskID  %i:, month: %i, year %i, training years: %i - %i, lag %i",
              task_id, month, forecast_year, train_start_year, train_end_year, lag))
eps <- .01

#info on observations
obs_file_path <- '.../Data/bootstrapV3_1/' 
bs_version <- 3.1
dat_type_obs <- "bootstrap"

#info on dynamic model
dyn_path <- ".../Data/ecmwfsipn/ecmwfsipn_sip" 
dyn_mod <- "ecmwfsipn"
dyn_start_year <- 1993
dat_type_pred <- "simple"
train_bc_start_year <- 1993

#misc fixed constants
level <- 15
n_iter <- 55000
burn_in <- 5000

#load package
library("IceCast")

#-------------------------------------------------------------------------------
#Load observations
#-------------------------------------------------------------------------------
#read in obs
obs <- read_monthly_BS(start_year = train_bc_start_year,
                       end_year = train_end_year, version = bs_version,
                       file_folder = obs_file_path)
print("loaded observations")

#-------------------------------------------------------------------------------
#Load dynamic model predictions and convert to binary field
#-------------------------------------------------------------------------------
init_month <- get_init_month(month, lag)
#determine if dyn_model has predictions for first year in specified month
if (init_month > month) {
  train_bc_start_year <-  train_bc_start_year + 1
}
load(sprintf("%s/initMonth%i.rda", dyn_path, init_month))
dyn_bin <- sip
dyn_bin[sip >= 0.5] <- 1
dyn_bin[sip < 0.5] <- 0
print("loaded dynamic model forecasts")

#compute lengths ice extends y
y_bc <- find_y(start_year = train_bc_start_year, end_year = train_end_year,
               obs_start_year = train_bc_start_year,
               pred_start_year = dyn_start_year, observed = obs[,month,,],
               predicted = dyn_bin[,month,,], reg_info, month, level,
               dat_type_obs, dat_type_pred)

#------------------------------------------------------------------------------
#Contour-shift for bias correction
#------------------------------------------------------------------------------
#apply contour-shifting
cont_bin_poly <- contour_shift(y = y_bc,
                             predicted = dyn_bin[length(dyn_start_year:forecast_year),
                                                 month,,],
                             bc_year = forecast_year,
                             pred_start_year = dyn_start_year, reg_info, level,
                             dat_type_pred)

cont_bin <- conv_to_grid(cont_bin_poly)
save(cont_bin, file = sprintf(".../results/cont_bin/cont_bin_Task%i_Month%i_Year%i_Train%i_%i_Init%i.rda",
                              task_id,  month, forecast_year, train_start_year,  train_end_year, init_month))


#------------------------------------------------------------------------------
#model fitting
#------------------------------------------------------------------------------
#find y for training period for observations
y_years_bc <- y_bc$start_year:y_bc$end_year
n_bc_years <- length(y_years_bc)
train_ind <- (1:n_bc_years)[y_years_bc %in% train_start_year:train_end_year]

y_train <- find_y(start_year = train_start_year, end_year = train_end_year,
                  obs_start_year = train_start_year,
                  pred_start_year = NULL, observed = obs[train_ind, month,,],
                  predicted = NULL, reg_info, month, level,
                  dat_type_obs, dat_type_pred, obs_only = TRUE)

#Determine which regions to fit, and the 'full' polygon
temp <- to_fit(y_obs = y_train$obs, reg_info)
regs_to_fit <- temp$regs_to_fit
full <- temp$full

#convert lengths to proportions 
prop_train <- y_to_prop(y = y_train$obs, regs_to_fit, reg_info)
prop_train <- lapply(prop_train, function(y){sapply(y, function(x){x})})

#convert observations to proportions and logit of proportions
prop_train_tilde <- list()
for (r in regs_to_fit) {
  ub_ind <- which(prop_train[[r]] >= 1 - eps)
  prop_train[[r]][ub_ind] <- 1 - eps
  lb_ind <- which(prop_train[[r]] <= eps)
  prop_train[[r]][lb_ind] <- eps
  prop_train_tilde[[r]] <- logit(prop_train[[r]])
}

#convert bias-corrected lengths to proportions
y_bc <- find_y_1(ice = cont_bin_poly, reg_info)
prop_bc <- y_to_prop(y = y_bc, regs_to_fit, reg_info)
for (r in regs_to_fit) {
  prop_bc[[r]][prop_bc[[r]] >= 1 - eps] <- 1 - eps
  prop_bc[[r]][prop_bc[[r]] <= eps] <- eps
}

#sigma bounds
ub_props <- c(.99, .99, .99, .99, .72)
lb_props <- c(.15, .01, .01, .01, .01)

#Run MCMC chains
print("starting MCMC chains")
res <- list()
for (r in regs_to_fit)  {
  start_time <- proc.time()
  res[[r]] <- fit_cont_pars(r = r, n_iter = n_iter,
                            prop_tilde = prop_train_tilde[[r]],
                            reg_info = reg_info, prop0 = prop_bc[[r]],
                            ub_prop = ub_props[r], lb_prop = lb_props[r])
  end_time <- proc.time()
  elapse_time <- end_time - start_time
  print(sprintf("MCMC for region %i finished, elapsed time %f", r, elapse_time[3]))
}

#save all chain info in a few cases
if (forecast_year == 2005 & lag  == 1) {
	save(res, file = sprintf(".../results/cont_fits/cont_fit_Task%i_Month%i_Year%i_Train%i_%i_Init%i.rda",
							 task_id, month, forecast_year, train_start_year, train_end_year, init_month))
}
#Compute mu and sigma for each region
pars <- list()
for (r in regs_to_fit) {
  pars[[r]] <- calc_pars(res_r = res[[r]], burn_in, r = r)
}

print("fitting complete")

#-------------------------------------
#generate contours
#-------------------------------------
#generate contours
indiv_conts <- list()
for (r in regs_to_fit) {
  indiv_conts[[r]] <- gen_cont(r, pars_r = pars[[r]], reg_info, n_gen,
                               rand = res[[r]]$rand, start = res[[r]]$start,
                               start_inds = res[[r]]$start_inds,
                               end = res[[r]]$end, end_inds = res[[r]]$end_inds)
  print(sprintf("contours for region %i generated", r))
}

conts <- merge_conts(conts = indiv_conts, full = full)
save(conts, file = sprintf(".../results/cont_conts/cont_conts_Task%i_Month%i_Year%i_Train%i_%i_Init%i.rda",
                           task_id, month, forecast_year, train_start_year, train_end_year, init_month))
cont_prob <- prob_map(merged = conts)
save(cont_prob, file = sprintf(".../results/cont_prob/cont_prob_Task%i_Month%i_Year%i_Train%i_%i_Init%i.rda",
                                task_id, month, forecast_year, train_start_year, train_end_year, init_month))
print("completed cont_prob")

