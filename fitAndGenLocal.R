#-------------------------------------------------------------------------------
#General set up
#-------------------------------------------------------------------------------
rm(list = ls())
set.seed(103)
n_gen <- 10

#load task information
task_id <- 244
task_path <-"/Users/hdirector/Dropbox/SeaIce_InProgress/probContours_ECMWF/exper_design/ecmwfExper.rda"
load(task_path)
task <- task_table[task_id,]
attach(task)
print(sprintf("Executing fitContours taskID  %i:, month: %i, training years: %i - %i",
              task_id, month, train_start_year, train_end_year))
eps <- .01

#info on observations
obs_file_path <- '/Users/hdirector/Dropbox/SeaIce_InProgress/probContours_ECMWF/Data/bootstrapV3_1/'
bs_version <- 3.1
dat_type_obs <- "bootstrap"

#info on dynamic model
dyn_path <-"/Users/hdirector/Dropbox/SeaIce_InProgress/probContours_ECMWF/Data/ecmwfsipn/forecast/ecmwfsipn_sip"
dyn_mod <- "ecmwfsipn"
dyn_start_year <- 1993
dat_type_pred <- "simple"
train_bc_start_year <- 1993

#misc fixed constants
level <- 15
n_iter <- 5000#150000
burn_in <- 2500#50000

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
bc_bin_poly <- contour_shift(y = y_bc,
                             predicted = dyn_bin[length(dyn_start_year:forecast_year),
                                                     month,,],
                             bc_year = forecast_year,
                             pred_start_year = dyn_start_year, reg_info, level,
                             dat_type_pred)

bc_bin <- conv_to_grid(bc_bin_poly)
#save(bc_bin, file = sprintf(".../%s/bc_bin/bc_bin_month%i_year%i_train%i_%i_init%i.rda",
#                                dyn_mod, month, forecast_year, train_start_year,
#                                 train_end_year, init_month))


#------------------------------------------------------------------------------
#model fitting
#------------------------------------------------------------------------------
#find y for training period for observations
y_years_bc <- y_bc$start_year:y_bc$end_year
n_bc_years <- length(y_years_bc)
train_ind <- (1:n_bc_years)[y_years_bc %in% train_start_year:train_end_year]
y_train_all <- find_y(start_year = train_start_year, end_year = train_end_year,
                      obs_start_year = train_start_year,
                      pred_start_year = NULL, observed = obs[train_ind,month,,],
                      predicted = NULL, reg_info, month, level,
                      dat_type_obs, dat_type_pred, obs_only = TRUE)

#Determine which regions to fit, and the 'full' polygon
temp <- to_fit(y_obs = y_train_all$obs, reg_info)
regs_to_fit <- temp$regs_to_fit
full <- temp$full

#For reg 1, keep only y that correspond to not touching land
y_train <- y_train_all$obs

#convert lengths to proportions and transformed porportions
prop_train <- y_to_prop(y = y_train, regs_to_fit, reg_info)
prop_train <- lapply(prop_train, function(y){sapply(y, function(x){x})})

#remove fixed lines in region 1
prop_train[[1]] <- prop_train[[1]][reg_info$line_fit1,]

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
y_bc <- find_y_1(ice = bc_bin_poly, reg_info)
prop_bc <- y_to_prop(y = y_bc, regs_to_fit, reg_info)
for (r in regs_to_fit) {
  prop_bc[[r]][prop_bc[[r]] >= 1 - eps] <- 1 - eps
  prop_bc[[r]][prop_bc[[r]] <= eps] <- eps
}
prop_bc[[1]] <- prop_bc[[1]][reg_info$line_fit1]

#start_year_prior <- 1981
#end_year_prior <- 1994
#prior info
#load data
# all <- IceCast::read_monthly_BS(start_year, end_year,
#                                 "/Users/hdirector/Dropbox/SeaIce_InProgress/probContours_ECMWF/Data/bootstrapV3_1/",
#                                 version = 3.1)
# march <- apply(all[,3,,], 1, function(x){get_region(x, dat_type = "bootstrap",
#                                            level = 15)})
# #save(march, file = "/users/hdirector/desktop/march.rda")
# load(file = "/users/hdirector/desktop/march.rda")
# y <- lapply(march, function(x){find_y_1(ice = x, reg_info)})
# props <- lapply(y, function(x){y_to_prop(x, regs_to_fit = 1:5, reg_info)})
# prop_max <- sapply(props, function(x){sapply(x, max)})
# ub_prop <- apply(prop_max, 1, max)
# ub_prop[ub_prop >= 1] <- 1 - eps
# lb_prop <- c(.35, rep(eps, 4))

#sigma bounds
ub_prop <- c(.99, .99, .99, .99, .73)
lb_prop <- c(.15, .01, .01, .01, .01)

#Run MCMC chains
print("starting MCMC chains")
res <- list()
for (r in regs_to_fit) {
    start_time <- proc.time()
    res[[r]] <- fit_cont_pars(r = r, n_iter = n_iter,
                              y_obs_r = prop_train_tilde[[r]],
                              reg_info = reg_info, prop0_r = prop_bc[[r]],
                              ub_prop = ub_prop[r], lb_prop = lb_prop[r])
    end_time <- proc.time()
    elapse_time <- end_time - start_time
    print(sprintf("MCMC for region %i finished, elapsed time %f", r, elapse_time[3]))
}
# save(res, file = sprintf("/Users/hdirector/Desktop/fit_Month%i_Train%i_%i.rda",
#                            month, train_start_year, train_end_year))

#Compute mu and sigma for each region
pars <- list()
for (r in regs_to_fit) {
	pars[[r]] <- calc_pars(res_r = res[[r]], burn_in)
}

print("fitting complete")

#-------------------------------------
#generate contours
#-------------------------------------
#generate contours
indiv_conts <- list()
for (r in regs_to_fit) {
  indiv_conts[[r]] <- gen_cont(r, pars_r = pars[[r]], reg_info, n_gen)
  print(sprintf("contours for region %i generated", r))
}

conts <- merge_conts(conts = indiv_conts, full = full)
prob_cont <- prob_map(merged = conts)
# save(hybrid_prob, file = sprintf(".../%s/hybrid_prob/hybrid_prob_month%i_year%i_train%i_%i_init%i.rda",
#                                 dyn_mod, month, forecast_year, train_start_year, train_end_year, init_month))
print("completed hybrid_prob")
print("job complete")
