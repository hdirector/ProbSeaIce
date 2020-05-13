#Script to produce all reliability diagrams in the paper

#Packages
library("IceCast")
library('RcppCNPy')

#Constants
months <- 1:12
years <- 2008:2016
lags <- 1:2
nX <- 304; nY <- 448
stat_train <- 10
sip_filepath <- "/Users/hdirector/Dropbox/SeaIce_InProgress/probContours_ECMWF/Data/ecmwfsipn/forecast/ecmwfsipn_sip"
sip_start_year <- 1993
n_ens <- 25

#Load observations,
obs_file_path <-  "Data/bootstrapV3_1/"
obs <- read_monthly_BS(start_year = years[1],
                       end_year = years[length(years)],
                       version = 3.1,
                       file_folder = obs_file_path)


#identify non-regional ocean
all_regions_mask <- conv_to_grid(all_regions)
non_reg_ocean <- matrix(nrow = nX, ncol = nY, data = 0)
non_reg_ocean[all_regions_mask == 0] <- 1

#basic info
n_months <- length(months)
n_years <- length(years)
n_lags <- length(lags)

###load one observation and one prediction to identify differences in NA patterns
#between dynamic model and post-processed
load(sprintf("%s/initMonth1.rda", sip_filepath))
dyn_prob <- sip[1,1,,]
#observation for comparison
obs_temp <- obs[1, 1,,]
obs_temp[obs_temp == 120] <- NA #land index
obs_temp[obs_temp == 110] <- 100 #satellite hole is assumed to be ice
obs_curr <- matrix(nrow = nX, ncol = nY)
obs_curr[obs_temp >= .15] <- 1
obs_curr[obs_temp < .15] <- 0
obs_curr[land_mat == 1] <- NA
NA_in_dyn <- which(is.na(dyn_prob) & !is.na(obs_curr))
rm(sip)

#Function to figure out the probability bins giveen the number of ensembles
get_cat_info <- function(n_ens) {
  cat <- 0:n_ens/n_ens; n_cat <- length(cat)
  diff <- (cat[2] - cat[1])/2
  #rounding is to ensure numeric precision
  lb <- round(c(0, cat[2:n_cat] - diff), digits = 2)
  ub <- round(c(cat[1:(n_cat - 1)] + diff, 1), digits = 2)
  mean_cat <- (lb + ub)/2
  return(list("lb" = lb, "ub" = ub, "mean_cat" = mean_cat))
}

#Function to convert probabilities to prediction categories
conv_to_cat <- function(pred, n_ens, cat_info) {
  lb <- cat_info$lb; ub <- cat_info$ub
  n_cat <- length(lb)
  cat_map <- matrix(nrow = nX, ncol = nY)
  for (i in 1:n_cat) {
    if (i != 1) {
      cat_map[pred > lb[i] & pred <= ub[i]] <- i
    } else {
      cat_map[pred >= lb[i] & pred <= ub[i]] <- i
    }
  }
  stopifnot(sum(is.na(cat_map)) == sum(is.na(pred)))
  cat_map[non_reg_ocean == 1] <- NA
  return(cat_map)
}

#Function to calculate mean of observations within a category
calc_mean_cat <- function(pred, n_ens, cat_info, obs_curr) {
  obs_curr[non_reg_ocean == 1] <-  NA
  cat_map <- conv_to_cat(pred, n_ens, cat_info)
  stopifnot(sum(is.na(obs_curr)) == sum(is.na(cat_map)))
  n_tot <- sum(!is.na(cat_map))
  n_cat <- n_ens + 1
  mean_prop <- rep(NA, n_cat) #mean proportion of times ice observed for each category
  for (i in 1:n_cat) {
    curr <- which(cat_map == i)
    mean_prop[i] <- mean(obs_curr[curr])
  }
  return(mean_prop)
}

#Make results data frame
mod <- c("dyn_prob", "mcf_prob", "taqm")
calib <- expand.grid("mod" = mod, "lag" = lags, "year" = years, "month" = months)
cat_info <- get_cat_info(n_ens)
temp <- matrix(nrow = nrow(calib), ncol = length(cat_info$mean_cat))
colnames(temp) <- cat_info$mean_cat
calib <- cbind(calib, temp)
c1_prop <- which(colnames(calib) == "0.01")
c2_prop <- which(colnames(calib) == "0.99")


#region map for analysis
reg_cat <- matrix(nrow = nX, ncol = nY, data = 0)
reg_cat[land_mat == 1] <- NA
for (i in 1:18) {
  temp <- conv_to_grid(reg_info$regions[[i]])
  curr <- which(temp == 1)
  stopifnot(all(reg_cat[curr] == 0)) #confirm no overlappping pixels
  reg_cat[curr] <- i
}


#Main loop
for (y in 1:n_years) {
  for (m in 1:n_months) {

    #observation for comparison
    obs_temp <- obs[length(years[1]:years[y]), months[m],,]
    obs_temp[obs_temp == 120] <- NA #land index
    obs_temp[obs_temp == 110] <- 100 #satellite hole is ice
    obs_curr <- matrix(nrow = nX, ncol = nY)
    obs_curr[obs_temp >= .15] <- 1
    obs_curr[obs_temp < .15] <- 0
    obs_curr[land_mat == 1] <- NA
    obs_curr[NA_in_dyn] <- NA

    #identify years where training begins for stat model
    train_start_year <- years[y] - stat_train
    train_end_year <- years[y] - 1

    #climatological probabilistic "clim_prob"
    load(sprintf("Results/clim_prob/clim_prob_month%i_year%i.rda",
                 months[m], years[y], train_start_year, train_end_year))
    clim_prob[NA_in_dyn] <- NA

    ####compute calibration info for forecasts where the lag matters
    for (j in 1:n_lags) {
      init_month <- get_init_month(months[m], lags[j])

      #post-processed ensemble probabilistic forecast ("dyn_prob")
      load(sprintf("%s/initMonth%i.rda", sip_filepath, init_month))
      dyn_prob <- sip[length(sip_start_year:years[y]), months[m],,]
      dyn_ind <- which(calib$mod == "dyn_prob" & calib$lag == lags[j] &
                       calib$year == years[y] & calib$month == months[m])
      stopifnot(length(dyn_ind) == 1)
      calib[dyn_ind, c1_prop:c2_prop]<- calc_mean_cat(pred = dyn_prob, n_ens, 
                                                      cat_info, obs_curr)
      rm(dyn_prob)

      #mixture contour forecast probabilistic
      load(sprintf("Results/mcf_prob/mcf_prob_month%i_year%i_train%i_%i_init%i.rda",
                   months[m], years[y], train_start_year, train_end_year, init_month))
      mcf_ind <- which(calib$mod == "mcf_prob" & calib$lag == lags[j] &
                        calib$year == years[y] & calib$month == months[m])
      stopifnot(length(mcf_ind) == 1)
      calib[mcf_ind, c1_prop:c2_prop] <- calc_mean_cat(pred = mcf_prob, n_ens,
                                                        cat_info, obs_curr)
      rm(mcf_prob)

      #trend adjusted quantile mapping ("TAQM")
      taqm <- npyLoad(sprintf("Results/taqm/taqm_month%i_year%i_init%i.npy",
                              months[m], years[y], init_month))
      taqm_ind <- which(calib$mod == "taqm" & calib$lag == lags[j] &
                        calib$year == years[y] & calib$month == months[m])
      stopifnot(length(taqm_ind) == 1)
      calib[taqm_ind, c1_prop:c2_prop] <- calc_mean_cat(pred = taqm, n_ens, 
                                                        cat_info, obs_curr)
      rm(taqm)
    }
    print(c(m, y))
  }
}
calib$month <- as.factor(calib$month)
#save(calib, file = "Results/summaries/calib.rda")
#load("Results/summaries/calib.rda")


library("tidyverse")

#rearrange table
calib_prop <- calib[,1:c2_prop] %>% gather(cat, prop, -month, -lag, -mod, -year)

#assign long or short category
lag_lab <- c(">= 2 month lag", "<2 month lag")
calib_prop <- calib_prop %>% mutate("lag_l2" = factor(lag_lab[as.numeric(lag >= 2) + 1],
                                                      levels = c("<2 month lag", ">= 2 month lag")))

#add seaoson category
calib_prop<- calib_prop %>%  
  mutate(season = ifelse(month %in% 1:3, "Jan/Feb/Mar",
                         ifelse(month %in% 4:6, "Apr/May/Jun", 
                                ifelse(month %in% 7:9, "Jul/Aug/Sep", "Oct/Nov/Dec"))))
calib_prop$season <- factor(calib_prop$season, levels = c("Jan/Feb/Mar", 
                                                          "Apr/May/Jun", 
                                                          "Jul/Aug/Sep", 
                                                          "Oct/Nov/Dec"))

#make category numeric
calib_prop$cat <- as.numeric(calib_prop$cat)

#compute means for each category
calib_sum <- calib_prop %>%
  group_by(mod, month, cat, lag_l2, season) %>%
  summarize(mean_prop = mean(prop))

#make month index numbers into month labels
calib_sum$month <- factor(month.abb[calib_sum$month],
                          levels  =  c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                       "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))

#Formal names
calib_sum$mod <-  recode(calib_sum$mod, "dyn_prob" = "Ensemble",
                         "mcf_prob" = "Mixture Contour Forecast",
                         "taqm" = "TAQM")

#------------------------------------------------------------
#Reliability Diagrams for all months, supplement
#----------------------------------------------------------
p_l2 <- ggplot(filter(calib_sum, lag_l2 == "<2 month lag"),
             aes(x = cat, y = mean_prop, color = month)) +
  geom_point(size = .95) +
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "Shorter Lead Times",
       caption = "Probabilities grouped to nearest 1/25") +
  xlab("Forecast Probability") + ylab("Observed Proportion") +
  scale_x_continuous(breaks = c(0, .25, .5, .75, 1),
                     labels=c("0", "0.25", "0.5", "0.75", "01"))+
  scale_y_continuous(breaks = c(0, .25, .5, .75, 1),
                     labels=c("0", "0.25", "0.5", "0.75", "01")) +
  theme(strip.text = element_text(colour = 'navy', face = "bold"))

p_g2 <- ggplot(filter(calib_sum, lag_l2 == ">= 2 month lag"),
               aes(x = cat, y = mean_prop, color = month)) +
  geom_point(size = .95)  +
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "Longer Lead Times",
       caption = "Probabilities grouped to nearest 1/25") +
  xlab("Forecast Probability") + ylab("Observed Proportion") +
  scale_x_continuous(breaks = c(0, .25, .5, .75, 1),
                     labels=c("0", "0.25", "0.5", "0.75", "01")) +
  scale_y_continuous(breaks = c(0, .25, .5, .75, 1),
                     labels=c("0", "0.25", "0.5", "0.75", "01")) +
  theme(strip.text = element_text(colour = 'navy', face = "bold"))

#pdf("Paper/Figures/calib_l2.pdf", width = 8.5, height = 8.5)
p_l2 + facet_grid(rows = vars(mod), cols = vars(season), switch = "y")
#dev.off()

#pdf("Paper/Figures/calib_g2.pdf", width = 8.5, height = 8.5)
p_g2 + facet_grid(rows = vars(mod), cols = vars(season), switch = "y")
#dev.off()


#------------------------------------------------------------
#Reliability Diagrams for ASO, results section
#----------------------------------------------------------
calib_ASO <- filter(calib_sum, month == "Aug"|month == "Sep"|month == "Oct")
calib_ASO$lag_l2 <-  recode(calib_ASO$lag_l2, 
                            "<2 month lag" ="Shorter Lead Times",
                            ">= 2 month lag" = "Longer Lead Times")

p_ASO <- ggplot(filter(calib_ASO),
                   aes(x = cat, y = mean_prop, color = month)) +
  geom_point(size = .95)  +
  geom_abline(intercept = 0, slope = 1) +
  labs(caption = "Probabilities grouped to nearest 1/25") +
  xlab("Forecast Probability") + ylab("Observed Proportion") +
  scale_x_continuous(breaks = c(0, .25, .5, .75, 1),
                     labels=c("0", "0.25", "0.5", "0.75", "01")) +
  scale_y_continuous(breaks = c(0, .25, .5, .75, 1),
                     labels=c("0", "0.25", "0.5", "0.75", "01")) +
  theme(strip.text = element_text(colour = 'navy', face = "bold"))

#pdf("Paper/Figures/calib_ASO.pdf", width = 7, height = 4)
p_ASO + facet_grid(cols = vars(mod), rows= vars(lag_l2), switch = "y")
#dev.off()


#----------------------------------------------------------
#Reliability Diagrams for September, introduction section
#----------------------------------------------------------
calib_sep <- filter(calib_sum, month == "Sep" & lag_l2 == "<2 month lag" &
                     mod != "TAQM")
p_sep <- ggplot(filter(calib_sep),
                aes(x = cat, y = mean_prop, color = month)) +
  geom_point(size = .95)  +
  geom_abline(intercept = 0, slope = 1) +
  labs(caption = "Probabilities grouped to nearest 1/25") +
  xlab("Forecast Probability") + ylab("Observed Proportion") +
  scale_x_continuous(breaks = c(0, .25, .5, .75, 1),
                     labels=c("0", "0.25", "0.5", "0.75", "01")) +
  scale_y_continuous(breaks = c(0, .25, .5, .75, 1),
                     labels=c("0", "0.25", "0.5", "0.75", "01")) +
  theme(strip.text = element_text(colour = 'navy', face = "bold")) + 
  theme(legend.position = "none")
#pdf("Paper/Figures/calib_sep.pdf", height = 3, width = 5)
p_sep + facet_grid(cols = vars(mod), switch = "y")
#dev.off()



