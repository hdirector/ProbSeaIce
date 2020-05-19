#Script to compute the damped persistence reference code for a particular month
#and lag

#month/lag combination to compute
months <- 1
lag <- 0

#general set up
library("IceCastV2")
fYears <- 2008:2016
nX <- 304
nY <- 448

#Read in observations
allYears <- 1981:2016; nYears <- length(allYears)
obs <- read_monthly_BS(allYears[1], allYears[nYears],
                       file_folder = "/homes/direch/bootstrapV3_1/",
                       version = 3.1)
obs[obs == 120] <- NA
obs[obs == 110] <- 100 #pole hole assumed to be ice
obs <- obs/100

#Function to bound a forecasted concentration that is outside [0, 1]
bdFcast <- function(x) {
  if (x < 0) {
    x <- 0
  } else if (x > 1) {
    x <- 1
  }
  return(x)
}

#function to compute linear trend and detrend data
#' @param x time series of concentration values
trendFcast <- function(x) {
  year_ind <- 1:length(x)
  lmTemp <- lm(x ~ year_ind)
  fitted <- rep(NA, length(year_ind))
  fitted[!is.na(x)] <- lmTemp$fitted.values
  dTrend <- x - fitted
  fCast <- predict(lmTemp, newdata = data.frame(year_ind = max(year_ind) + 1))
  fCast <- bdFcast(fCast)
  return(list("fCast" = fCast, "dTrend" = dTrend))
}

#function to compute dampened persistence forecast
#' @param x_v validation time concentration values
#' @param x_i initialization time concncetration values
#' @param t forecast year
#' @param m forecast month
#' @param init_month initialization month
#' @param allYears years of data available
dPersisFcast <- function(x_v, x_i, t, m, init_month, allYears) {
  #no need to do anything if whole time seris is NA's
  if (all(is.na(x_v))) {
    return(NA)
  }

  #pull out training data
  if (init_month > m) {
    #offset of year indices and one less training year
    x_i_train <- x_i[1:max(which(allYears < t - 1))]
    x_v_train <- x_v[2:max(which(allYears < t))]
  } else {
    x_i_train <- x_i[which(allYears < t)]
    x_v_train <- x_v[which(allYears < t)]
  }

  #compute linear regression trends for validation and intiialization month
  temp <- trendFcast(x_v_train)
  trend_v <- temp$fCast
  dTrend_train_v <- temp$dTrend
  rm(temp)
  temp <- trendFcast(x_i_train)
  trend_i <- temp$fCast
  dTrend_train_i <- temp$dTrend
  rm(temp)

  #remove cases where either validation or intialization is NA
  if (any(is.na(dTrend_train_i))|any(is.na(dTrend_train_v))) {
    rm_ind <- which(is.na(dTrend_train_i) | is.na(dTrend_train_v))
    dTrend_train_v <- dTrend_train_v[-rm_ind]
    dTrend_train_i <- dTrend_train_i[-rm_ind]
  }

  #compute alpha
  if ((sd(dTrend_train_v, na.rm  = T) != 0) &
      (sd(dTrend_train_i, na.rm  = T) != 0)) {
    alpha <- cor(dTrend_train_v, dTrend_train_i)
  } else { #sd of at least one time series = 0, just use linear regression
    return(trend_v)
  }

  #compute forecast
  D <- (x_i[which(allYears == t)] - trend_i)*alpha
  dPersisFcast <- trend_v + D

  #bound if above 1 or below 0
  dPersisFcast <- bdFcast(dPersisFcast)
  return(dPersisFcast)
}

#compute damped persistence forecast
for (m in months) {
  init_month <- get_init_month(m, lag)
  for (t in fYears) {
    dPersisInit <- matrix(nrow = nX, ncol = nY)
    for (i in 1:nX) {
      for (j in 1:nY) {
        #pull out data in validation month and initialization month
        #for all years
        x_i <- obs[,init_month, i, j]
        x_v <- obs[,m, i, j]
        dPersisInit[i,j] <- dPersisFcast(x_v, x_i, t, m, init_month, allYears)
      }
    }
    dPersis <- matrix(nrow = nX, ncol = nY, data = 0)
    dPersis[dPersisInit >= .15] <- 1
    dPersis[land_mat == 1] <- NA
    save(dPersis, file = sprintf("~/probForecast/results/dPersis/dPersis_month%i_year%i_lag%i.rda",
                                 m, t, lag))
    print(c(m ,t))
  }
}

