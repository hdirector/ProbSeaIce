#-------------------------------------------------------------------------------
# Script to produce the task_table (a .rda file giving all years, months, lead 
# times, and training lengths to forecast with contour model and MCF)
#-------------------------------------------------------------------------------

n_train <- 10
forecast_years <- 2005:2016
n_years <- length(forecast_years)
lag <- 0:6
n_lags <- length(lag)
months <- 1:12
n_months <- length(months)
task_table <- data.frame("month" = rep(rep(months, each = n_years), n_lags),
                       "train_start_year" = rep(rep(forecast_years - n_train,
                                                  n_months), n_lags),
                       "train_end_year" = rep(rep(forecast_years - 1, n_months),
                                            n_lags),
                       "forecast_year" = rep(rep(forecast_years, n_months),
                                            n_lags),
                       "lag" = rep(lag, each = n_months*n_years))


task_table$n_train_years <- task_table$train_end_year - task_table$train_start_year + 1
save(task_table,
     file = "exper_design/ecmwfExper.rda")
