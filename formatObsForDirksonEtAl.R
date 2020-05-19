library("IceCastV2")
obs <- read_monthly_BS(1993, 2017, 
              file_folder = "/Users/hdirector/Documents/SeaIce/SeaIce/bootstrap/bootstrapV3_1/",
              version = 3.1)
obs[obs == 120] <- NA
obs[obs == 110] <- 100
obs <- obs/100

save(obs, file = "/Users/hdirector/Documents/SeaIce/SeaIce/bootstrap/bootstrapV3_1/bootstrapObs1993_2017.rda")

