#Script to compute 10-year climatology sea ice probability (proportion of times
#in last ten years grid box contained sea ice)


library("IceCastV2")

#general set up
months <- 1:12
fYears <- 2005:2016
nX <- 304
nY <- 448

#Read in observations
allYears <- 1981:2016; nYears <- length(allYears)
obs <- read_monthly_BS(allYears[1], allYears[nYears],
                       file_folder = "/Users/hdirector/Dropbox/SeaIce_InProgress/probContours_ECMWF/Data/bootstrapV3_1/",
                       version = 3.1)
obs[obs == 120] <- NA
obs[obs == 110] <- 100
obs <- obs/100


#compute climatology forecast (proportion years in last ten years ice covered)
for (m in months) {
  for (t in fYears) {
    curr_ind <- t - allYears[1]

    #probabilistic climatology
    sic <- obs[(curr_ind - 9):curr_ind, m,,] #last ten years' observations
    use <- which(apply(sic, 1, function(x){any(!is.na(x))})) #only consider years with observed data
    clim_prob <- apply(sic[use,,], 2:3, function(x){sum(x >= .15)/length(use)})
    clim_prob[land_mat == 1] <- NA
    save(clim_prob, file = sprintf("/Users/hdirector/Dropbox/SeaIce_InProgress/probContours_ECMWF/Results/clim_prob/clim_prob_month%i_year%i.rda",
                                   m, t))

    #binary climatology
    clim_bin <- matrix(nrow = nX, ncol = nY, data = 0)
    clim_bin[is.na(clim_prob)] <- NA
    clim_bin[clim_prob >= .5] <- 1
    save(clim_bin, file = sprintf("/Users/hdirector/Dropbox/SeaIce_InProgress/probContours_ECMWF/Results/clim_bin/clim_bin_month%i_year%i.rda",
                                  m, t))
    print(c(m, t))
  }
}




