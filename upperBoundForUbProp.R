start_year <- 1981
end_year <- 1994

#load data
all <- IceCast::read_monthly_BS(start_year, end_year,
                                "/Users/hdirector/Dropbox/SeaIce_InProgress/probContours_ECMWF/Data/bootstrapV3_1/",
                                version = 3.1)
march <- apply(all[,3,,], 1, function(x){get_region(x, dat_type = "bootstrap",
                                                    level = 15)})
#save(march, file = "/users/hdirector/desktop/march.rda")
load(file = "/users/hdirector/desktop/march.rda")
y <- lapply(march, function(x){find_y_1(ice = x, reg_info)})
props <- lapply(y, function(x){y_to_prop(x, regs_to_fit = 1:5, reg_info)})
prop_max <- sapply(props, function(x){sapply(x, max)})
ub_prop <- apply(prop_max, 1, max)
ub_prop[ub_prop >= 1] <- 1 - eps
