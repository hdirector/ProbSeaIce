#Converting specified proportion to polygon back
# to proportion should approximately retain proportion
# Slight differences expected because of rounding to land
# at 12.5 and iterated application of Douglas-Peuker
# algorithm when merging

############# r = 5 case
set.seed(109)
r <- 5
n_lines <- length(reg_info$lines[[r]])
angs_r <- reg_info$angs[[r]]

#convert to proportions
prop_gen <- c(rep(0, 15), runif(19, .3, .8))

#convert back to y, includes land in lengths
y_gen <- prop_to_y(prop = prop_gen, r = r, reg_info = reg_info,
                   inds = 1:n_lines)

#convert to (x, y) points and SpatialPolygons object
new_pts <- reg_info$start_coords[[r]] +
  cbind(y_gen*cos(angs_r), y_gen*sin(angs_r))
new_poly <- make_polygons(r = r, my_end = new_pts)

#compute y, excluding land, from new_poly
comp_y <- find_y_1(new_poly, reg_info)[[r]]

#calculate proportions, excluding land
prop_calc <- comp_y/reg_info$max_lengths[[r]]

#confirm original proportion matches compute proportion
#(or known reason for difference)
plot(prop_gen, type = 'l')
points(prop_calc, type= 'l', col = 'red')
which(abs(prop_gen - prop_calc) > 1e-5)

############# r = 1 test case
set.seed(109)
r <- 1
n_lines <- length(reg_info$lines[[r]])
angs_r <- reg_info$angs[[r]]
#convert to proportions
prop_gen <- c(rep(1, 15), runif(57, .6, .8), rep(1, 18))

#convert back to y, includes land in lengths
y_gen <- prop_to_y(prop = prop_gen, r = r, reg_info = reg_info,
                   inds = 1:n_lines)

#convert to (x, y) points and SpatialPolygons object
new_pts <- reg_info$start_coords[[r]] +
           cbind(y_gen*cos(angs_r), y_gen*sin(angs_r))
new_poly <- make_polygons(r = r, my_end = new_pts)

#compute y, excluding land, from new_poly
comp_y <- find_y_1(new_poly, reg_info)[[r]]

#calculate proportions, excluding land
prop_calc <- comp_y/reg_info$max_lengths[[r]]

#confirm original proportion matches compute proportion
#(or known reason for difference)
plot(prop_gen, type = 'l')
points(prop_calc, type= 'l', col = 'red')
which(abs(prop_gen - prop_calc) > 1e-5)


############# r = 4 case
set.seed(109)
r <- 4
n_lines <- length(reg_info$lines[[r]])
angs_r <- reg_info$angs[[r]]

#convert to proportions
prop_gen <- c(rep(0, 10), runif(16, .4, .9), 1)

#convert back to y, includes land in lengths
y_gen <- prop_to_y(prop = prop_gen, r = r, reg_info = reg_info,
                   inds = 1:n_lines)

#convert to (x, y) points and SpatialPolygons object
new_pts <- reg_info$start_coords[[r]] +
  cbind(y_gen*cos(angs_r), y_gen*sin(angs_r))
new_poly <- make_polygons(r = r, my_end = new_pts)

#compute y, excluding land, from new_poly
comp_y <- find_y_1(new_poly, reg_info)[[r]]

#calculate proportions, excluding land
prop_calc <- comp_y/reg_info$max_lengths[[r]]

#confirm original proportion matches compute proportion
#(or known reason for difference)
plot(prop_gen, type = 'l')
points(prop_calc, type= 'l', col = 'red')
which(abs(prop_gen - prop_calc) > 1e-5)






