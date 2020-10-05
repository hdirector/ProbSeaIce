#-------------------------------------------------------------------------------
# Script to make Figure 2 in paper diagramming the components of the contour 
# model (lines, boundary points, etc.) 
#-------------------------------------------------------------------------------

library("IceCast")
nReg <- length(reg_info$regions)

#----------------------------------------
#make diagram of contour model components
#----------------------------------------
#pdf("Paper/Figures/boundLines.pdf",height = 4.5, width = 8.5)
layout(matrix(c(1, 2, 2), byrow = T, ncol = 3))
#par(oma = c(2,1,1,1), mar = c(2, 1, 1, 1))

#Bering sea 
#region to plot
bbBerPts <- rbind(c(-3600, 1000), c(-3600, 3900), c(-1300, 3900), c(-1300, 1000))
bbBer <- SpatialPolygons(list(Polygons(list(Polygon(bbBerPts)), "bbBer")))

#make up a contour as a demo
angsR <- reg_info$angs[[2]]
yLength <- rep(rep(c(.75, .65, .55), each = 4), 2)[1:20]*reg_info$max_lengths[[2]]
ep <- reg_info$start_coords[[2]] + cbind(yLength*cos(angsR), yLength*sin(angsR))
newReg <- make_polygons(r = 2, my_end = ep)

#plot everything
plot(bbBer, col = "grey5", border = "black")
mtext("Bering Sea Region", side = 3)
plot(gIntersection(bbBer, land), add = T, col = "grey55", border = "grey55")
plot(reg_info$regions[[2]], add = T, lwd = 1, col = 'white', border = 'white')
plot(rm_holes(newReg), col = 'lightblue', add = T, lwd = 2, border = "lightblue", 
     lty = 2)
plot(rm_holes(newReg), add = T, lwd = 2, border = "blue", lty = 2)
for (i in 1:length(reg_info$lines[[2]])) {
  plot(reg_info$lines[[2]][[i]], add = T, col = "mediumpurple1")
}
plot(bbBer, add = T)
points(reg_info$start_coords[[2]], col = 'purple4', pch = 20, cex = 1.5, 
       lwd = 2)
points(rbind(reg_info$start_coords[[2]][14,], ep[14,]), type = "l",
       col = "darkgreen", lwd = 3)


legend(-3500, 3800, ncol = 2, cex = 1.8, bg = "white", text.font = 2,
       legend = c(expression('y'[14]), "B"), pch = c(NA, 20), lty = c(1, NA),
       col = c("darkgreen", 'purple4'), lwd = c(5, NA), pt.cex = c(NA, 1.5))


##Central Arctic
#region to plot
bbCAPts <- rbind(c(-2900, -1200), c(2400, -1200), c(2400, 2100), c(-2900, 2100))
bbCA <- SpatialPolygons(list(Polygons(list(Polygon(bbCAPts)), "bbCA")))

#make up a contour as a demo
angsR <- reg_info$angs[[1]]
tot_lengths <- unlist(lapply(reg_info$sec_lengths[[1]], sum))
yLength <- c(rep(1, 15), seq(.9, .5, length = 7), rep(.6, 9), rep(.6, 37),
             seq(.75, 1, length = 6), rep(1, 11), rep(1, 5))*tot_lengths
ep <- reg_info$start_coords[[1]] + cbind(yLength*cos(angsR), yLength*sin(angsR))
newPts <- interp_new_pts(r = 1, new_pts = ep, reg_info = reg_info)
newReg <- SpatialPolygons(list(Polygons(list(Polygon(newPts)), "newCent")))
if (suppressWarnings(!gIsValid(newReg))) {
  newReg <- untwist(newReg)
  newReg <- gDifference(newReg, land)
}




#plot everything
plot(bbCA, col = "grey5", border = "white")
mtext("Central Arctic Region", side = 3)
plot(gIntersection(bbCA, land), add = T, col = "grey55", border = "grey55")
plot(reg_info$regions[[1]], add = T, lwd = 1, col = 'white', border = 'white')

plot(rm_holes(newReg), col = 'lightblue', add = T, lwd = 2, border = "lightblue", 
     lty = 2)
plot(rm_holes(newReg), add = T, lwd = 2, border = "blue", lty = 2)
plot(keep_poly(gIntersection(land, rm_holes(bbCA))), add = T, col = 'grey55', 
     border = 'grey55')

for (i in 1:length(reg_info$lines[[1]])) {
  plot(reg_info$lines[[1]][[i]], add = T, col = "mediumpurple")
}
points(rbind(reg_info$start_coords[[1]][33,], ep[33,]), type = "l",
       col = "darkgreen", lwd = 3)
points(reg_info$start_coords[[1]], col = 'purple4', pch = 3, cex = 2, lwd = 2)
plot(rm_holes(newReg), border = "blue", add = T, lwd = 2, lty = 2)
plot(bbCA, border = 'black', add = T)
legend(-2800,2000, ncol = 2, cex = 1.8, bg = "white",
       legend = c(expression('y'[33]), "B"), pch = c(NA, 3), lty = c(1, NA),
       col = c("darkgreen", 'purple4'),lwd = c(5, 2), pt.cex = c(NA, 1.5), 
       text.font = 2)
par(xpd = TRUE)


legend(-4500, -1350, ncol = 3, cex = 2,
       legend = c("S", "L", "ice", "not ice", "land", "outside region"),
       pch = c(NA, NA, 15, 15, 15, 15), lty = c(2, 1, NA, NA, NA, NA), 
       col = c("blue", "mediumpurple", NA, NA, NA, NA),
       lwd = c(2, 1, NA, NA, NA, NA), pt.cex = c(NA, NA, 2, 2, 2, 2), 
       border = c("white", "white", rep("black", 4)), text.font = 2,
       fill = c(NA, NA,  "lightblue", "white", "grey55", "grey5"),
       xpd = NA, text.width = 800, bty = "n")

#dev.off()




