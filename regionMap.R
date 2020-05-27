library("IceCast")
pdf("Paper/Figures/regionMap.pdf")
xBdInd <- 55:225; xBdN <- length(xBdInd)
yBdInd <- 130:310; yBdN <- length(yBdInd)
par(oma = rep(0, 4), mar = rep(0, 4))
#colors
lineCols <- c("darkblue", "purple4", "darkgreen",
              "hotpink", "darkred")
regCols <- c("dodgerblue1", "mediumpurple1", "lightgreen",
             "lightpink", "coral")


#Read in NSIDC regions
fileName <- "/Users/hdirector/Dropbox/SeaIce_InProgress/probContours_ECMWF/Data/grids/region_n.msk"
to.read <- file(fileName, "rb")
nX <- 304; nY <- 448
to.read = file(fileName, "rb")
skip <- readBin(to.read, "raw", 300) #skip 300 byte header
dat <- readBin(to.read, "int", size = 1, signed = FALSE, n = nX*nY)
close(to.read)
reg <- matrix(dat, nrow = nX)
reg <- reg[, nY:1]
non_reg <- matrix(nrow = nX, ncol = nY, data = 0)
non_reg[reg == 1] <- 1
non_reg <- get_region(non_reg, dat_type = "simple", use_all = TRUE)
xmn = -3850; xmx = 3750; ymn = -5350; ymx = 5850 #min and max
bbPts <- rbind(c(xmn, ymn), c(xmn, ymx), c(xmx, ymx), c(xmx, ymn))
bb <- SpatialPolygons(list(Polygons(list(Polygon(bbPts)), "box")))


#make map
nReg <- length(reg_info$regions)
plot(bb, col = 'beige')
for (i in 1:nReg) {
  plot(reg_info$regions[[i]], add = T, col = regCols[i], border = regCols[i])
}
plot(non_reg, col = "white", add = T, border = "white")
plot(land, add = T, col = 'grey', border = 'grey')
plot(bb, border = 'white', add = T)
plot(reg_info$start_lines[[1]], add = T, col = lineCols[1],
     lwd = 2)
for (i in 2:nReg) {
  points(reg_info$start_lines_coords[[i]], type= 'l', col = lineCols[i])
}

text(1000, 2800, "Russia", font = 2, cex = .65, col= 'navy')
text(-2700, 1000, "Alaska", font = 2, cex = .65, col= 'navy')
text(-2800, -1700, "Canada", font = 2, cex = .65, col= 'navy')
text(-500, 500, "Central Arctic", font = 2, cex = .65)
text(-2700, 2300, "Bering Sea", font = 2, cex = .65)
text(-600, -3000, "Baffin \n Bay", font = 2, cex = .65)
text(1250, -1500, "Greenland \n Sea", font = 2, cex = .65)
text(-950, 4000, "Sea of \n Okhotsk", font = 2, cex = .65)

dev.off()
