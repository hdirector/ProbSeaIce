#-------------------------------------------------------------------------------
# Script to make Figure 1 in paper showing regions of analysis. 
#
# Requires Data: region_n.msk (region mask for the Arctic produced by the 
#                              National Snow and Ice Data Center)
#-------------------------------------------------------------------------------

library("IceCast")
#pdf("Paper/Figures/regionMap.pdf", height = 3, width = 2)
par(oma = c(2, 0, 0, 0), mar = c(0, 0, 0, 0))
#---------------------------------------
#region maps
#---------------------------------------
xBdInd <- 55:225; xBdN <- length(xBdInd)
yBdInd <- 130:310; yBdN <- length(yBdInd)
par(oma = rep(0, 4), mar = rep(0, 4))
#colors
lineCols <- rep("grey5", 5)
regCols <- c("khaki1", "lavender", "lightgreen",
             "mistyrose", "lightsalmon")


#Read in NSIDC regions
fileName <- "Data/grids/region_n.msk"
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
plot(bb, col = 'white', border = 'white')
for (i in 1:nReg) {
  plot(reg_info$regions[[i]], add = T, col = regCols[i], border = regCols[i])
}
plot(non_reg, col = "white", add = T, border = "white", )
plot(land, add = T, col = 'grey55', border = 'grey55')
plot(reg_info$start_lines[[1]], add = T, col = lineCols[1],
     lwd = 2)


for (r in 2:nReg) {
  points(reg_info$start_coords[[r]], pch = 20, cex = .6, col = lineCols[[r]])
}


text(1000, 2800, "Russia", font = 2, cex = .6, col= 'black')
text(-2600, 1000, "Alaska", font = 2, cex = .6, col= 'black')
text(-2800, -1400, "Canada", font = 2, cex = .6, col= 'black')
#dev.off()
