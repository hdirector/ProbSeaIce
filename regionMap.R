pdf("Paper/Figures/regionMap.pdf")
xBdInd <- 55:225; xBdN <- length(xBdInd) 
yBdInd <- 130:310; yBdN <- length(yBdInd)
par(oma = rep(0, 4), mar = rep(0, 4))
#colors
lineCols <- c("darkblue", "purple4", "darkgreen",
              "hotpink", "darkred")
regCols <- c("dodgerblue1", "mediumpurple1", "lightgreen",
             "lightpink", "coral")

#make map
plot(land, col = 'grey', border = "grey")
for (i in 1:nReg) {
  plot(reg_info$regions[[i]], add = T, col = regCols[i], border = regCols[i])
} 
plot(reg_info$start_lines[[1]], add = T, col = lineCols[1],
     lwd = 2)
for (i in 2:nReg) {
  plot(reg_info$start_lines[[i]], add = T, col = lineCols[i])
}

text(1000, 2800, "Russia", font = 2, cex = .65, col= 'navy')
text(-2700, 1000, "Alaska", font = 2, cex = .65, col= 'navy')
text(-2800, -1700, "Canada", font = 2, cex = .65, col= 'navy')
text(-500, 500, "Central Arctic", font = 2, cex = .65)
text(-2700, 2300, "Bering Sea", font = 2, cex = .65)
text(-600, -3000, "Baffin \n Bay", font = 2, cex = .65)
text(1250, -1500, "Greenland \n Sea", font = 2, cex = .65)
dev.off()