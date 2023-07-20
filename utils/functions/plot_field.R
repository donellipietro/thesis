plot_field <- function(nodes, locations, data, data.range, title, UNIQUE = TRUE, ONLY_KEY = FALSE) {
  
  if(!ONLY_KEY){
    cols <- hcl.colors(100, "YlOrRd", rev = FALSE)
    
    if(UNIQUE)
      par(mfrow = c(1,1), mar = c(3, 0, 2, 0), oma = c(3,1,1,1))
    
    xlim = range(nodes[,1])
    ylim = range(nodes[,2])
    
    xlim[1] = floor(xlim[1])
    ylim[1] = floor(ylim[1])
    
    xlim[2] = ceiling(xlim[2])
    ylim[2] = ceiling(ylim[2])
    
    data.range[1] = floor(data.range[1])
    data.range[1] = floor(data.range[1])
    
    plot(NA, xlim = xlim, ylim = ylim, asp = 1,
         main = title, bty = "n", xaxt = "n", yaxt = "n",
         ylab = "", xlab = "")
    points(nodes, cex = 0.8, pch = 20, col = "grey") 
    scatter2D(x = locations[, 1], y = locations[, 2], colvar = data, col = cols,
              pch = 16, colkey = FALSE, add = TRUE)
    
    if(UNIQUE)
      colkey(clim = data.range, col = cols, side = 1, add = TRUE, 
             width = 0.5, length = 0.6,
             dist = -0.1)
  }
  else{
    colkey(clim = data.range, col = cols, side = 2, add = FALSE, 
           width = 0.5, length = 0.5,
           dist = 0)
  }
  
}

save(plot_field, file = "utils/functions/plot_field.RData")
