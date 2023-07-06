load("scripts/functions/import_fdaPDE_mesh.RData")

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

show_tests <- function(nodes){
  
  par(mfrow = c(3,2), mar = c(1,1,1,1), oma = c(0,0,2,0))
  for(i in 1:5){
    temp_finer <- generate_B(nodes, i)
    data.range <- round(range(temp_finer))
    plot_field(nodes, nodes, temp_finer, data.range, paste("B - Test ", i, sep = ""), FALSE)
  }
  
  par(mfrow = c(2,1), mar = c(1,1,1,1), oma = c(0,0,2,0))
  for(i in 1:2){
    temp_finer <- generate_X(nodes, i)
    data.range <- round(range(temp_finer))
    plot_field(nodes, nodes, temp_finer, data.range, paste("X - Test ", i, sep = ""), FALSE)
  }
  
}

save(import_fdaPDE_mesh, plot_field, show_tests, file = "scripts/functions/help_functions.RData")