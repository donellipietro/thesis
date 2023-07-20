
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

save(show_tests, file = "utils/functions/help_functions.RData")