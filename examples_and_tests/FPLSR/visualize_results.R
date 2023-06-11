# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Utility to visualize results of FPLSR testing %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()


# ||||||||||||||||||||||||||
# Environment variables ----
# ||||||||||||||||||||||||||

# Where to get the data
SIMPLS = "_SIMPLS"
SELECTION = ""
test_name <- "2D_test4"
tests_dir <- paste("../../fdaPDE/test/data/models/FPLSR", SIMPLS ,"/",test_name,"/", sep = '')
at_locations <- TRUE
is_test0 <- FALSE

path_images <- paste("images", SIMPLS, SELECTION, "/", sep = '')



# ||||||||||||
# Results ----
# ||||||||||||

if (!file.exists(path_images)){
  dir.create(path_images)
}
path_images_test_suite <- paste(path_images, test_name, "/", sep = '')
if (!file.exists(path_images_test_suite)){
  dir.create(path_images_test_suite)
}
for(i in 1:6){
  path_images_test <- paste(path_images_test_suite, "test", i, sep = '')
  if (!file.exists(path_images_test)){
    dir.create(path_images_test)
  }
}


if(!is_test0){
  result <- read.csv(paste(tests_dir,"errors.csv", sep = ''), header=TRUE)
  print(result)
}

for(i in 1:6){
  
  path_images_test <- paste(path_images_test_suite, "test", i, "/", sep = '')
  jpeg(file=paste(path_images_test, "Y.jpg", sep = ''))
  
  par(mfrow = c(1,1))
  
  path <- paste(tests_dir, "test", i, "/", sep = '')
  
  Y_clean <- read.csv(paste(path,"Y_clean.csv", sep = ''), header = TRUE)[,2]
  plot(1:length(Y_clean), Y_clean, main = "Y-space", xlab = "samples")
  Y <- read.csv(paste(path,"Y.csv", sep = ''), header = TRUE)[,2]
  points(1:length(Y), Y, pch = 3,  col = "blue")
  if(is_test0)
    Y_hat <- as.vector(read.csv(paste(path,"Y_hat.csv", sep = ''), header = TRUE))$V1
  else
    Y_hat <- as.vector(read.csv(paste(path,"results", SELECTION, "/Y_hat.csv", sep = ''), header = FALSE))$V1
    
  points(1:length(Y_hat), as.vector(Y_hat), pch = 4,  col = "green")
  legend("topleft", c("Y_clean", "Y", "Y_hat"), pch = c(1,3,4), col = c("black", "blue", "green"))
  
  dev.off();
}

for(i in 1:6){
  
  path_images_test <- paste(path_images_test_suite, "test", i, "/", sep = '')
  jpeg(file=paste(path_images_test, "X.jpg", sep = ''))
  
  path <- paste(tests_dir, "test", i, "/", sep = '')
  
  par(mfrow = c(2,2))
  X_clean <- read.csv(paste(path, "X_clean.csv", sep = ''), header = TRUE)[,2:3601]
  image(matrix(as.numeric(X_clean[1,]), 60, 60), main = "X_clean")
  
  if(!at_locations){
    X <- read.csv(paste(path, "X.csv", sep = ''), header = TRUE)[,2:3601]
    image(matrix(as.numeric(X[1,]), 60, 60), main = "X")
  }
  
  if(is_test0)
    X_hat <- read.csv(paste(path,"X_hat.csv", sep = ''), header = TRUE)[,2:3601]
  else
    X_hat <- read.csv(paste(path,"results", SELECTION, "/X_hat.csv", sep = ''), header = FALSE)
  image(matrix(as.numeric(X_hat[1,]), 60, 60), main = paste("X_hat", sep = ''))
  
  error <- X_clean[1,] - X_hat[1,]
  image(matrix(as.numeric(error), 60, 60), main = paste("error", sep = ''))
  
  if(at_locations){
    locs <- read.csv(paste(path, "locations.csv", sep = ''), header = TRUE)[,2:3]
    plot(locs, pch = 8, col = "red", main = "Locations")
  }

  dev.off();
}

for(i in 1:6){
  
  path_images_test <- paste(path_images_test_suite, "test", i, "/", sep = '')
  jpeg(file=paste(path_images_test, "B.jpg", sep = ''))
  
  path <- paste(tests_dir, "test", i, "/", sep = '')

  par(mfrow = c(2,2))
  B <- read.csv(paste(path,"B.csv", sep = ''), header = TRUE)[,2]
  image(matrix(as.numeric(B), 60, 60), main = "B")
  
  if(is_test0)
    B_hat <- read.csv(paste(path,"B_hat.csv", sep = ''), header = TRUE)$V1
  else
    B_hat <- read.csv(paste(path,"results", SELECTION, "/B_hat.csv", sep = ''), header = FALSE)$V1
  
  image(matrix(as.numeric(B_hat), 60, 60), main = paste("B_hat", sep = ''))
  
  error <- B - B_hat
  image(matrix(as.numeric(error), 60, 60), main = paste("error", sep = ''))

  dev.off();
}

# |||||||||||||||||||||||||||||||||
# Comparison lambda for test 2 ----
# |||||||||||||||||||||||||||||||||
# 
# exp_range <- 1:12
# 
# results = list()
# for(x in exp_range){
#   results[[paste(x, sep = '')]] <- read.csv(paste(tests_dir,"errors_lambda1e-",x,".csv", sep = ''), header=TRUE)
# }
# 
# choosen <- 3
# 
# 
# jpeg(file=paste(path_images, "lambda_selection.jpg", sep = ''), quality = 100, width = 1100, height = 800, units = 'px')
# par(mfrow = c(2,2))
# 
# Y_error = NULL
# for(x in exp_range){
#   Y_error = rbind(Y_error, results[[paste(x, sep = '')]]$Y_error)
# }
# matplot(log(t(Y_error)), type = "l", lty = 1,
#         main = "MSE(Y, lambda)",
#         xlab = "Tests", col = rainbow(length((results))),
#         ylab = "log(MSE)")
# matlines(log(Y_error[choosen,]), type = "l", lty = 2, lwd = 2,
#         col = rainbow(length((results)))[choosen])
# 
# X_error = NULL
# for(x in exp_range){
#   X_error = rbind(X_error, results[[paste(x, sep = '')]]$X_error)
# }
# matplot(log(t(X_error)), type = "l", lty = 1,
#         main = "MSE(X, lambda)",
#         xlab = "Tests", col = rainbow(length((results))),
#         ylab = "log(MSE)")
# matlines(log(X_error[choosen,]), type = "l", lty = 2, lwd = 2,
#          col = rainbow(length((results)))[choosen])
# 
# 
# B_error = NULL
# for(x in exp_range){
#   B_error = rbind(B_error, results[[paste(x, sep = '')]]$B_error)
# }
# matplot(log(t(B_error)), type = "l", lty = 1,
#         main = "MSE(B, lambda)",
#         xlab = "Tests", col = rainbow(length((results))),
#         ylab = "log(MSE)")
# matlines(log(B_error[choosen,]), type = "l", lty = 2, lwd = 2,
#          col = rainbow(length((results)))[choosen])
# 
# 
# plot.new()
# legend("center", c(paste("lambda = 1e-", exp_range, sep = "")), col = rainbow(length((results))), lty = 1)
# 
# dev.off()
# 
# 
