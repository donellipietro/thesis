# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Utility to visualize results of FPLSR testing %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# ||||||||||||||
# Libraries ----
# ||||||||||||||

library(rstudioapi)


# ||||||||||||||||
# Environment ----
# ||||||||||||||||

rm(list = ls())
graphics.off()

# script.dir <- dirname(sys.frame(1)$ofile)
# setwd(paste(script.dir, "/..", sep = ""))

print(getwd())


# ||||||||||||||
# Functions ----
# ||||||||||||||

create_directories <- function(path_images, test_name) {
  if (!file.exists("images/")){
    dir.create("images/")
  }
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
  
  return(path_images_test_suite)
}

plot_Y <- function(path_images_test_suite, tests_dir, is_test0, SELECTION) {
  for(i in 1:6){
    
    path_images_test <- paste(path_images_test_suite, "test", i, "/", sep = '')
    jpeg(file=paste(path_images_test, "Y.jpg", sep = ''))
    
    par(mfrow = c(1,1))
    
    path <- paste(tests_dir, "test", i, "/", sep = '')
    
    Y_clean <- read.csv(paste(path,"Y_clean.csv", sep = ''), header = TRUE)[,2]
    
    Y <- read.csv(paste(path,"Y.csv", sep = ''), header = TRUE)[,2]
    plot(Y_clean, Y, main = "Y-space", xlab = "Y_clean", ylab = "Y")
    if(is_test0)
      Y_hat <- as.vector(read.csv(paste(path,"Y_hat.csv", sep = ''), header = TRUE))$V1
    else
      Y_hat <- as.vector(read.csv(paste(path,"results", SELECTION, "/Y_hat.csv", sep = ''), header = FALSE))$V1
    points(Y_clean, as.vector(Y_hat), pch = 4,  col = "blue")
    abline(a = 0, b = 1, col = "red", lty = 2)
    
    legend("topleft", c("Y", "Y_hat"), pch = c(1,4), col = c("black", "blue"))
    
    dev.off();
  }
}


plot_X <- function(path_images_test_suite, tests_dir, at_locations, is_test0, SELECTION) {
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
}

plot_B <- function(path_images_test_suite, tests_dir, is_test0, SELECTION) {
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
}


# ||||||||||||
# Results ----
# ||||||||||||

# ||||||||||||||||||||||||||
# Environment variables ----
# ||||||||||||||||||||||||||

# Methods
METHODS_vec <- c("NIPALS", "SIMPLS")
METHODS_PATH_test_vec <- c("", "_SIMPLS")
METHODS_PATH_results_vec <- c("_NIPALS", "_SIMPLS")
SELECTION_vec <- list(c("fixed", "KCV"), c("fixed", "GCV"))
SELECTION_PATH_test_vec <- list(c("", "_KCV"), c("", "_GCV"))
SELECTION_PATH_results_vec <- list(c("_fixed", "_KCV"), c("_fixed", "_GCV"))

# Tests
test_name_vec <- c("2D_test0", "2D_test1", "2D_test2", "2D_test3", "2D_test4")
at_locations_vec <- c(FALSE, FALSE, TRUE, TRUE, TRUE)
is_test0_vec <- c(TRUE, FALSE, FALSE, FALSE, FALSE) 


ACTIVE <- list()

# NIPALS
ACTIVE[[1]] <- matrix(c(TRUE,  TRUE,  TRUE,  TRUE,  TRUE,
                        FALSE, TRUE,  FALSE, FALSE, FALSE),
                        nrow = length(SELECTION_vec[[1]]),
                        ncol = length(test_name_vec),
                        byrow = TRUE)

# SIMPLS 
ACTIVE[[2]] <- matrix(c(TRUE,  TRUE,  TRUE,  TRUE,  TRUE,
                        FALSE, TRUE,  FALSE, FALSE, FALSE),
                        nrow = length(SELECTION_vec[[2]]),
                        ncol = length(test_name_vec),
                        byrow = TRUE)

print(" ")
print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print("% FPLSR tests results visualization %")
print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print(" ")

for(m in 1:length(METHODS_vec)){
  
  METHOD <- METHODS_vec[m]
  METHOD_PATH_test <- METHODS_PATH_test_vec[m]
  METHOD_PATH_results <- METHODS_PATH_results_vec[m]
  
  for(s in 1:length(SELECTION_vec[[m]])){
    
    SELECTION <- SELECTION_vec[[m]][s]
    SELECTION_PATH_test <- SELECTION_PATH_test_vec[[m]][s]
    SELECTION_PATH_results <- SELECTION_PATH_results_vec[[m]][s]
    
    print("###############")
    print(paste(METHOD, SELECTION, sep = " "))
    print("###############")
  
    for(t in 1:length(test_name_vec)){
      
      if(ACTIVE[[m]][s,t]){
      
        # Test
        test_name <- test_name_vec[t]
        tests_dir <- paste("../../fdaPDE/test/data/models/FPLSR", METHOD_PATH_test ,"/",test_name,"/", sep = '')
        
        print(test_name)
        
        # Options
        at_locations <- at_locations_vec[t]
        is_test0 <- is_test0_vec[t]
        
        # Directories
        path_images <- paste("images/images", METHOD_PATH_results, SELECTION_PATH_results, "/", sep = '')
        path_images_test_suite <- create_directories(path_images, test_name)
        
        time_last_result <- max(tail(file.info(list.files(tests_dir, full.names = TRUE, recursive = TRUE))$ctime))
        time_last_image <- max(tail(file.info(list.files(path_images_test_suite, full.names = TRUE, recursive = TRUE))$ctime))
        
        # Results
        if(time_last_result > time_last_image){
          plot_Y(path_images_test_suite, tests_dir, is_test0, SELECTION_PATH_test)
          plot_X(path_images_test_suite, tests_dir, at_locations, is_test0, SELECTION_PATH_test)
          plot_B(path_images_test_suite, tests_dir, is_test0, SELECTION_PATH_test)
        }
        
      }
    }
    print("")
  }
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
