# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Utility to compare quantitatively the results of fPLSR with Harold %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()

# ||||||||||||||
# Functions ----
# ||||||||||||||

plot_results <- function(n_test_options, tests_to_be_shown, n_tests, errors_Y, errors_X, errors_B, TEST, type, colors, names) {
  
  cols_to_be_shown <- matrix(rep(tests_to_be_shown, n_tests), ncol = length(tests_to_be_shown), byrow = TRUE)
  cols_to_be_shown <- as.vector(t(cols_to_be_shown + 0:(n_tests-1)*n_test_options))
  
  #jpeg(file=paste(path_images, "comparison_Y", type,".jpg", sep = ''))
  par(mfrow = c(1,1))
  boxplot(errors_Y[, cols_to_be_shown], xlab = "Tests", ylab = "MSE",
          col = colors[tests_to_be_shown],
          main = paste("MSE", type, "Y"), xaxt = 'n')
  axis(1, at=(1:n_tests)*(length(tests_to_be_shown))-(length(tests_to_be_shown)-1)/2, labels=1:n_tests)
  abline(v = (0:n_tests)*(length(tests_to_be_shown))+1/2, lty = 2, col = "grey")
  grid()
  #dev.off()
  
  jpeg(file=paste(path_images, "comparison_X", type,".jpg", sep = ''))
  boxplot(errors_X[,cols_to_be_shown], xlab = "Tests", ylab = "MSE",
          col = colors[tests_to_be_shown],
          main = paste("MSE", type, "X"), xaxt = 'n')
  axis(1, at=(1:n_tests)*(length(tests_to_be_shown))-(length(tests_to_be_shown)-1)/2, labels=1:n_tests)
  abline(v = (0:n_tests)*(length(tests_to_be_shown))+1/2, lty = 2, col = "grey")
  grid()
  #dev.off()
  
  #jpeg(file=paste(path_images, "comparison_X_zoom", type,".jpg", sep = ''))
  boxplot(errors_X[,cols_to_be_shown], xlab = "Tests", ylab = "MSE",
          col = colors[tests_to_be_shown], 
          main = "MSE X (zoom)", xaxt = 'n',
          ylim = c(0.003, 0.004))
  axis(1, at=(1:n_tests)*(length(tests_to_be_shown))-(length(tests_to_be_shown)-1)/2, labels=1:n_tests)
  abline(v = (0:n_tests)*(length(tests_to_be_shown))+1/2, lty = 2, col = "grey")
  grid()
  #dev.off()
  
  if(TEST){
    #jpeg(file=paste(path_images, "comparison_B", type,".jpg", sep = ''))
    par(mfrow = c(1,1))
    boxplot((errors_B[,cols_to_be_shown]), xlab = "Tests", ylab = "MSE",
            col = colors[tests_to_be_shown], 
            main = "MSE B", xaxt = 'n')
    axis(1, at=(1:n_tests)*(length(tests_to_be_shown))-(length(tests_to_be_shown)-1)/2, labels=1:n_tests)
    abline(v = (0:n_tests)*(length(tests_to_be_shown))+1/2, lty = 2, col = "grey")
    grid()
    #dev.off()
  }
  
  plot.new()
  legend("center", names[tests_to_be_shown], col = colors[tests_to_be_shown], pch = 15)
}

plot_results_divided <- function(n_test_options, tests_to_be_shown, n_tests, errors_Y, errors_X, errors_B, TEST, type, colors, names, path_images) {
  
  cols_to_be_shown <- matrix(rep(tests_to_be_shown, n_tests), ncol = length(tests_to_be_shown), byrow = TRUE)
  cols_to_be_shown <- cols_to_be_shown + 0:(n_tests-1)*n_test_options
  
  type = paste('_', type, sep = '')
  
  jpeg(file=paste(path_images, "comparison_Y", type,".jpg", sep = ''), quality = 100, width = 1100, height = 800, units = 'px')
  par(mfrow = c(2,3))
  for(i in 1:dim(cols_to_be_shown)[1]){
    boxplot(errors_Y[, cols_to_be_shown[i,]], ylab = "MSE", # xlab = "Test",
            col = colors[tests_to_be_shown],
            main = paste("MSE", type, "Y, Test", i), xaxt = 'n')
    # axis(1, at=(length(tests_to_be_shown))-(length(tests_to_be_shown)-1)/2, labels=i)
    grid()
  }
  dev.off()
  
  jpeg(file=paste(path_images, "comparison_X", type,".jpg", sep = ''), quality = 100, width = 1100, height = 800, units = 'px')
  par(mfrow = c(2,3))
  for(i in 1:dim(cols_to_be_shown)[1]){
    boxplot(errors_X[, cols_to_be_shown[i,]],  ylab = "MSE", # xlab = "Test",
            col = colors[tests_to_be_shown],
            main = paste("MSE", type, "X, Test", i), xaxt = 'n')
    # axis(1, at=(length(tests_to_be_shown))-(length(tests_to_be_shown)-1)/2, labels=i)
    grid()
  }
  dev.off()
  
  if(TEST){
    jpeg(file=paste(path_images, "comparison_B", type,".jpg", sep = ''), quality = 100, width = 1100, height = 800, units = 'px')
    par(mfrow = c(2,3))
    for(i in 1:dim(cols_to_be_shown)[1]){
      boxplot(errors_B[, cols_to_be_shown[i,]], ylab = "MSE", # xlab = "Test",
              col = colors[tests_to_be_shown],
              main = paste("MSE", type, "B, Test", i), xaxt = 'n')
      # axis(1, at=(length(tests_to_be_shown))-(length(tests_to_be_shown)-1)/2, labels=i)
      grid()
    }
    dev.off()
  }
  
  jpeg(file=paste(path_images, "comparison_legend.jpg", sep = ''), quality = 100, width = 1100, height = 800, units = 'px')
  par(mfrow = c(1,1))
  plot.new()
  legend("center", names[tests_to_be_shown], col = colors[tests_to_be_shown], pch = 15)
  dev.off()
}

merge_errors <- function(errors1, errors2, n_test, n_test_options_1, n_test_options_2){
  n_batches <- min(dim(errors1)[1], dim(errors2)[1])
  n_test_options <- n_test_options_1 + n_test_options_2
  errors_merged <- matrix(0, nrow = n_batches, ncol = 0)
  for(i in 1:n_tests){
    errors_merged <- cbind(errors_merged, errors1[1:n_batches, n_test_options_1*(i-1) + 1:n_test_options_1])
    errors_merged <- cbind(errors_merged, errors2[1:n_batches, n_test_options_2*(i-1) + 1:n_test_options_2])
  }
  
  return(errors_merged)
}

create_directories <- function(METHOD) {
  path_images <- "images/"
  if (!file.exists(path_images)){
    dir.create(path_images)
  }
  path_images <- paste("images/images",METHOD,"_comparison/", sep ="")
  if (!file.exists(path_images)){
    dir.create(path_images)
  }
  return(path_images)
}


plot_testVStrain <- function(path_images, temp_Y_train, temp_Y_test, temp_X_train, temp_X_test, functional_methods, multivariate_methods, colors) {
  jpeg(file=paste(path_images, "trainVStest",".jpg", sep = ''), quality = 100, width = 1100, height = 800, units = 'px')
  par(mfrow = c(2,2))
  y_min <- log(min(min(temp_Y_train), min(temp_Y_test)))
  y_max <- log(max(max(temp_Y_train), max(temp_Y_test)))
  matplot(log(temp_Y_test[,functional_methods]), type = "l", lty = 1,
          main = "Functional methods", xlab = "Test", ylab = "log(mean MSE Y)",
          ylim = c(y_min, y_max), col = colors[functional_methods])
  matplot(log(temp_Y_train[,functional_methods]), type = "l", lty = 2,  add = TRUE,
          col = colors[functional_methods])
  legend("topright", c("Test", "Train"), lty = c(1, 2))
  grid()
  matplot(log(temp_Y_test[,multivariate_methods]), type = "l", lty = 1,
          main = "Multivariate methods", xlab = "Test", ylab = "log(mean MSE Y)",
          ylim = c(y_min, y_max), col = colors[multivariate_methods])
  matplot(log(temp_Y_train[,multivariate_methods]), type = "l", lty = 2,  add = TRUE,
          col = colors[multivariate_methods])
  legend("topright", c("Test", "Train"), lty = c(1, 2))
  grid()
  
  y_min <- log(min(min(temp_X_train), min(temp_X_test)))
  y_max <- log(max(max(temp_X_train), max(temp_X_test)))
  matplot(log(temp_X_test[,functional_methods]), type = "l", lty = 1,
          main = "Functional methods", xlab = "Test", ylab = "log(mean MSE X)",
          ylim = c(y_min, y_max), col = colors[functional_methods])
  matplot(log(temp_X_train[,functional_methods]), type = "l", lty = 2,  add = TRUE,
          col = colors[functional_methods])
  legend("topright", c("Test", "Train"), lty = c(1, 2))
  grid()
  matplot(log(temp_X_test[,multivariate_methods]), type = "l", lty = 1,
          main = "Multivariate methods", xlab = "Test", ylab = "log(mean MSE X)",
          ylim = c(y_min, y_max), col = colors[multivariate_methods])
  matplot(log(temp_X_train[,multivariate_methods]), type = "l", lty = 2,  add = TRUE,
          col = colors[multivariate_methods])
  legend("topright", c("Test", "Train"), lty = c(1, 2))
  grid()
  dev.off()
}



# ||||||||||||||||||||||
# Visualize results ----
# ||||||||||||||||||||||

# Methods
METHODS_vec <- c("NIPALS", "SIMPLS")
METHODS_PATH_test_vec <- c("", "_SIMPLS")
METHODS_PATH_results_vec <- c("_NIPALS", "_SIMPLS")

# Test options
test_name_options_vec <- list(c("hcpp_l0", "ns_l0", "hcpp", "ns", "sr", "sri", "hcpp-KCV", "sri-GCV"),
                              c("ns_l0", "ns", "sr", "sri"))
n_batches_vec <- c(20, 20)
n_tests_vec <- c(6, 6)

# Visualization
tests_to_be_shown_vec <- list(
  # NIPALS
  c(3, 4:6, 7:8, 1, 2, 9),
  # SIMPLS
  c(2:3, 1, 5, 7)
)
colors_vec <- list(
  # NIPALS
  c("red", "orange", # 1, 2
    "darkgrey", # 3
    "purple", "blue", "lightblue", # 4, 5, 6
    "darkred", "darkblue", # 7, 8
    "darkgreen", "green", "lightgreen"), # 9, 10, 11
  # SIMPLS
  c("orange", # 1
    "purple", "blue", "lightblue", # 2, 3, 4
    "darkgreen", "green", "lightgreen") # 5, 6, 7
)
names_vec <- list(
  # NIPALS
  c("Harold C++, lambda = 0", "fPLSR no smoothing, lambda = 0",
    "Harold C++",
    "fPLSR no smoothing", "fPLSR smoothing regression", "fPLSR smoothing reg. + init.",
    "Harold C++, KCV", "fPLSR smoothing reg. + init., GCV",
    "M-PLSR (NIPALS)", "M-PLSR (NIPALS no Y deflation)", "M_PLSR (SIMPLS)"),
  # SIMPLS
  c("fPLSR_SIMPLS no smoothing, lambda = 0",
    "fPLSR_SIMPLS no smoothing", "fPLSR_SIMPLS smoothing reg.", "fPLSR_SIMPLS smoothing reg. + init.", 
    "M-PLSR (NIPALS)", "M-PLSR (NIPALS no Y deflation)", "M_PLSR (SIMPLS)")
)
functional_methods_vect <- list(
  c(3:6,7:8),
  c(2:3)
)
multivariate_methods_vect <- list(
  c(1:2,9:11),
  c(1,5:7)
)


for(m in 1:length(METHODS_vec)){
  
  # Method
  METHOD <- METHODS_vec[m]
  METHOD_PATH_test <- METHODS_PATH_test_vec[m]
  METHOD_PATH_results <- METHODS_PATH_results_vec[m]
  
  # Test dimensions
  n_batches <- n_batches_vec[m]
  n_tests <- n_tests_vec[m]
  
  # Test options
  test_name_options <- test_name_options_vec[[m]]
  n_test_options <- length(test_name_options)
  
  # Visualization options
  tests_to_be_shown <- tests_to_be_shown_vec[[m]]
  colors <- colors_vec[[m]]
  names <- names_vec[[m]]
  
  # Test vs Train options
  functional_methods <- functional_methods_vect[[m]]
  multivariate_methods <- multivariate_methods_vect[[m]]
  
  # Directories
  tests_dir = paste("../../fdaPDE/test/data/models/FPLSR",METHOD_PATH_test,"/2D_test_comparison/", sep ="")
  path_images <- create_directories(METHOD_PATH_results)
    
  # Load errors
  errors_Y <- read.csv(paste(tests_dir, "errors_Y.csv", sep = ''), header = FALSE)
  errors_X <- read.csv(paste(tests_dir, "errors_X.csv", sep = ''), header = FALSE)
  errors_B <- read.csv(paste(tests_dir, "errors_B.csv", sep = ''), header = FALSE)
  errors_Y_multivariate <- read.csv(paste(tests_dir, "errors_Y_multivariate.csv", sep = ''), header = TRUE)[,-1]
  errors_X_multivariate <- read.csv(paste(tests_dir, "errors_X_multivariate.csv", sep = ''), header = TRUE)[,-1]
  errors_B_multivariate <- read.csv(paste(tests_dir, "errors_B_multivariate.csv", sep = ''), header = TRUE)[,-1]
  errors_Y_train <- read.csv(paste(tests_dir, "errors_Y_train.csv", sep = ''), header = FALSE)
  errors_X_train <- read.csv(paste(tests_dir, "errors_X_train.csv", sep = ''), header = FALSE)
  errors_Y_train_multivariate <- read.csv(paste(tests_dir, "errors_Y_train_multivariate.csv", sep = ''), header = TRUE)[,-1]
  errors_X_train_multivariate <- read.csv(paste(tests_dir, "errors_X_train_multivariate.csv", sep = ''), header = TRUE)[,-1]
  
  # Merge dataframes
  n_test_option_1 <- dim(errors_Y)[2]/n_tests
  n_test_option_2 <- dim(errors_Y_multivariate)[2]/n_tests
  errors_Y_merged_test = merge_errors(errors_Y, errors_Y_multivariate, n_tests, n_test_option_1, n_test_option_2)
  errors_X_merged_test = merge_errors(errors_X, errors_X_multivariate, n_tests, n_test_option_1, n_test_option_2)
  errors_B_merged_test = merge_errors(errors_B, errors_B_multivariate, n_tests, n_test_option_1, n_test_option_2)
  errors_Y_merged_train = merge_errors(errors_Y_train, errors_Y_train_multivariate, n_tests, n_test_option_1, n_test_option_2)
  errors_X_merged_train = merge_errors(errors_X_train, errors_X_train_multivariate, n_tests, n_test_option_1, n_test_option_2)
  
  # Plot
  n_test_options_merged <- dim(errors_Y_merged_test)[2]/n_tests
  # Test
  plot_results_divided(n_test_options_merged, tests_to_be_shown, n_tests, errors_Y_merged_test, errors_X_merged_test, errors_B_merged_test, TRUE, "test", colors, names, path_images) 
  # Train
  plot_results_divided(n_test_options_merged, tests_to_be_shown, n_tests, errors_Y_merged_train, errors_X_merged_train, NULL, FALSE, "train", colors, names, path_images)
  
  # Comparison 
  temp_Y_test = matrix(colSums(errors_Y_merged_test)/n_batches, ncol = n_test_options_merged, byrow = TRUE)
  temp_Y_train = matrix(colSums(errors_Y_merged_train)/n_batches, ncol = n_test_options_merged, byrow = TRUE)
  temp_X_test = matrix(colSums(errors_X_merged_test)/n_batches, ncol = n_test_options_merged, byrow = TRUE)
  temp_X_train = matrix(colSums(errors_X_merged_train)/n_batches, ncol = n_test_options_merged, byrow = TRUE)
  plot_testVStrain(path_images, temp_Y_train, temp_Y_test, temp_X_train, temp_X_test, functional_methods, multivariate_methods, colors) 
  
}







































# ||||||||||||||||||
# Test vs train ----
# ||||||||||||||||||







# |||||||||||||||||||
# Compute errors ----
# |||||||||||||||||||

# for(i in 1:n_tests){
#     
#   path <- paste(tests_dir, "test", i, "/", sep = '')
#     
#   B <- read.csv(paste(path,"B.csv", sep = ''), header = TRUE)[,2]
#   n_nodes = length(B)
#     
#   for(j in 1:n_batches){
#       
#     if(TEST)
#       Y_clean <- read.csv(paste(path,"Y_clean_", j, ".csv", sep = ''), header = TRUE)[,2]
#     else
#       Y_clean <- read.csv(paste(path,"Y_", j, ".csv", sep = ''), header = TRUE)[,2]
#     
#     batch_size = length(Y_clean)
#     Y_hat_harold <- read.csv(paste(path,"Y_hat_", j, ".csv", sep = ''), header = TRUE)[,2]
#     errors_Y[j, (n_test_options+1)*i - n_test_options] <- sum((Y_clean-Y_hat_harold)^2)/batch_size
#     
#     if(TEST)
#       X_clean <- read.csv(paste(path,"X_clean_", j, ".csv", sep = ''), header = TRUE)[,-1]
#     else
#       X_clean <- read.csv(paste(path,"X_", j, ".csv", sep = ''), header = TRUE)[,-1]
#     
#     X_hat_harold <- read.csv(paste(path,"X_hat_", j, ".csv", sep = ''), header = TRUE)[,-1]
#     errors_X[j, (n_test_options+1)*i - n_test_options] <- sum((X_clean-X_hat_harold)^2)/(n_nodes*batch_size)
#     
#     if(TEST){
#       B_hat_harold <- read.csv(paste(path,"B_hat_", j, ".csv", sep = ''), header = TRUE)[,2]
#       errors_B[j, (n_test_options+1)*i - n_test_options] <- sum((B-B_hat_harold)^2)/n_nodes
#       # errors_B_min[j, (n_test_options+1)*i - n_test_options] <- min(abs(B-B_hat_harold))
#       # errors_B_max[j, (n_test_options+1)*i - n_test_options] <- max(abs(B-B_hat_harold))
#     }
#       
#     for(t in 1:n_test_options){
# 
#       test_name_option <- test_name_option_vect[t]
# 
#       # Y
#       Y_hat_test <- read.csv(paste(path,"results/Y_hat_", test_name_option, "_", j, ".csv", sep = ''), header = FALSE)$V1
#       errors_Y[j, (n_test_options+1)*i - n_test_options + t] <- sum((Y_clean-Y_hat_test)^2)/batch_size
#       
#       # X
#       X_hat_test <- read.csv(paste(path,"results/X_hat_", test_name_option, "_", j, ".csv", sep = ''), header = FALSE)
#       errors_X[j, (n_test_options+1)*i - n_test_options + t] <- sum((X_clean-X_hat_test)^2)/(n_nodes*batch_size)
# 
#       # B
#       if(TEST){
#         B_hat_test <- read.csv(paste(path,"results/B_hat_", test_name_option, "_", j, ".csv", sep = ''), header = FALSE)[,1]
#         errors_B[j, (n_test_options+1)*i - n_test_options + t] <- sum((B-B_hat_test)^2)/n_nodes
#         # errors_B_min[j, (n_test_options+1)*i- n_test_options + t] <- min(abs(B-B_hat_test))
#         # errors_B_max[j, (n_test_options+1)*i- n_test_options + t] <- max(abs(B-B_hat_test))
#       }
#       
#     }
#   }
# }
# 
# if(TEST){
#   save(errors_Y, errors_X, errors_B, file = paste(tests_dir, "errors_test.RData", sep = ''))
# } else{
#   save(errors_Y, errors_X, errors_B, file = paste(tests_dir, "errors_train.RData", sep = ''))
# }




