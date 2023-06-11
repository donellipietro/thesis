# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Utility to compare quantitatively the results of fPLSR with Harold %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()


# ||||||||||||||||||||||||||
# Environment variables ----
# ||||||||||||||||||||||||||

# Where to get the data
tests_dir = "../../fdaPDE/test/data/models/FPLSR/2D_test_comparison/"
# tests_dir <- "2D_test_comparison/"
n_batches <- 20
n_tests <- 6

test_name_option_vect <- c("hcpp_l0", "ns_l0", "hcpp", "ns", "sr", "sri", "hcpp-KCV", "sri-GCV")
n_test_options <- length(test_name_option_vect)

TEST <- FALSE
# TRUE: test error is computed
# FALSE: training error is computed

path_images <- "images/" # _SIMPLS
if (!file.exists(path_images)){
  dir.create(path_images)
}
path_images <- "images/comparison/" # _SIMPLS
if (!file.exists(path_images)){
  dir.create(path_images)
}

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

plot_results_divided <- function(n_test_options, tests_to_be_shown, n_tests, errors_Y, errors_X, errors_B, TEST, type, colors, names) {
  
  cols_to_be_shown <- matrix(rep(tests_to_be_shown, n_tests), ncol = length(tests_to_be_shown), byrow = TRUE)
  cols_to_be_shown <- cols_to_be_shown + 0:(n_tests-1)*n_test_options
  
  if(type == "train")
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


# ||||||||||||
# Results ----
# ||||||||||||

# Legend:
# - 0: Harold R
# - hcpp_l0: Harold c++, lambda = 0 
# - ns_lo: No smoothing at all (Only diff with H. is R0), lambda = 0 
# - hcpp: Harold c++
# - ns: No smoothing at all (Only diff with H. is R0)
# - sr: Smoothing only in regression
# - sri: Smoothing in the estimation of the mean and in regression
# - PLS: Multivataite PLS (NIPALS algorithm)
# - SIMPLS: Multivataite PLS (SIMPLS algorithm)

# errors_Y <- matrix(0, nrow = n_batches, ncol = (n_test_options+1)*n_tests)
# errors_X  <- matrix(0, nrow = n_batches, ncol = (n_test_options+1)*n_tests)
# errors_B  <- matrix(0, nrow = n_batches, ncol = (n_test_options+1)*n_tests)
# errors_B_min  <- matrix(0, nrow = n_batches, ncol = (n_test_options+1)*n_tests)
# errors_B_max  <- matrix(0, nrow = n_batches, ncol = (n_test_options+1)*n_tests)

# ||||||||||||||||||||||
# Visualize results ----
# ||||||||||||||||||||||

# load errors
if(TEST){
  type = ""
  errors_Y <- read.csv(paste(tests_dir, "errors_Y.csv", sep = ''), header = FALSE)
  errors_X <- read.csv(paste(tests_dir, "errors_X.csv", sep = ''), header = FALSE)
  errors_B <- read.csv(paste(tests_dir, "errors_B.csv", sep = ''), header = FALSE)
  errors_Y_multivariate <- read.csv(paste(tests_dir, "errors_Y_multivariate.csv", sep = ''), header = TRUE)[,-1]
  errors_X_multivariate <- read.csv(paste(tests_dir, "errors_X_multivariate.csv", sep = ''), header = TRUE)[,-1]
  errors_B_multivariate <- read.csv(paste(tests_dir, "errors_B_multivariate.csv", sep = ''), header = TRUE)[,-1]
} else{
  type = "train"
  errors_Y <- read.csv(paste(tests_dir, "errors_Y_train.csv", sep = ''), header = FALSE)
  errors_X <- read.csv(paste(tests_dir, "errors_X_train.csv", sep = ''), header = FALSE)
  errors_Y_multivariate <- read.csv(paste(tests_dir, "errors_Y_train_multivariate.csv", sep = ''), header = TRUE)[,-1]
  errors_X_multivariate <- read.csv(paste(tests_dir, "errors_X_train_multivariate.csv", sep = ''), header = TRUE)[,-1]
}

errors_Y_merged = merge_errors(errors_Y, errors_Y_multivariate, n_tests, n_test_options, 3)
errors_X_merged = merge_errors(errors_X, errors_X_multivariate, n_tests, n_test_options, 3)
if(TEST)
  errors_B_merged = merge_errors(errors_B, errors_B_multivariate, n_tests, n_test_options, 3)

colors <- c("red", "orange",
            "darkgrey",
            "purple", "blue", "lightblue",
            "darkred", "darkblue",
            "darkgreen", "green", "lightgreen") # "red"
names <- c("Harold C++, lambda = 0", "fPLSR no smoothing, lambda = 0",
           "Harold C++",
           "fPLSR no smoothing", "fPLSR smoothing regression", "fPLSR smoothing reg. + init.",
           "Harold C++, KCV", "fPLSR smoothing reg. + init., GCV",
           "M-PLSR (NIPALS)", "M-PLSR (NIPALS no Y deflation)", "M_PLSR (SIMPLS)") # "Harold R"

tests_to_be_shown <- c(3:6, 7:8, 1, 2, 9)
n_test_options_merged <- dim(errors_Y_merged)[2]/n_tests
plot_results_divided(n_test_options_merged, tests_to_be_shown, n_tests, errors_Y_merged, errors_X_merged, errors_B_merged, TEST, type, colors, names) 
# plot_results(n_test_options_merged, tests_to_be_shown, n_tests, errors_Y_merged, errors_X_merged, errors_B_merged, TEST, type, colors, names) 


# ||||||||||||||||||
# Test vs train ----
# ||||||||||||||||||

# Load errors
errors_Y <- read.csv(paste(tests_dir, "errors_Y.csv", sep = ''), header = FALSE)
errors_X <- read.csv(paste(tests_dir, "errors_X.csv", sep = ''), header = FALSE)
errors_Y_multivariate <- read.csv(paste(tests_dir, "errors_Y_multivariate.csv", sep = ''), header = TRUE)[,-1]
errors_X_multivariate <- read.csv(paste(tests_dir, "errors_X_multivariate.csv", sep = ''), header = TRUE)[,-1]
errors_Y_train <- read.csv(paste(tests_dir, "errors_Y_train.csv", sep = ''), header = FALSE)
errors_X_train <- read.csv(paste(tests_dir, "errors_X_train.csv", sep = ''), header = FALSE)
errors_Y_train_multivariate <- read.csv(paste(tests_dir, "errors_Y_train_multivariate.csv", sep = ''), header = TRUE)[,-1]
errors_X_train_multivariate <- read.csv(paste(tests_dir, "errors_X_train_multivariate.csv", sep = ''), header = TRUE)[,-1]

# merge dataframe
n_test_option_1 <- dim(errors_Y)[2]/n_tests
n_test_option_2 <- dim(errors_Y_multivariate)[2]/n_tests
errors_X_merged_test = merge_errors(errors_X, errors_X_multivariate, n_tests, n_test_option_1, n_test_option_2)
errors_Y_merged_test = merge_errors(errors_Y, errors_Y_multivariate, n_tests, n_test_option_1, n_test_option_2)
errors_X_merged_train = merge_errors(errors_X_train, errors_X_train_multivariate, n_tests, n_test_option_1, n_test_option_2)
errors_Y_merged_train = merge_errors(errors_Y_train, errors_Y_train_multivariate, n_tests, n_test_option_1, n_test_option_2)

temp_Y_test = matrix(colSums(errors_Y_merged_test)/n_batches, ncol = n_test_options_merged, byrow = TRUE)
temp_Y_train = matrix(colSums(errors_Y_merged_train)/n_batches, ncol = n_test_options_merged, byrow = TRUE)
temp_X_test = matrix(colSums(errors_X_merged_test)/n_batches, ncol = n_test_options_merged, byrow = TRUE)
temp_X_train = matrix(colSums(errors_X_merged_train)/n_batches, ncol = n_test_options_merged, byrow = TRUE)


functional_methods <- c(3:6,7:8)
multivariate_methods <- c(1:2,9:11)


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




