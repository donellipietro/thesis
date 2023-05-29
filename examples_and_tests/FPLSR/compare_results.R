# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Utility to compare quantitatively the results of fPLSR with Harold %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()


# ||||||||||||||||||||||||||
# Environment variables ----
# ||||||||||||||||||||||||||

# Where to get the data
# tests_dir = "../../fdaPDE/test/data/models/FPLSR/2D_test_comparison/"
tests_dir <- "2D_test_comparison/"
n_batches <- 20
n_tests <- 6

TEST <- TRUE
# TRUE: test error is computed
# FALSE: training erro is computed

# ||||||||||||
# Results ----
# ||||||||||||

path_images <- "images/"
if (!file.exists(path_images)){
  dir.create(path_images)
}

# Legend:
# - 0: Harold R
# - hcoo: Harold c++
# - ns: No smoothing at all (Only diff with H. is R0)
# - sr: Smoothing only in regression
# - sri: Smoothing in the estimation of the mean and in regression
# - PLS: Multivataite PLS (NIPALS algorithm)
# - SIMPLS: Multivataite PLS (SIMPLS algorithm)

test_name_option_vect <- c("hcpp", "ns", "sr", "sri", "PLS", "SIMPLS")
n_test_options <- length(test_name_option_vect)

errors_Y <- matrix(0, nrow = n_batches, ncol = (n_test_options+1)*n_tests)
errors_X  <- matrix(0, nrow = n_batches, ncol = (n_test_options+1)*n_tests)
errors_B  <- matrix(0, nrow = n_batches, ncol = (n_test_options+1)*n_tests)
errors_B_min  <- matrix(0, nrow = n_batches, ncol = (n_test_options+1)*n_tests)
errors_B_max  <- matrix(0, nrow = n_batches, ncol = (n_test_options+1)*n_tests)



for(i in 1:n_tests){
    
  path <- paste(tests_dir, "test", i, "/", sep = '')
    
  B <- read.csv(paste(path,"B.csv", sep = ''), header = TRUE)[,2]
  n_nodes = length(B)
    
  for(j in 1:n_batches){
      
    if(TEST)
      Y_clean <- read.csv(paste(path,"Y_clean_", j, ".csv", sep = ''), header = TRUE)[,2]
    else
      Y_clean <- read.csv(paste(path,"Y_", j, ".csv", sep = ''), header = TRUE)[,2]
    
    batch_size = length(Y_clean)
    Y_hat_harold <- read.csv(paste(path,"Y_hat_", j, ".csv", sep = ''), header = TRUE)[,2]
    errors_Y[j, (n_test_options+1)*i - n_test_options] <- sum((Y_clean-Y_hat_harold)^2)/batch_size
    
    if(TEST)
      X_clean <- read.csv(paste(path,"X_clean_", j, ".csv", sep = ''), header = TRUE)[,-1]
    else
      X_clean <- read.csv(paste(path,"X_", j, ".csv", sep = ''), header = TRUE)[,-1]
    
    X_hat_harold <- read.csv(paste(path,"X_hat_", j, ".csv", sep = ''), header = TRUE)[,-1]
    errors_X[j, (n_test_options+1)*i - n_test_options] <- sum((X_clean-X_hat_harold)^2)/(n_nodes*batch_size)
    
    if(TEST){
      B_hat_harold <- read.csv(paste(path,"B_hat_", j, ".csv", sep = ''), header = TRUE)[,2]
      errors_B[j, (n_test_options+1)*i - n_test_options] <- sum((B-B_hat_harold)^2)/n_nodes
      # errors_B_min[j, (n_test_options+1)*i - n_test_options] <- min(abs(B-B_hat_harold))
      # errors_B_max[j, (n_test_options+1)*i - n_test_options] <- max(abs(B-B_hat_harold))
    }
      
    for(t in 1:n_test_options){

      test_name_option <- test_name_option_vect[t]

      # Y
      Y_hat_test <- read.csv(paste(path,"results/Y_hat_", test_name_option, "_", j, ".csv", sep = ''), header = FALSE)$V1
      errors_Y[j, (n_test_options+1)*i - n_test_options + t] <- sum((Y_clean-Y_hat_test)^2)/batch_size
      
      # X
      X_hat_test <- read.csv(paste(path,"results/X_hat_", test_name_option, "_", j, ".csv", sep = ''), header = FALSE)
      errors_X[j, (n_test_options+1)*i - n_test_options + t] <- sum((X_clean-X_hat_test)^2)/(n_nodes*batch_size)

      # B
      if(TEST){
        B_hat_test <- read.csv(paste(path,"results/B_hat_", test_name_option, "_", j, ".csv", sep = ''), header = FALSE)[,1]
        errors_B[j, (n_test_options+1)*i - n_test_options + t] <- sum((B-B_hat_test)^2)/n_nodes
        # errors_B_min[j, (n_test_options+1)*i- n_test_options + t] <- min(abs(B-B_hat_test))
        # errors_B_max[j, (n_test_options+1)*i- n_test_options + t] <- max(abs(B-B_hat_test))
      }
      

    }
    
  }
  
}

if(TEST){
  type = ""
  save(errors_Y, errors_X, errors_B, file = paste(tests_dir, "errors_test.RData", sep = ''))
  load(paste(tests_dir, "errors_test.RData", sep = ''))
} else{
  type = "train"
  save(errors_Y, errors_X, errors_B, file = paste(tests_dir, "errors_train.RData", sep = ''))
  load(paste(tests_dir, "errors_train.RData", sep = ''))
}
  
tests_to_be_shown <- 1:7
cols_to_be_shown <- rep(1:7 %in% tests_to_be_shown, n_tests)

colors <- c("red", "orange", "purple", "blue", "lightblue", "darkgreen", "green")
names <- c("Harld R", "Harold C++",
           "fPLSR no smoothing", "fPLSR smoothing regression", "fPLSR smoothing reg. + init.", 
           "M-PLSR (NIPALS)", "M_PLSR (SIMPLS)")

#jpeg(file=paste(path_images, "comparison_Y.jpg", sep = ''))
par(mfrow = c(1,1))
boxplot(errors_Y[, cols_to_be_shown], xlab = "Tests", ylab = "MSE",
        col = colors[tests_to_be_shown],
        main = paste("MSE", type, "- Y"), xaxt = 'n')
axis(1, at=(1:n_tests)*(length(tests_to_be_shown))-(length(tests_to_be_shown)-1)/2, labels=1:6)
abline(v = (0:n_tests)*(length(tests_to_be_shown))+1/2, lty = 2, col = "grey")
grid()
#dev.off()

#jpeg(file=paste(path_images, "comparison_X.jpg", sep = ''))
boxplot(errors_X[,cols_to_be_shown], xlab = "Tests", ylab = "MSE",
        col = colors[tests_to_be_shown],
        main = paste("MSE", type, "- X"), xaxt = 'n')
axis(1, at=(1:n_tests)*(length(tests_to_be_shown))-(length(tests_to_be_shown)-1)/2, labels=1:6)
abline(v = (0:n_tests)*(length(tests_to_be_shown))+1/2, lty = 2, col = "grey")
grid()
#dev.off()

#jpeg(file=paste(path_images, "comparison_X_zoom.jpg", sep = ''))
boxplot(errors_X[,cols_to_be_shown], xlab = "Tests", ylab = "MSE",
        col = colors[tests_to_be_shown], 
        main = "MSE - X (zoom)", xaxt = 'n',
        ylim = c(0.003, 0.004))
axis(1, at=(1:n_tests)*(length(tests_to_be_shown))-(length(tests_to_be_shown)-1)/2, labels=1:6)
abline(v = (0:n_tests)*(length(tests_to_be_shown))+1/2, lty = 2, col = "grey")
grid()
#dev.off()

if(TEST){
  #jpeg(file=paste(path_images, "comparison_B.jpg", sep = ''))
  par(mfrow = c(1,1))
  boxplot((errors_B[,cols_to_be_shown]), xlab = "Tests", ylab = "MSE",
          col = colors[tests_to_be_shown], 
          main = "MSE - B", xaxt = 'n')
  axis(1, at=(1:n_tests)*(length(tests_to_be_shown))-(length(tests_to_be_shown)-1)/2, labels=1:6)
  abline(v = (0:n_tests)*(length(tests_to_be_shown))+1/2, lty = 2, col = "grey")
  grid()
  #dev.off()
}

plot.new()
legend("center", names[tests_to_be_shown], col = colors[tests_to_be_shown], pch = 15)


# Test vs train
load(paste(tests_dir, "errors_test.RData", sep = ''))
errors_Y_test <- errors_Y
errors_X_test <- errors_X
errors_B_test <- errors_B
load(paste(tests_dir, "errors_train.RData", sep = ''))
errors_Y_train <- errors_Y
errors_X_train <- errors_X
errors_B_train <- errors_B


temp_test = matrix(colSums(errors_Y_test)/n_batches, ncol = n_test_options+1, byrow = TRUE)
temp_train = matrix(colSums(errors_Y_train)/n_batches, ncol = n_test_options+1, byrow = TRUE)

par(mfrow = c(1,2))
matplot(temp_train[,1:5], type = "l", lty = 2,
      main = "Functional methods", xlab = "Test", ylab = "mean MSE")
matplot(temp_test[,1:5], type = "l", lty = 1,  add = TRUE)
legend("topright", c("Test", "Train"), lty = c(1, 2))
grid()
matplot(temp_test[,6:7], type = "l", lty = 1,
        main = "Multivaraite methods", xlab = "Test", ylab = "mean MSE")
matplot(temp_train[,6:7], type = "l", lty = 2,  add = TRUE)
legend("topright", c("Test", "Train"), lty = c(1, 2))
grid()



