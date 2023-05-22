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
n_batches <- 20
n_tests <- 6

# ||||||||||||
# Results ----
# ||||||||||||

path_images <- "images/"
if (!file.exists(path_images)){
  dir.create(path_images)
}

test_name_option_vect <- c("hcpp", "sr", "sri")
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
      
    Y_clean <- read.csv(paste(path,"Y_clean_", j, ".csv", sep = ''), header = TRUE)[,2]
    batch_size = length(Y_clean)
    Y_hat_harold <- read.csv(paste(path,"Y_hat_", j, ".csv", sep = ''), header = TRUE)[,2]
    errors_Y[j, (n_test_options+1)*i - n_test_options] <- sum((Y_clean-Y_hat_harold)^2)/batch_size
    
    X_clean <- read.csv(paste(path,"X_clean_", j, ".csv", sep = ''), header = TRUE)[,2]
    X_hat_harold <- read.csv(paste(path,"X_hat_", j, ".csv", sep = ''), header = TRUE)[,2]
    errors_X[j, (n_test_options+1)*i - n_test_options] <- sum((X_clean-X_hat_harold)^2)/n_nodes
    
    B_hat_harold <- read.csv(paste(path,"B_hat_", j, ".csv", sep = ''), header = TRUE)[,2]
    errors_B[j, (n_test_options+1)*i - n_test_options] <- sum((B-B_hat_harold)^2)/n_nodes
    errors_B_min[j, (n_test_options+1)*i - n_test_options] <- min(abs(B-B_hat_harold))
    errors_B_max[j, (n_test_options+1)*i - n_test_options] <- max(abs(B-B_hat_harold))
      
    for(t in 1:n_test_options){

      test_name_option <- test_name_option_vect[t]

      # Y
      Y_hat_test <- read.csv(paste(path,"results/Y_hat_", test_name_option, "_", j, ".csv", sep = ''), header = FALSE)$V1
      errors_Y[j, (n_test_options+1)*i - n_test_options + t] <- sum((Y_clean-Y_hat_test)^2)/batch_size
      
      # X
      X_hat_test <- read.csv(paste(path,"results/X_hat_", test_name_option, "_", j, ".csv", sep = ''), header = FALSE)$V1
      errors_X[j, (n_test_options+1)*i - n_test_options + t] <- sum((X_clean-X_hat_test)^2)/n_nodes

      # B
      B_hat_test <- read.csv(paste(path,"results/B_hat_", test_name_option, "_", j, ".csv", sep = ''), header = FALSE)$V1
      errors_B[j, (n_test_options+1)*i - n_test_options + t] <- sum((B-B_hat_test)^2)/n_nodes
      errors_B_min[j, (n_test_options+1)*i- n_test_options + t] <- min(abs(B-B_hat_test))
      errors_B_max[j, (n_test_options+1)*i- n_test_options + t] <- max(abs(B-B_hat_test))

    }
    
  }
  
}

jpeg(file=paste(path_images, "comparison_Y.jpg", sep = ''))
par(mfrow = c(1,1))
boxplot(errors_Y, xlab = "Tests", ylab = "MSE", col = c("red", "orange", "blue", "lightblue"), main = "MSE - Y", xaxt = 'n')
axis(1, at=(1:n_tests)*(n_test_options+1)-(n_test_options)/2, labels=1:6)
#legend("topright", c("Harold", "fPLSR"), col = c("red", "blue"), pch = c(15, 15))
grid()
dev.off()

jpeg(file=paste(path_images, "comparison_B.jpg", sep = ''))
par(mfrow = c(1,1))
boxplot((errors_B), xlab = "Tests", ylab = "MSE", col = c("red", "orange", "blue", "lightblue"), main = "MSE - B", xaxt = 'n')
axis(1, at=(1:n_tests)*(n_test_options+1)-(n_test_options)/2, labels=1:6)
# legend("topright", c("Harold", "fPLSR"), col = c("red", "blue"), pch = c(15, 15))
grid()
dev.off()

jpeg(file=paste(path_images, "comparison_X.jpg", sep = ''))
boxplot(errors_X, xlab = "Tests", ylab = "MSE", col = c("red", "orange", "blue", "lightblue"), main = "MSE - X", xaxt = 'n')
axis(1, at=(1:n_tests)*(n_test_options+1)-(n_test_options)/2, labels=1:6)
# legend("topright", c("Harold", "fPLSR"), col = c("red", "blue"), pch = c(15, 15))
grid()
dev.off()

# boxplot(errors_B_min, xlab = "Tests", ylab = "min error", col = c("red", "orange", "blue", "lightblue"), main = "min error - B", xaxt = 'n')
# axis(1, at=(1:n_tests)*(n_test_options+1)-(n_test_options)/2, labels=1:6)
# # legend("topright", c("Harold", "fPLSR"), col = c("red", "blue"), pch = c(15, 15))
# grid()
# 
# boxplot(errors_B_max, xlab = "Tests", ylab = "max error", col = c("red", "orange", "blue", "lightblue"), main = "max error - B", xaxt = 'n')
# axis(1, at=(1:n_tests)*(n_test_options+1)-(n_test_options)/2, labels=1:6)
# # legend("topright", c("Harold", "fPLSR"), col = c("red", "blue"), pch = c(15, 15))
# grid()

# dev.off()

