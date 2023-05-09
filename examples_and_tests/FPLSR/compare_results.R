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

# ||||||||||||||
# Functions ----
# ||||||||||||||

errors_Y <- matrix(0, nrow = n_batches, ncol = 2*n_tests)
errors_B  <- matrix(0, nrow = n_batches, ncol = 2*n_tests)

for(i in 1:n_tests){
  
  path <- paste(tests_dir, "test", i, "/", sep = '')
  
  B <- read.csv(paste(path,"B.csv", sep = ''), header = TRUE)[,2]
  
  for(j in 1:n_batches){
    
    # Y
    Y_clean <- read.csv(paste(path,"Y_clean_", j, ".csv", sep = ''), header = TRUE)[,2]
    Y_hat_harold <- read.csv(paste(path,"Y_hat_", j, ".csv", sep = ''), header = TRUE)[,2]
    Y_hat_fPLSR <- read.csv(paste(path,"results/Y_hat_", j, ".csv", sep = ''), header = FALSE)$V1
    batch_size = length(Y_clean)
    errors_Y[j, 2*i - 1] <- sum((Y_clean-Y_hat_harold)^2)/batch_size
    errors_Y[j, 2*i] <- sum((Y_clean-Y_hat_fPLSR)^2)/batch_size
    
    # B
    
    B_hat_harold <- read.csv(paste(path,"B_hat_", j, ".csv", sep = ''), header = TRUE)[,2]
    B_hat_fPLSR <- read.csv(paste(path,"results/B_hat_", j, ".csv", sep = ''), header = FALSE)$V1
    n_nodes = length(B)
    errors_B[j, 2*i - 1] <- sum((B-B_hat_harold)^2)/n_nodes
    errors_B[j, 2*i] <- sum((B-B_hat_fPLSR)^2)/n_nodes
    
  }
  
}

boxplot(errors_Y, xlab = "Tests", ylab = "MSE", col = c("red", "blue"), main = "MSE - Y", xaxt = 'n')
axis(1, at=(1:6)*2-1/2, labels=1:6)
legend("topright", c("Harold", "fPLSR"), col = c("red", "blue"), pch = c(15, 15))
grid()

boxplot(log(errors_B), xlab = "Tests", ylab = "log(MSE)", col = c("red", "blue"), main = "MSE - B")
axis(1, at=(1:6)*2-1/2, labels=1:6)
legend("topright", c("Harold", "fPLSR"), col = c("red", "blue"), pch = c(15, 15))
grid()

