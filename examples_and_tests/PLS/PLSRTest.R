# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# %% Multivariate Principal Least Squares Regression %% #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

rm(list = ls())
graphics.off()


# ||||||||||||||
# Functions ----
# ||||||||||||||

load("PLSR.RData")
load("SIMPLS.RData")

# ||||||||||
# Tests ----
# ||||||||||

tests_dir <- "../../fdaPDE/test/data/models/FPLSR/2D_test_comparison/"
# tests_dir <- "../FPLSR/2D_test_comparison/"
n_samples <- 50
n_batches <- 10
n_tests <- 6

test_name_option_vect <- c("PLS", "PLS_noYDefl", "SIMPLS")
n_test_options <- length(test_name_option_vect)

errors_Y <- matrix(0, nrow = n_batches, ncol = n_test_options*n_tests)
errors_X  <- matrix(0, nrow = n_batches, ncol = n_test_options*n_tests)
errors_B  <- matrix(0, nrow = n_batches, ncol = n_test_options*n_tests)
errors_Y_train <- matrix(0, nrow = n_batches, ncol = n_test_options*n_tests)
errors_X_train <- matrix(0, nrow = n_batches, ncol = n_test_options*n_tests)


for(i in 1:n_tests){
  
  
  print("##########")
  print(paste("# Test", i, "#"))
  print("##########")
  
  path <- paste(tests_dir, "test", i, "/", sep = '')
  sub_test_dir <- paste(path, "results/", sep = '')
  
  B_clean <- read.csv(paste(path,"B.csv", sep = ''), header = TRUE)[,-1]
  n_nodes = length(B_clean)
  
  for(j in 1:n_batches){
    
    print(paste("- Batch #",  j, sep = ''))
    
    # Original data
    Y <- matrix(read.csv(paste(path,"Y_", j, ".csv", sep = ''), header = TRUE)$V1, ncol = 1, nrow = n_samples)
    X <- read.csv(paste(path,"X_", j, ".csv", sep = ''), header = TRUE)[,-1]
    
    # Expected results
    Y_clean <- matrix(read.csv(paste(path,"Y_clean_", j, ".csv", sep = ''), header = TRUE)$V1, ncol = 1, nrow = n_samples)
    X_clean <- read.csv(paste(path,"X_clean_", j, ".csv", sep = ''), header = TRUE)[,-1]
    batch_size = length(Y_clean)
    
    plsr = PLSR(X, Y, 3, FALSE)
    # write.table(plsr$Y_hat, paste(sub_test_dir, "/Y_hat_", test_name_option_vect[1], "_",j,".csv", sep = ''), col.names = FALSE, row.names = FALSE, sep = ",")
    # write.table(plsr$X_hat, paste(sub_test_dir, "/X_hat_", test_name_option_vect[1], "_",j,".csv", sep = ''), col.names = FALSE, row.names = FALSE, sep = ",")
    # write.table(plsr$Beta, paste(sub_test_dir, "/B_hat_", test_name_option_vect[1], "_",j,".csv", sep = ''), col.names = FALSE, row.names = FALSE, sep = ",")
    errors_Y[j, n_test_options * (i - 1) + 1] <- sum((Y_clean-plsr$Y_hat)^2)/(batch_size)
    errors_X[j, n_test_options * (i - 1) + 1] <- sum((X_clean-plsr$X_hat)^2)/(n_nodes*batch_size)
    errors_B[j, n_test_options * (i - 1) + 1] <- sum((B_clean-plsr$Beta)^2)/(n_nodes)
    errors_Y_train[j, n_test_options * (i - 1) + 1] <- sum((Y-plsr$Y_hat)^2)/(batch_size)
    errors_X_train[j, n_test_options * (i - 1) + 1] <- sum((X-plsr$X_hat)^2)/(n_nodes*batch_size)
    
    plsr_noYDefl = PLSR(X, Y, 3, TRUE)
    # write.table(plsr$Y_hat, paste(sub_test_dir, "/Y_hat_", test_name_option_vect[2], "_",j,".csv", sep = ''), col.names = FALSE, row.names = FALSE, sep = ",")
    # write.table(plsr$X_hat, paste(sub_test_dir, "/X_hat_", test_name_option_vect[2], "_",j,".csv", sep = ''), col.names = FALSE, row.names = FALSE, sep = ",")
    # write.table(plsr$Beta, paste(sub_test_dir, "/B_hat_", test_name_option_vect[2], "_",j,".csv", sep = ''), col.names = FALSE, row.names = FALSE, sep = ",")
    errors_Y[j, n_test_options * (i - 1) + 2] <- sum((Y_clean-plsr_noYDefl$Y_hat)^2)/(batch_size)
    errors_X[j, n_test_options * (i - 1) + 2] <- sum((X_clean-plsr_noYDefl$X_hat)^2)/(n_nodes*batch_size)
    errors_B[j, n_test_options * (i - 1) + 2] <- sum((B_clean-plsr_noYDefl$Beta)^2)/(n_nodes)
    errors_Y_train[j, n_test_options * (i - 1) + 2] <- sum((Y-plsr_noYDefl$Y_hat)^2)/(batch_size)
    errors_X_train[j, n_test_options * (i - 1) + 2] <- sum((X-plsr_noYDefl$X_hat)^2)/(n_nodes*batch_size)
    
    simpls = SIMPLS(X, Y, 3)
    # write.table(plsr$Y_hat, paste(sub_test_dir, "/Y_hat_", test_name_option_vect[3], "_",j,".csv", sep = ''), col.names = FALSE, row.names = FALSE, sep = ",")
    # write.table(plsr$X_hat, paste(sub_test_dir, "/X_hat_", test_name_option_vect[3], "_",j,".csv", sep = ''), col.names = FALSE, row.names = FALSE, sep = ",")
    # write.table(plsr$Beta, paste(sub_test_dir, "/B_hat_", test_name_option_vect[3], "_",j,".csv", sep = ''), col.names = FALSE, row.names = FALSE, sep = ",")
    errors_Y[j, n_test_options * (i - 1) + 3] <- sum((Y_clean-simpls$Y_hat)^2)/(batch_size)
    errors_X[j, n_test_options * (i - 1) + 3] <- sum((X_clean-simpls$X_hat)^2)/(n_nodes*batch_size)
    errors_B[j, n_test_options * (i - 1) + 3] <- sum((B_clean-simpls$Beta)^2)/(n_nodes)
    errors_Y_train[j, n_test_options * (i - 1) + 3] <- sum((Y-simpls$Y_hat)^2)/(batch_size)
    errors_X_train[j, n_test_options * (i - 1) + 3] <- sum((X-simpls$X_hat)^2)/(n_nodes*batch_size)
    
    # results[["PLS"]] <- list(Y_hat = plsr$Y_hat,
    #                          X_hat = plsr$X_hat,
    #                          B_hat = plsr$Beta)
    # results[["PLS_noYDefl"]] <- list(Y_hat = plsr_noYDefl$Y_hat,
    #                                  X_hat = plsr_noYDefl$X_hat,
    #                                  B_hat = plsr_noYDefl$Beta)
    # results[["SIMPLS"]] <- list(Y_hat = simpls$Y_hat,
    #                             X_hat = simpls$X_hat,
    #                             B_hat = simpls$Beta)
    
  }
  
}

write.csv(errors_Y, paste(tests_dir, "errors_Y_multivariate.csv", sep = ''))
write.csv(errors_X, paste(tests_dir, "errors_X_multivariate.csv", sep = ''))
write.csv(errors_B, paste(tests_dir, "errors_B_multivariate.csv", sep = ''))
write.csv(errors_Y_train, paste(tests_dir, "errors_Y_train_multivariate.csv", sep = ''))
write.csv(errors_X_train, paste(tests_dir, "errors_X_train_multivariate.csv", sep = ''))

tests_dir <- "../../fdaPDE/test/data/models/FPLSR_SIMPLS/2D_test_comparison/"

write.csv(errors_Y, paste(tests_dir, "errors_Y_multivariate.csv", sep = ''))
write.csv(errors_X, paste(tests_dir, "errors_X_multivariate.csv", sep = ''))
write.csv(errors_B, paste(tests_dir, "errors_B_multivariate.csv", sep = ''))
write.csv(errors_Y_train, paste(tests_dir, "errors_Y_train_multivariate.csv", sep = ''))
write.csv(errors_X_train, paste(tests_dir, "errors_X_train_multivariate.csv", sep = ''))


