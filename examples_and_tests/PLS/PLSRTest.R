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

# tests_dir <- "../../fdaPDE/test/data/models/FPLSR/2D_test_comparison/"
tests_dir <- "../FPLSR/2D_test_comparison/"
n_samples <- 50
n_batches <- 20
n_tests <- 6

test_name_option_vect <- c("PLS", "PLS_noYDefl", "SIMPLS")
n_test_options <- length(test_name_option_vect)


for(i in 1:n_tests){
  
  path <- paste(tests_dir, "test", i, "/", sep = '')
  sub_test_dir <- paste(path, "results/", sep = '')
  
  for(j in 1:n_batches){
    
    # Original data
    Y <- matrix(read.csv(paste(path,"Y_", j, ".csv", sep = ''), header = TRUE)$V1, ncol = 1, nrow = n_samples)
    X <- read.csv(paste(path,"X_", j, ".csv", sep = ''), header = TRUE)[,-1]
    
    plsr = PLSR(X, Y, 3, FALSE)
    write.table(plsr$Y_hat, paste(sub_test_dir, "/Y_hat_", test_name_option_vect[1], "_",j,".csv", sep = ''), col.names = FALSE, row.names = FALSE, sep = ",")
    write.table(plsr$X_hat, paste(sub_test_dir, "/X_hat_", test_name_option_vect[1], "_",j,".csv", sep = ''), col.names = FALSE, row.names = FALSE, sep = ",")
    write.table(plsr$Beta, paste(sub_test_dir, "/B_hat_", test_name_option_vect[1], "_",j,".csv", sep = ''), col.names = FALSE, row.names = FALSE, sep = ",")
    
    plsr_noYDefl = PLSR(X, Y, 3, TRUE)
    write.table(plsr$Y_hat, paste(sub_test_dir, "/Y_hat_", test_name_option_vect[2], "_",j,".csv", sep = ''), col.names = FALSE, row.names = FALSE, sep = ",")
    write.table(plsr$X_hat, paste(sub_test_dir, "/X_hat_", test_name_option_vect[2], "_",j,".csv", sep = ''), col.names = FALSE, row.names = FALSE, sep = ",")
    write.table(plsr$Beta, paste(sub_test_dir, "/B_hat_", test_name_option_vect[2], "_",j,".csv", sep = ''), col.names = FALSE, row.names = FALSE, sep = ",")
    
    simpls = SIMPLS(X, Y, 3)
    write.table(plsr$Y_hat, paste(sub_test_dir, "/Y_hat_", test_name_option_vect[3], "_",j,".csv", sep = ''), col.names = FALSE, row.names = FALSE, sep = ",")
    write.table(plsr$X_hat, paste(sub_test_dir, "/X_hat_", test_name_option_vect[3], "_",j,".csv", sep = ''), col.names = FALSE, row.names = FALSE, sep = ",")
    write.table(plsr$Beta, paste(sub_test_dir, "/B_hat_", test_name_option_vect[3], "_",j,".csv", sep = ''), col.names = FALSE, row.names = FALSE, sep = ",")
    
    
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