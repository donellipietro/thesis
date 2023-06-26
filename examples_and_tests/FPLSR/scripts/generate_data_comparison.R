# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Utility to generate data for FPLSR testing %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()


# |||||||||||||||||||||||||||||||||||
# Import Libraries and functions ----
# |||||||||||||||||||||||||||||||||||

library(pls)
library(fda)
library(pR1FPLS)

load("generate_2d_data.RData")


# ||||||||||||||||||||||||||
# Environment variables ----
# ||||||||||||||||||||||||||

complete = FALSE
tests_dir <- "tests/"


# ||||||||||||||||||
# Generate data ----
# ||||||||||||||||||

if (!file.exists(tests_dir)){
  dir.create(tests_dir)
}

# Grid
x <- seq(0, 1, length.out = 60)
y <- seq(0, 1, length.out = 60)


## |||||||||||||
## Test 0.1 ----
## |||||||||||||

# Comparison with Harold's results

# Directory
test_dir <- paste(tests_dir, "2D_test_comparison/", sep = '')
if (!file.exists(test_dir)){
  dir.create(test_dir)
}

print(" ")
print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print("% FPLSR comparison data generation %")
print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print(" ")



for(i in 1:6){
  
  print("###########")
  print(paste("Test ", i, sep = ""))
  print("###########")
  
  n_batches <- 20
  batch_size <- 50
  n_samples <- n_batches*batch_size
  
  data <- generate_2d_data(x, y, 81, n_samples, i, 0.95)
  
  # Response
  Y <- data[["Y"]]
  Y_clean <- data[["Y_clean"]]
  
  # Coefficient function
  B <- data[["B"]]
  
  # At nodes
  X_nodes <- data[["X_nodes"]]
  X_clean_nodes <- data[["X_clean_nodes"]]
  X_clean <- data[["X_clean"]]
  
  # Create directory
  sub_test_dir <- paste(test_dir, "test", i, "/", sep = '')
  if (!file.exists(sub_test_dir)){
    dir.create(sub_test_dir)
  }
  
  # Export B
  write.csv(B, paste(sub_test_dir, "/B.csv", sep = ''))
  
  for(j in 1:n_batches){
    
    # Batch data
    picked <- ((j-1)*batch_size+1):(j*batch_size)
    X_clean_batch <- X_clean_nodes[picked,]
    X_batch <- X_nodes[picked,]
    Y_batch <- as.matrix(Y[picked,], nrow = batch_size, ncol = 1)
    Y_clean_batch <- as.matrix(Y_clean[picked,], nrow = batch_size, ncol = 1)
  
    # r1fpls_fem solution
    FEM_basis <- data[["basisobj"]]
    mesh <- data[["mesh"]]
    results_fpls_fem <- r1fpls_fem(X_batch, Y_batch, ncomp = 3, center = TRUE,
                                   basisobj = FEM_basis, penalty = 10,
                                   verbose = TRUE )
    
    # X field
    write.csv(X_clean_batch, paste(sub_test_dir, "/X_clean_",j,".csv", sep = ''))
    write.csv(X_batch, paste(sub_test_dir, "/X_",j,".csv", sep = ''))
    
    # Response
    Y_hat_batch <- as.matrix(results_fpls_fem[["fitted.values"]])
    write.csv(Y_clean_batch, paste(sub_test_dir, "/Y_clean_",j,".csv", sep = ''))
    write.csv(Y_batch, paste(sub_test_dir, "/Y_",j,".csv", sep = ''))
    write.csv(Y_hat_batch, paste(sub_test_dir, "/Y_hat_",j,".csv", sep = ''))
    
    # Field
    TT <- as.matrix(results_fpls_fem[["TT"]])
    C <- as.matrix(results_fpls_fem[["C"]])
    K <- dim(C)[1]
    X_mean <- matrix(as.matrix(results_fpls_fem[["X_mean"]]), ncol = K, nrow = batch_size, byrow = TRUE)
    X_hat_batch <- TT %*% t(C) + X_mean
    write.csv(X_hat_batch, paste(sub_test_dir, "/X_hat_",j,".csv", sep = ''))
    
    # Coefficient function
    B_hat_batch <- results_fpls_fem[["coefficient_function"]]
    write.csv(B_hat_batch, paste(sub_test_dir, "/B_hat_",j,".csv", sep = ''))
    
    # FPLS quanities of interest
    # W <- as.matrix(results_fpls_fem[["W"]])
    
    # D <- as.matrix(results_fpls_fem[["D"]])
    # write.csv(W, paste(sub_test_dir, "/W.csv", sep = ''))
    # write.csv(TT, paste(sub_test_dir, "/T.csv", sep = ''))
    # write.csv(C, paste(sub_test_dir, "/C.csv", sep = ''))
    # write.csv(D, paste(sub_test_dir, "/D.csv", sep = ''))
  
  }
  
}
