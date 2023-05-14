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

# Dimensions
L <- 1;
K <- 60*60;
S1 <- K;
S2 <- 9*9;
N <- 50;

# Grid
x <- seq(0, 1, length.out = sqrt(K))
y <- seq(0, 1, length.out = sqrt(K))

# Directories
test_dirs <- c()
for(i in 1:5){
  test_dirs[i] <- paste(tests_dir, "2D_test", i-1, "/", sep = '')
  if (!file.exists(test_dirs[i])){
    dir.create(test_dirs[i])
  }
}

for(i in 1:6){

  data <- generate_2d_data(x, y, S2, N, i, 0.95)
  
  # Response
  Y <- data[["Y"]]
  Y_clean <- data[["Y_clean"]]
  
  # Coefficient function
  B <- data[["B"]]

  # At nodes
  X_nodes <- data[["X_nodes"]]
  X_clean_nodes <- data[["X_clean_nodes"]]
  X_clean <- data[["X_clean"]]
  
  # At locations (#locations == #nodes)
  locations <- data[["locations"]]
  X_clean_locations <- data[["X_clean_locations"]]
  X_locations <- data[["X_locations"]]
  
  # At locations (#locations < #nodes)
  picked <- sample.int(dim(locations)[1], dim(locations)[1]*0.1)
  locations_less <- locations[picked,]
  X_clean_locations_less <- X_clean_locations[, picked]
  X_locations_less <- X_locations[, picked]
  
  # At locations (#locations << #nodes, equispaced subgrid)
  locations_sub <- data[["locations_sub"]]
  X_clean_locations_sub <- data[["X_clean_locations_sub"]]
  X_locations_sub <- data[["X_locations_sub"]]
  
  
  ## |||||||||||
  ## Test 0 ----
  ## |||||||||||
  
  # Comparison with Harold's results
  
  test_dir <- test_dirs[1]
  sub_test_dir <- paste(test_dir, "test", i, "/", sep = '')
  if (!file.exists(sub_test_dir)){
    dir.create(sub_test_dir)
  }
  
  # r1fpls_fem solution
  FEM_basis <- data[["basisobj"]]
  mesh <- data[["mesh"]]
  results_fpls_fem <- r1fpls_fem(X_nodes, Y, ncomp = 3, center = TRUE,
                                 basisobj = FEM_basis, penalty = 10,
                                 verbose = TRUE )
  
  # FPLS quanities of interest
  W <- as.matrix(results_fpls_fem[["W"]])
  V <- as.matrix(results_fpls_fem[["V"]])
  TT <- as.matrix(results_fpls_fem[["TT"]])
  C <- as.matrix(results_fpls_fem[["C"]])
  D <- as.matrix(results_fpls_fem[["D"]])
  write.csv(V, paste(sub_test_dir, "/V.csv", sep = ''))
  write.csv(W, paste(sub_test_dir, "/W.csv", sep = ''))
  write.csv(TT, paste(sub_test_dir, "/T.csv", sep = ''))
  write.csv(C, paste(sub_test_dir, "/C.csv", sep = ''))
  write.csv(D, paste(sub_test_dir, "/D.csv", sep = ''))

  # X field
  X_mean <- as.matrix(results_fpls_fem[["X_mean"]], ncol = K, nrow = 1, byrow = TRUE)
  X_hat <- TT %*% t(C) + matrix(X_mean, ncol = K, nrow = N, byrow = TRUE)
  write.csv(X_clean_nodes, paste(sub_test_dir, "/X_clean.csv", sep = ''))
  write.csv(X_nodes, paste(sub_test_dir, "/X.csv", sep = ''))
  write.csv(X_hat, paste(sub_test_dir, "/X_hat.csv", sep = ''))
  write.csv(X_mean, paste(sub_test_dir, "/X_mean.csv", sep = ''))
  

  # Response
  Y_mean <- as.matrix(results_fpls_fem[["Y_mean"]], ncol = L, nrow = 1, byrow = TRUE)
  Y_hat <- as.matrix(results_fpls_fem[["fitted.values"]])
  write.csv(Y_clean, paste(sub_test_dir, "/Y_clean.csv", sep = ''))
  write.csv(Y, paste(sub_test_dir, "/Y.csv", sep = ''))
  write.csv(Y_hat, paste(sub_test_dir, "/Y_hat.csv", sep = ''))
  write.csv(Y_mean, paste(sub_test_dir, "/Y_mean.csv", sep = ''))
  
  # Coefficient function
  B_hat <- results_fpls_fem[["coefficient_function"]]
  write.csv(B, paste(sub_test_dir, "/B.csv", sep = ''))
  write.csv(B_hat, paste(sub_test_dir, "/B_hat.csv", sep = ''))
  
  
  ## |||||||||||
  ## Test 1 ----
  ## |||||||||||
  
  test_dir <- test_dirs[2]
  sub_test_dir <- paste(test_dir, "test", i, "/", sep = '')
  if (!file.exists(sub_test_dir)){
    dir.create(sub_test_dir)
  }
  
  # X field
  write.csv(X_clean_nodes, paste(sub_test_dir, "/X_clean.csv", sep = ''))
  write.csv(X_nodes, paste(sub_test_dir, "/X.csv", sep = ''))
  
  # Response
  write.csv(Y_clean, paste(sub_test_dir, "/Y_clean.csv", sep = ''))
  write.csv(Y, paste(sub_test_dir, "/Y.csv", sep = ''))

  # Coefficient function
  write.csv(B, paste(sub_test_dir, "/B.csv", sep = ''))
  
  
  ## |||||||||||
  ## Test 2 ----
  ## |||||||||||
  
  test_dir <- test_dirs[3]
  sub_test_dir <- paste(test_dir, "test", i, "/", sep = '')
  if (!file.exists(sub_test_dir)){
    dir.create(sub_test_dir)
  }
  
  # Locations
  write.csv(locations, paste(sub_test_dir, "/locations.csv", sep = ''))
  
  # X field
  write.csv(X_clean_nodes, paste(sub_test_dir, "/X_clean.csv", sep = ''))
  write.csv(X_locations, paste(sub_test_dir, "/X_locations.csv", sep = ''))
  
  # Response
  write.csv(Y_clean, paste(sub_test_dir, "/Y_clean.csv", sep = ''))
  write.csv(Y, paste(sub_test_dir, "/Y.csv", sep = ''))
  
  # Coefficient function
  write.csv(B, paste(sub_test_dir, "/B.csv", sep = ''))
  
  
  ## |||||||||||
  ## Test 3 ----
  ## |||||||||||
  
  test_dir <- test_dirs[4]
  sub_test_dir <- paste(test_dir, "test", i, "/", sep = '')
  if (!file.exists(sub_test_dir)){
    dir.create(sub_test_dir)
  }
  
  # Locations
  write.csv(locations_less, paste(sub_test_dir, "/locations.csv", sep = ''))
  
  # X field
  write.csv(X_clean_nodes, paste(sub_test_dir, "/X_clean.csv", sep = ''))
  write.csv(X_locations_less, paste(sub_test_dir, "/X_locations.csv", sep = ''))
  
  # Response
  write.csv(Y_clean, paste(sub_test_dir, "/Y_clean.csv", sep = ''))
  write.csv(Y, paste(sub_test_dir, "/Y.csv", sep = ''))
  
  # Coefficient function
  write.csv(B, paste(sub_test_dir, "/B.csv", sep = ''))
  
  
  ## |||||||||||
  ## Test 4 ----
  ## |||||||||||
  
  test_dir <- test_dirs[5]
  sub_test_dir <- paste(test_dir, "test", i, "/", sep = '')
  if (!file.exists(sub_test_dir)){
    dir.create(sub_test_dir)
  }
  
  # Locations
  write.csv(locations_sub, paste(sub_test_dir, "/locations.csv", sep = ''))
  
  # X field
  write.csv(X_clean_nodes, paste(sub_test_dir, "/X_clean.csv", sep = ''))
  write.csv(X_locations_sub, paste(sub_test_dir, "/X_locations.csv", sep = ''))
  
  # Response
  write.csv(Y_clean, paste(sub_test_dir, "/Y_clean.csv", sep = ''))
  write.csv(Y, paste(sub_test_dir, "/Y.csv", sep = ''))
  
  # Coefficient function
  write.csv(B, paste(sub_test_dir, "/B.csv", sep = ''))


}
