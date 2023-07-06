# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Utility to generate data for FPLSR testing %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()


# |||||||||||||||||||||||||||||||||||
# Import Libraries and functions ----
# |||||||||||||||||||||||||||||||||||

library(pls)
library(fdaPDE)
library(pR1FPLS)

load("scripts/functions/generate_2D_data.RData")
load("scripts/functions/help_functions.RData")
load("../PLS/PLSR.RData")
load("../PLS/SIMPLS.RData")


# ||||||||||||||||||||||||||
# Environment variables ----
# ||||||||||||||||||||||||||

# Tests directory
tests_dir <- "tests/"

# Covariates
X.index_vect <- c(1)

# Responses
L <- 1;
B.indexes_list <- list(c(1), c(2), c(3), c(4), c(5))

# Number of samples
N <- 50;

# Locations
number.locations_per_element <- 1

# Mesh
# "c_shaped",  "unit_square", "unit_square_coarse", "unit_square_medium"
mesh.symbol <- c("s")
mesh.name_vec <- c("unit_square", "pippo")
mesh_finer.name_vec <- c("unit_square")
mesh.area_refine_vec <- c(0)
directory.test_functions_vec <- c("scripts/functions/tests_unit_square.RData") 

# Noise
STRATEGY = 1
# NSR.X <- (1/3)^2
NSR.Y <- (1/3)^2


# ||||||||||||||||||
# Generate data ----
# ||||||||||||||||||

# Directory
if (!file.exists(tests_dir)){
  dir.create(tests_dir)
}
test_dir <- paste(tests_dir, "2D_test/", sep = '')
if (!file.exists(test_dir)){
  dir.create(test_dir)
}

cat("\n")
cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
cat("% FPLSR tests data generation %\n")
cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
cat("\n")

for(m in 1:length(mesh.name_vec)){
  
  load(directory.test_functions_vec[m])
  mesh.name <- mesh.name_vec[m]
  mesh_finer.name <- mesh_finer.name_vec[m]
  mesh.area_refine <- mesh.area_refine_vec[m]
  
  for(i in 1:length(X.index_vect)){
    
    X.index <- X.index_vect[i]
    
    for(j in 1:length(B.indexes_list)){
      
      B.indexes <- B.indexes_list[[j]]
      
      set.seed(0)
      
      options_vec <- c()
      if(length(mesh.name_vec)>1)
        options_vec <- c(options_vec, mesh.symbol[m])
      if(length(X.index_vect)>1)
        options_vec <- c(options_vec, i)
      if(length(B.indexes_list)>1)
        options_vec <- c(options_vec, j)
      
      test.name <- paste("test", paste(options_vec, collapse="_"), sep = "_")
      cat("##########\n")
      cat(paste(test.name, "\n", sep = ''))
      cat("##########\n")
      
      sub_test_dir <- paste(test_dir, test.name, "/", sep = '')
      if (!file.exists(sub_test_dir)){
        dir.create(sub_test_dir)
      }
      
      
      # Generate data
      # |||||||||||||
      
      data <- generate_2D_data(mesh.name, mesh_finer.name, mesh.area_refine,    # mesh
                               number.locations_per_element,                    # locations
                               N,                                               # samples
                               generate_X, X.index, generate_B, B.indexes,      # tests
                               STRATEGY, NSR.Y,                                 # noise
                               FALSE, FALSE)                                    # plot
      
      # Response
      Y_clean <- data[["Y_clean"]]
      Y <- data[["Y"]]
      write.csv(Y_clean, paste(sub_test_dir, "Y_clean.csv", sep = ''))
      write.csv(Y, paste(sub_test_dir, "Y.csv", sep = ''))
      
      # Coefficient function
      B <- data[["B"]]
      write.csv(B, paste(sub_test_dir, "B.csv", sep = ''))
      
      # At nodes
      X_clean_nodes <- data[["X_clean_nodes"]]
      X_nodes <- data[["X_nodes"]]
      write.csv(X_clean_nodes, paste(sub_test_dir, "X_clean.csv", sep = ''))
      write.csv(X_nodes, paste(sub_test_dir, "X.csv", sep = ''))
      
      # At locations (1 node for each element)
      locations <- data[["locations"]]
      X_clean_locations <- data[["X_clean_locations"]]
      X_locations <- data[["X_locations"]]
      
      # Mesh
      FEM_basis <- data[["basisobj"]]
      mesh <- data[["mesh"]]
      
      # Dimensions
      K <- dim(mesh$nodes)[1]
      
      
      # Solution with r1fpls_fem
      # ||||||||||||||||||||||||
      
      # r1fpls_fem solution
      results_fpls_fem <- r1fpls_fem(X_nodes, Y, ncomp = 3, center = TRUE,
                                     basisobj = FEM_basis, penalty = 10,
                                     verbose = FALSE )
      
      method_dir <- paste(sub_test_dir, "r1fpls", "/", sep = '')
      if (!file.exists(method_dir)){
        dir.create(method_dir)
      }
      
      # FPLS quanities of interest
      W <- as.matrix(results_fpls_fem[["W"]])
      V <- as.matrix(results_fpls_fem[["V"]])
      TT <- as.matrix(results_fpls_fem[["TT"]])
      C <- as.matrix(results_fpls_fem[["C"]])
      D <- as.matrix(results_fpls_fem[["D"]])
      write.csv(V, paste(method_dir, "V.csv", sep = ''))
      write.csv(W, paste(method_dir, "W.csv", sep = ''))
      write.csv(TT, paste(method_dir, "T.csv", sep = ''))
      write.csv(C, paste(method_dir, "C.csv", sep = ''))
      write.csv(D, paste(method_dir, "D.csv", sep = ''))
      
      # X field
      X_mean <- as.matrix(results_fpls_fem[["X_mean"]], ncol = K, nrow = 1, byrow = TRUE)
      X_hat <- TT %*% t(C) + matrix(X_mean, ncol = K, nrow = N, byrow = TRUE)
      write.csv(X_hat, paste(method_dir, "X_hat.csv", sep = ''))
      write.csv(X_mean, paste(method_dir, "X_mean.csv", sep = ''))
      
      # Response
      Y_mean <- as.matrix(results_fpls_fem[["Y_mean"]], ncol = L, nrow = 1, byrow = TRUE)
      Y_hat <- as.matrix(results_fpls_fem[["fitted.values"]])
      write.csv(Y_hat, paste(method_dir, "Y_hat.csv", sep = ''))
      write.csv(Y_mean, paste(method_dir, "Y_mean.csv", sep = ''))
      
      # Coefficients function
      B_hat <- results_fpls_fem[["coefficient_function"]]
      write.csv(B_hat, paste(method_dir, "B_hat.csv", sep = ''))
      
      # Lambda selection with KCV
      penalty_vec <- c()
      for(x in seq(-10, 1, 0.5))
        penalty_vec <- c(penalty_vec, 10^x)
      
      cv_fem <- cv_seq(X = X_nodes, Y = Y, penalty_vec = penalty_vec,
                       ncomp = 3, folds = 10, basisobj = FEM_basis,
                       method = "r1fpls_fem",
                       verbose = FALSE, stripped = FALSE)
      best_penalties <- cv_fem$best_penalties
      write.csv(best_penalties, paste(method_dir, "lambdas_KCV.csv", sep = ''))
      
      mse <- sum((Y_clean - cv_fem$final_model$fitted.values)^2)/N
      write.csv(mse, paste(method_dir, "mse_KCV.csv", sep = ''))
      
      
      # Solution with Multvariate-NIPALS
      # ||||||||||||||||||||||||||||||||
      
      # Multivariate NIPALS solution
      results_nipals <- PLSR(X_nodes, Y, 3)
      
      method_dir <- paste(sub_test_dir, "nipals", "/", sep = '')
      if (!file.exists(method_dir)){
        dir.create(method_dir)
      }
      
      # FPLS quanities of interest
      W <- as.matrix(results_nipals[["W"]])
      V <- as.matrix(results_nipals[["V"]])
      TT <- as.matrix(results_nipals[["TT"]])
      C <- as.matrix(results_nipals[["C"]])
      D <- as.matrix(results_nipals[["D"]])
      write.csv(V, paste(method_dir, "V.csv", sep = ''))
      write.csv(W, paste(method_dir, "W.csv", sep = ''))
      write.csv(TT, paste(method_dir, "T.csv", sep = ''))
      write.csv(C, paste(method_dir, "C.csv", sep = ''))
      write.csv(D, paste(method_dir, "D.csv", sep = ''))
      
      # X field
      X_mean <- as.matrix(results_nipals[["X.mean"]], ncol = K, nrow = 1, byrow = TRUE)
      X_hat <- as.matrix(results_nipals[["X_hat"]])
      write.csv(X_hat, paste(method_dir, "X_hat.csv", sep = ''))
      write.csv(X_mean, paste(method_dir, "X_mean.csv", sep = ''))
      
      # Response
      Y_mean <- as.matrix(results_nipals[["Y.mean"]], ncol = L, nrow = 1, byrow = TRUE)
      Y_hat <- as.matrix(results_nipals[["Y_hat"]])
      write.csv(Y_hat, paste(method_dir, "Y_hat.csv", sep = ''))
      write.csv(Y_mean, paste(method_dir, "Y_mean.csv", sep = ''))
      
      # Coefficients function
      B_hat <- results_nipals[["Beta"]]
      write.csv(B_hat, paste(method_dir, "B_hat.csv", sep = ''))
      
      
      # Solution with Multvariate-SIMPLS
      # ||||||||||||||||||||||||||||||||
      
      # Multivariate SIMPLS solution
      results_simpls <- SIMPLS(X_nodes, Y, 3, FALSE)
      
      method_dir <- paste(sub_test_dir, "simpls", "/", sep = '')
      if (!file.exists(method_dir)){
        dir.create(method_dir)
      }
      
      # FPLS quanities of interest
      R <- as.matrix(results_simpls[["R"]])
      Q <- as.matrix(results_simpls[["Q"]])
      TT <- as.matrix(results_simpls[["TT"]])
      UU <- as.matrix(results_simpls[["UU"]])
      P <- as.matrix(results_simpls[["P"]])
      V <- as.matrix(results_simpls[["V"]])
      write.csv(R, paste(method_dir, "R.csv", sep = ''))
      write.csv(Q, paste(method_dir, "Q.csv", sep = ''))
      write.csv(TT, paste(method_dir, "T.csv", sep = ''))
      write.csv(UU, paste(method_dir, "U.csv", sep = ''))
      write.csv(P, paste(method_dir, "P.csv", sep = ''))
      write.csv(V, paste(method_dir, "V.csv", sep = ''))
      
      # X field
      X_mean <- as.matrix(results_simpls[["X.mean"]], ncol = K, nrow = 1, byrow = TRUE)
      X_hat <- as.matrix(results_simpls[["X_hat"]])
      write.csv(X_hat, paste(method_dir, "X_hat.csv", sep = ''))
      write.csv(X_mean, paste(method_dir, "X_mean.csv", sep = ''))
      
      # Response
      Y_mean <- as.matrix(results_simpls[["Y.mean"]], ncol = L, nrow = 1, byrow = TRUE)
      Y_hat <- as.matrix(results_simpls[["Y_hat"]])
      write.csv(Y_hat, paste(method_dir, "Y_hat.csv", sep = ''))
      write.csv(Y_mean, paste(method_dir, "Y_mean.csv", sep = ''))
      
      # Coefficients function
      B_hat <- results_simpls[["Beta"]]
      write.csv(B_hat, paste(method_dir, "B_hat.csv", sep = ''))

    }
  }
}
