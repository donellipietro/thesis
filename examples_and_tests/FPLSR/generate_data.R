# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Utility to generate data for FPLSR testing %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()


# |||||||||||||||||||||
# Import Libraries ----
# |||||||||||||||||||||

library(pls)
library(fda)
library(pR1FPLS)


# ||||||||||||||||||||||||||
# Environment variables ----
# ||||||||||||||||||||||||||

complete = FALSE
tests_dir <- "tests/"


# ||||||||||||||
# Functions ----
# ||||||||||||||

generate_2d_data <- function(x, y, num_samples = 100, beta_num = 3, Rsq = 0.95) {
  
  set.seed(0)
  
  # Nodes and 2D mesh:
  nodes <- expand.grid(x = x, y = y)
  mesh <- fdaPDE::create.mesh.2D(nodes = nodes)
  
  # FEM basis:
  FEM_basis  <-  fdaPDE::create.FEM.basis(mesh)
  
  # Mass matrix:
  R0 <- fdaPDE:::CPP_get.FEM.Mass.Matrix(FEMbasis = FEM_basis)
  
  # Subgrid
  x_sub <- seq(0.1, 0.9, length.out = 9)
  y_sub <- seq(0.1, 0.9, length.out = 9)
  locations_sub <- expand.grid(x = x_sub, y = y_sub)
  
  # Room for data:
  locations = NULL
  X_clean_nodes = NULL
  X_clean_locations = NULL
  X_clean_locations_sub = NULL
  noise = NULL
  noise_sub = NULL
  
  # random locations generation
  ds = 1/length(x)
  Sigma = (ds/3)^2*matrix(c(1,0,0,1), 2, 2)
  dx = mvrnorm(dim(mesh$nodes)[1], mu = c(0,0), Sigma = Sigma);
  locations = mesh$nodes + dx
  # manage points outside the domain
  locations = abs(locations)
  locations[locations > 1] = 1-(locations[locations > 1]-1)
  
  # generate X field
  for(ii in 1:num_samples){
    a1 = stats::rnorm(1, mean = 1, sd = 0.2)
    a2 = stats::rnorm(1, mean = 1, sd = 0.2)
    
    func_evaluation_nodes = numeric(nrow(mesh$nodes))
    func_evaluation_locations = numeric(nrow(locations))
    func_evaluation_locations_sub = numeric(nrow(locations_sub))
    
    for (i in 1:nrow(mesh$nodes)){
      func_evaluation_nodes[i] = a1* cos(2*pi*mesh$nodes[i,1]) +
        a2* cos(2*pi*mesh$nodes[i,2]) + 1
      func_evaluation_locations[i] = a1* cos(2*pi*locations[i,1]) +
        a2* cos(2*pi*locations[i,2]) + 1
    }
    
    for (i in 1:nrow(locations_sub)){
      func_evaluation_locations_sub[i] = a1* cos(2*pi*locations_sub[i,1]) +
        a2* cos(2*pi*locations_sub[i,2]) + 1
    }
    
    noise = rbind(noise, stats::rnorm(nrow(mesh$nodes), mean = 0, sd = 0.2))
    noise_sub = rbind(noise_sub, stats::rnorm(nrow(locations_sub), mean = 0, sd = 0.2))
    X_clean_nodes = rbind(X_clean_nodes, func_evaluation_nodes)
    X_clean_locations = rbind(X_clean_locations, func_evaluation_locations)
    X_clean_locations_sub = rbind(X_clean_locations_sub, func_evaluation_locations_sub)
  }
  
  # adding noise to X-data
  X_nodes = X_clean_nodes + noise
  X_locations = X_clean_locations + noise
  X_locations_sub = X_clean_locations_sub + noise_sub
  
  # Generate beta(x, y):
  if (beta_num == 1) {
    # centered at (0.5, 0.5):
    r  <-  0.4 # set the r parameter
    z  <-  5*exp(-((nodes[, 1] - 0.5 )^2 + (nodes[, 2] - 0.5 )^2)/( 2*r^2 ))
  }else if (beta_num == 2) {
    # top right corner
    z <- 5*exp(-((nodes[, 1] - 0.75)^2 + (nodes[, 2] - 0.75)^2)/( 2*0.2^2 ))
  }else if (beta_num == 3) {
    # bottom left + top right corner:
    z <- 5*exp(-((nodes[, 1] - 0.75)^2 + (nodes[, 2] - 0.75)^2)/( 2*0.25^2 )) +
      5*exp(-((nodes[, 1] - 0.1)^2 + (nodes[, 2] - 0.1)^2)/( 2*0.25^2 ))
  }else if (beta_num == 4) {
    # semi-circumference:
    r  <-  0.2 # set the r parameter
    z  <-  5*exp(-((nodes[, 1] - 0.5 )^2 + (nodes[, 2] )^2)/( 2*r^2 ))
  }else if (beta_num == 5) {
    # monkey saddle
    z = ((nodes[, 1]*4 - 2)^3 - 3*(nodes[, 1]*4 - 2)*((nodes[, 2]*4 - 2)^2))
  }else if (beta_num == 6) {
    # Test 1 - fdaPDE
    f = function(x, y, z = 1)
    {
      coe = function(x,y) 1/2*sin(5*pi*x)*exp(-x^2)+1
      sin(2*pi*(coe(y,1)*x*cos(z-2)-y*sin(z-2)))*cos(2*pi*(coe(y,1)*x*cos(z-2+pi/2)+coe(x,1)*y*sin((z-2)*pi/2)))
    }
    # Exact solution (pointwise at nodes)
    z = f(nodes[,1], nodes[,2])
  }
  
  B <- as.matrix(z)
  
  # X centered
  X_clean_center_nodes <- scale(X_clean_nodes, scale = F)
  X_clean_center_locations <- scale(X_clean_locations, scale = F)
  X_clean_center_locations_sub <- scale(X_clean_locations_sub, scale = F)
  X_center_nodes <- scale(X_nodes, scale = F)
  X_center_locations <- scale(X_locations, scale = F)
  X_center_locations_sub <- scale(X_locations_sub, scale = F)
  
  # No-noise Y:
  Y_clean <- as.matrix(X_center_nodes %*% R0 %*% B)
  
  # Variance of errors:
  var_e <- (1/Rsq - 1)*stats::var(Y_clean)
  
  # Noisy Y:
  Y <- Y_clean + as.matrix(stats::rnorm(length(Y_clean), mean = 0, sd = sqrt(var_e)))
  
  return(list(X_clean_nodes = X_clean_nodes, 
              X_nodes = X_nodes,
              locations = locations,
              X_clean_locations = X_clean_locations, 
              X_locations = X_locations,
              locations_sub = locations_sub,
              X_clean_locations_sub = X_clean_locations_sub, 
              X_locations_sub = X_locations_sub,
              Y_clean = Y_clean,
              Y = Y,
              B = B,
              basisobj = FEM_basis,
              mesh = mesh))
  
}


# ||||||||||||||||||
# Generate data ----
# ||||||||||||||||||

if (!file.exists(tests_dir)){
dir.create(tests_dir)
}

# Grid
x <- seq(0, 1, length.out = 60)
y <- seq(0, 1, length.out = 60)

# Directories
test_dirs <- c()
for(i in 1:5){
  test_dirs[i] <- paste(tests_dir, "2D_test", i-1, "/", sep = '')
  if (!file.exists(test_dirs[i])){
    dir.create(test_dirs[i])
  }
}

for(i in 1:6){

  data <- generate_2d_data(x, y, 60, i, 0.95)
  
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

  # X field
  write.csv(X_clean_nodes, paste(sub_test_dir, "/X_clean.csv", sep = ''))
  write.csv(X_nodes, paste(sub_test_dir, "/X.csv", sep = ''))

  # Response
  Y_hat <- as.matrix(results_fpls_fem[["fitted.values"]])
  write.csv(Y_clean, paste(sub_test_dir, "/Y_clean.csv", sep = ''))
  write.csv(Y, paste(sub_test_dir, "/Y.csv", sep = ''))
  write.csv(Y_hat, paste(sub_test_dir, "/Y_hat.csv", sep = ''))
  
  # Coefficient function
  B_hat <- results_fpls_fem[["coefficient_function"]]
  write.csv(B, paste(sub_test_dir, "/B.csv", sep = ''))
  write.csv(B_hat, paste(sub_test_dir, "/B_hat.csv", sep = ''))
  
  # FPLS quanities of interest
  W <- as.matrix(results_fpls_fem[["W"]])
  TT <- as.matrix(results_fpls_fem[["TT"]])
  C <- as.matrix(results_fpls_fem[["C"]])
  D <- as.matrix(results_fpls_fem[["D"]])
  write.csv(W, paste(sub_test_dir, "/W.csv", sep = ''))
  write.csv(TT, paste(sub_test_dir, "/T.csv", sep = ''))
  write.csv(C, paste(sub_test_dir, "/C.csv", sep = ''))
  write.csv(D, paste(sub_test_dir, "/D.csv", sep = ''))
  
  
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
