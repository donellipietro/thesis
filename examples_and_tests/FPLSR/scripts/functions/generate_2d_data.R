# # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # %% Function to generate data for FPLSR testing %%
# # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# rm(list = ls())
# graphics.off()
# 
# 
# # ||||||||||||||
# # Libraries ----
# # ||||||||||||||
# 
# library(MASS)
# library(plot3D)
# 
# # ||||||||||||||
# # Functions ----
# # ||||||||||||||
# 
# load("scripts/functions/generate_2D_data.RData")
# load("scripts/functions/help_functions.RData")
# 
# # ||||||||||||
# # Example ----
# # ||||||||||||
# 
# set.seed(0)
# 
# # Options
# PLOTS = TRUE
# # Covariates
# X.index <- 1   # covariates to be used
# 
# # Responses
# B.indexes <- c(1, 2)      # coefficients functions to be used
# L <- length(B.indexes);   # number of responses
# 
# # Number of samples
# N <- 50;
# 
# # Locations
# number.locations_per_element <- 1
# 
# # Mesh
# # "c_shaped",  "unit_square", "unit_square_coarse", "unit_square_medium"
# mesh.name <- "c_shaped"
# mesh_finer.name <- ""
# mesh.area_refine <- 0.004
# # "scripts/functions/tests_unit_square.RData", "scripts/functions/tests_c_shaped.RData"
# directory.test_functions <- "scripts/functions/tests_c_shaped.RData"
# load(directory.test_functions)
# 
# # Noise
# STRATEGY = 1
# NSR.X <- (1/3)^2
# NSR.Y <- (1/3)^2
#
# data <- generate_2D_data(mesh.name, mesh_finer.name, mesh.area_refine,          # mesh
#                          number.locations_per_element,                          # locations
#                          N,                                                     # samples
#                          generate_X, X.index, generate_B, B.indexes,            # tests
#                          STRATEGY, NSR.Y,                                       # noise
#                          PLOTS, FALSE)                                           # plots
  
  
# |||||||||||||
# Function ----
# |||||||||||||

generate_2D_data <- function(mesh.name, mesh_finer.name, mesh.area_refine,      # mesh
                             number.locations_per_element,                      # locations
                             N,                                                 # samples
                             generate_X, X.index, generate_B, B.indexes,            # tests
                             STRATEGY, NSR.Y,                                   # noise
                             PLOTS, SHOW)                                       # plots
{
  
  def.par = par()
  
  # |||||||||
  # Mesh ----
  # |||||||||
  
  cat("\n")
  cat("# |||||||||\n")
  cat("# Mesh ----\n")
  cat("# |||||||||\n")
  cat("\n")
  
  meshes.directory <- "../../fdaPDE/test/data/mesh/"
  
  # Mesh
  mesh.directory <- paste(meshes.directory, mesh.name, "/", sep = "")
  mesh_import <- import_fdaPDE_mesh(mesh.directory, FALSE)
  mesh <- mesh_import$mesh
  if(mesh_import$check){
    cat("Mesh imported correctly\n")
  } else{ 
    cat("Error in importing the mesh\n")
  }
  nodes <- mesh$nodes
  
  # Mesh finer
  mesh_finer <- mesh
  if(mesh_finer.name != ""){
    mesh.directory <- paste(meshes.directory, mesh_finer.name, "/", sep = "")
    mesh_import <- import_fdaPDE_mesh(mesh.directory, FALSE)
    mesh_finer <- mesh_import$mesh
    if(mesh_import$check){
      cat("Mesh imported correctly\n")
    } else{
      cat("Error in importing the mesh\n")
    }
  }
  if(mesh.area_refine != 0) {
    mesh_finer <- refine.mesh.2D(mesh_finer, minimum_angle = 30, maximum_area = mesh.area_refine, delaunay = TRUE)
  } else {
    mesh_finer <- refine.by.splitting.mesh.2D(mesh_finer)
  }
  nodes_finer <- mesh_finer$nodes
   
  # Grids
  # x <- seq(0, 1, length.out = nodes.number_on_side)
  # y <- seq(0, 1, length.out = nodes.number_on_side)
  # x_finer <- seq(0, 1, length.out = nodes.number_on_side_finer)
  # y_finer <- seq(0, 1, length.out = nodes.number_on_side_finer)
  
  # Nodes and 2D meshes:
  # nodes <- expand.grid(x = x, y = y)
  # mesh <- fdaPDE::create.mesh.2D(nodes = nodes)
  # nodes_finer <- expand.grid(x = x_finer, y = y_finer)
  # mesh_finer <- fdaPDE::create.mesh.2D(nodes = nodes_finer)
  
  if(PLOTS){
    par(mfrow = c(1,1))
    plot(mesh, asp = 1, main = "Mesh")
    par(mfrow = c(1,1))
    plot(mesh_finer, asp = 1, main = "Mesh, finer")
  }
  
  # FEM basis:
  FEM_basis  <-  fdaPDE::create.FEM.basis(mesh)
  FEM_basis_finer  <-  fdaPDE::create.FEM.basis(mesh_finer)
  
  # Mass matrix:
  R0_finer <- fdaPDE:::CPP_get.FEM.Mass.Matrix(FEMbasis = FEM_basis_finer)
  
  # |||||||||
  # Show ----
  # |||||||||
  
  cat("\n")
  cat("# |||||||||||||||\n")
  cat("# Show tests ----\n")
  cat("# |||||||||||||||\n")
  cat("\n")
  
  if(SHOW){
    show_tests(nodes_finer)
    return()
  }
  
  
  # ||||||||||||||
  # Locations ----
  # ||||||||||||||
  
  cat("\n")
  cat("# |||||||||||||||\n")
  cat("# Locations ----\n")
  cat("# |||||||||||||||\n")
  cat("\n")
  
  # Random locations generation
  # For each element number.locations_per_element locations are randomly picked 
  # inside each element
  number.elements <- dim(mesh$triangles)[1]
  locations <- matrix(0, nrow = number.elements*number.locations_per_element, ncol = 2)
  
  for(elem in 1:number.elements){
    nodes.element <- nodes[mesh$triangles[elem,],]
    a1 <- runif(number.locations_per_element, 0, 1)
    a2 <- runif(number.locations_per_element, 0, 1-a1)
    indexes <- elem*number.locations_per_element - seq(number.locations_per_element-1, 0, length = number.locations_per_element)
    locations[indexes,] <- a1%*%t(as.numeric(nodes.element[1,])) +
                        a2%*%t(as.numeric(nodes.element[2,])) +
                        (1-a1-a2)%*%t(as.numeric(nodes.element[3,]))
  }
  
  if(PLOTS){
    plot(mesh, asp = 1, main = "Locations")
    points(locations, pch = 8, cex = 0.5, col = "blue")
  }
  
  
  # |||||||||
  # Data ----
  # |||||||||
  
  cat("\n")
  cat("# |||||||||\n")
  cat("# Data ----\n")
  cat("# |||||||||\n")
  cat("\n")
  
  # Room for data:
  X_clean_nodes = NULL
  X_clean_nodes_finer = NULL
  X_clean_locations = NULL
  noise_nodes = NULL
  noise_nodes_finer = NULL
  noise_locations = NULL
  
  ## X data ----
  ## |||||||||||
  
  cat("- X field\n")
  
  for(n in 1:N){
    
    func_evaluation_nodes = numeric(nrow(nodes))
    func_evaluation_nodes_finer = numeric(nrow(nodes_finer))
    func_evaluation_locations = numeric(nrow(locations))
    
    func_evaluation_nodes = generate_X(nodes, X.index)
    func_evaluation_nodes_finer = generate_X(nodes_finer, X.index)
    func_evaluation_locations = generate_X(locations, X.index)
    
    noise_nodes = rbind(noise_nodes, stats::rnorm(nrow(nodes), mean = 0, sd = 0.2))
    noise_nodes_finer = rbind(noise_nodes_finer, stats::rnorm(nrow(nodes_finer), mean = 0, sd = 0.2))
    noise_locations = rbind(noise_locations, stats::rnorm(nrow(locations), mean = 0, sd = 0.2))
    X_clean_nodes = rbind(X_clean_nodes, func_evaluation_nodes)
    X_clean_nodes_finer = rbind(X_clean_nodes_finer, func_evaluation_nodes_finer)
    X_clean_locations = rbind(X_clean_locations, func_evaluation_locations)
  }
  
  # Rename cols
  rownames(noise_nodes) <- 1:N
  rownames(noise_nodes_finer) <- 1:N
  rownames(noise_locations) <- 1:N
  rownames(X_clean_nodes) <- 1:N
  rownames(X_clean_nodes_finer) <- 1:N
  rownames(X_clean_locations) <- 1:N
  
  # Add noise to X-data
  X_nodes = X_clean_nodes + noise_nodes
  X_nodes_finer = X_clean_nodes_finer + noise_nodes_finer
  X_locations = X_clean_locations + noise_locations
  
  # X centered
  X_clean_center_nodes <- scale(X_clean_nodes, scale = F)
  X_clean_center_nodes_finer <- scale(X_clean_nodes_finer, scale = F)
  X_clean_center_locations <- scale(X_clean_locations, scale = F)
  X_center_nodes <- scale(X_nodes, scale = F)
  X_center_nodes_finer <- scale(X_nodes_finer, scale = F)
  X_center_locations <- scale(X_locations, scale = F)
  
  data.range = range(X_center_nodes_finer)
  
  
  if(PLOTS){
  
    par(mfrow = c(3, 2), mar = c(1,1,1,1), oma = c(0,0,2,0))
    
    plot_field(nodes, nodes_finer, X_clean_nodes_finer[1,], data.range, "Data at nodes, finer", FALSE)
    plot_field(nodes, nodes_finer, X_nodes_finer[1,], data.range, "Data at nodes, finer, noisy", FALSE)
    
    plot_field(nodes, nodes, X_clean_nodes[1,], data.range, "Data at nodes", FALSE)
    plot_field(nodes, nodes, X_nodes[1,], data.range, "Data at nodes, noisy", FALSE)
    
    plot_field(nodes, locations, X_clean_locations[1,], data.range, "Data at locations", FALSE)
    plot_field(nodes, locations, X_locations[1,], data.range, "Data at locations, noisy", FALSE)
  
  }
  
  
  ## B coefficients function ----
  ## ||||||||||||||||||||||||||||
  
  cat("- B coefficients function\n")
  
  # Room for data:
  B = NULL
  B_finer = NULL
  
  for(i in 1:L){
    B <- cbind(B, generate_B(nodes, B.indexes[i]))
    B_finer <- cbind(B_finer, generate_B(nodes_finer, B.indexes[i]))
  }
  
  if(PLOTS){
    for(i in 1:L){
      data.range = range(B_finer[,i])
      par(mfrow = c(1,1), mar = def.par$mar, oma = def.par$oma)
      plot_field(nodes, nodes_finer, B_finer[,i], data.range,
                 bquote(paste('Coefficients function - B'[.(i)])), TRUE)
    }
  }
    
  ## Y data ----
  ## |||||||||||
  
  cat("- Y data\n")
  
  # No-noise Y:
  if(STRATEGY == 1)
    Y_clean <- as.matrix(X_clean_center_nodes_finer %*% R0_finer %*% B_finer)
  if(STRATEGY == 2)
    Y_clean <-  as.matrix(X_center_nodes_finer %*% R0_finer %*% B_finer)
  
  # Noise:
  # NSR.Y = sigma_noise^2/Var(Y)
  # => sigma_noise^2 = NSR.Y*Var(Y)
  sigma_noise <- sqrt(NSR.Y*diag(var(Y_clean)))
  
  # Noisy Y:
  noise_Y = NULL
  for(i in 1:L){
    noise_Y <- cbind(noise_Y, stats::rnorm(N, mean = 0, sd = sigma_noise[i]))
  }
  Y <- Y_clean + noise_Y
  
  if(PLOTS){
    par(mfrow = c(1, L), mar = def.par$mar, oma = def.par$oma)
    for(i in 1:L){
      plot(Y_clean[,i], Y[,i],
           main =  bquote(paste('Y'[.(i)]*' - NSR = ', .(round(NSR.Y, digits = 2)), sep = "")),
           ylab = bquote(paste('Y'[.(i)])), xlab = bquote(paste('Y'[.(i)]^' clean')),
           asp = 1)
      abline(a = 0, b = 1, col = "red", lty = 2)
      grid()
    }
  }
  
  
  # |||||||||||
  # Return ----
  # |||||||||||
  
  return(list(X_clean_nodes = X_clean_nodes, 
              X_nodes = X_nodes,
              locations = locations,
              X_clean_locations = X_clean_locations, 
              X_locations = X_locations,
              Y_clean = Y_clean,
              Y = Y,
              B = B,
              basisobj = FEM_basis,
              mesh = mesh))
  
}

save(generate_2D_data, file = "scripts/functions/generate_2D_data.RData")


