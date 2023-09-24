
generate_data <- function(N, data.path, options, options.global)                                
{
  
  n.nodes = options$n.nodes
  
  ## Domain ----
  ## |||||||||||
  
  cat("\n- mesh")
  
  # Set square domain
  K <- n.nodes^2
  x <- seq(0, 1, length.out = n.nodes)
  y <- x
  nodes <- expand.grid(x, y)
  mesh <- fdaPDE::create.mesh.2D(nodes)
  
  # Prepare list data structure
  mesh_data <- list(
    "nodes"    = mesh$nodes,
    "edges"    = mesh$edges,
    "elements" = mesh$triangles,
    "neigh"    = mesh$neighbors,
    "boundary" = mesh$nodesmarkers
  )
  
  cat("\n- mesh_finer")
  
  # Mesh finer
  x <- seq(0, 1, length.out = 100)
  y <- x
  nodes_finer <- expand.grid(x, y)
  mesh_finer <- fdaPDE::create.mesh.2D(nodes_finer)
  
  # FEM basis:
  FEM_basis  <-  fdaPDE::create.FEM.basis(mesh)
  FEM_basis_finer  <-  fdaPDE::create.FEM.basis(mesh_finer)
  
  # Mass matrix:
  R0_finer <- fdaPDE:::CPP_get.FEM.Mass.Matrix(FEMbasis = FEM_basis_finer)
  
  
  # ||||||||||||||
  # Locations ----
  # ||||||||||||||
  
  cat("\n- locations")
  
  # Random locations generation
  # For each element number.locations_per_element locations are randomly picked 
  # inside each element
  number.elements <- dim(mesh$triangles)[1]
  number.locations_per_element <- options.global$number.locations_per_element
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
  
  
  # |||||||||
  # Data ----
  # |||||||||
  
  # Room for data:
  X_clean_nodes = matrix(0, nrow = N, ncol = nrow(nodes))
  X_clean_nodes_finer = matrix(0, nrow = N, ncol = nrow(nodes_finer))
  X_clean_locations = matrix(0, nrow = N, ncol = nrow(locations))
  
  ## X data ----
  ## |||||||||||
  
  cat("\n- X")
  
  set.seed(0)
  X.index <- options$X.index
  
  for(n in 1:N){
    X_clean_nodes[n,] <- generate_X(nodes, X.index)
    X_clean_nodes_finer[n,] <- generate_X(nodes_finer, X.index)
    X_clean_locations[n,] <- generate_X(locations, X.index)
  }
  
  noise_nodes = matrix(stats::rnorm(N*nrow(nodes), mean = 0, sd = 0.2), nrow = N, byrow = TRUE)
  noise_nodes_finer = matrix(stats::rnorm(N*nrow(nodes_finer), mean = 0, sd = 0.2), nrow = N, byrow = TRUE)
  noise_locations = matrix(stats::rnorm(N*nrow(locations), mean = 0, sd = 0.2), nrow = N, byrow = TRUE)
  
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

  
  ## B coefficients function ----
  ## ||||||||||||||||||||||||||||
  
  cat("\n- B")
  
  B.index <- options$B.index
  B <- matrix(generate_B(nodes, B.index), ncol = 1)
  B_finer <- matrix(generate_B(nodes_finer, B.index), ncol = 1)
  
    
  ## Y data ----
  ## |||||||||||
  
  cat("\n- Y\n")
  
  STRATEGY <- options.global$STRATEGY
  NSR.Y <- options.global$NSR.Y
  L <- options.global$L

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

  
  # ||||||||||||||
  # Save data ----
  # ||||||||||||||
  
  data <- list(X_clean_nodes = X_clean_nodes, 
               X_nodes = X_nodes,
               locations = locations,
               X_clean_locations = X_clean_locations, 
               X_locations = X_locations,
               Y_clean = Y_clean,
               Y = Y,
               B = B,
               basisobj = FEM_basis,
               mesh = mesh)
  
  save(data,
       file = data.path)
  
}

save(generate_data, file = "scripts/functions/generate_data_RvsCpp.RData")


