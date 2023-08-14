generate_data <- function(N, data.path, options) {
  
  n.nodes = options$n.nodes
  
  ## Domain ----
  ## |||||||||||
  
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
  
  
  ## Locations ----
  ## ||||||||||||||
  
  n.locations <- n.nodes*2
  S <- n.locations^2
  x <- seq(0, 1, length.out = n.locations)
  y <- x
  locations <- as.matrix(expand.grid(x, y))
  mesh_locations <- fdaPDE::create.mesh.2D(locations)
  
  # Prepare list data structure
  mesh_data_locations <- list(
    "nodes"    = mesh_locations$nodes,
    "edges"    = mesh_locations$edges,
    "elements" = mesh_locations$triangles,
    "neigh"    = mesh_locations$neighbors,
    "boundary" = mesh_locations$nodesmarkers
  )
  
  
  ## Eigenfunctions ----
  ## |||||||||||||||||||
  
  C1 <- c(1, 2, 1)
  C2 <- c(1, 1, 2)
  
  load_ex <- list()
  for(i in 1:3){
    load_ex[[i]] <- sin(C1[i]*pi*locations[,1])*sin(C2[i]*pi*locations[,2])
    # plot_field(nodes, locations, eigenfunctions[[i]], range(eigenfunctions[[i]]), paste("Eigenfunction", i), TRUE)
  }
  
  
  # Initialization of a fdaPDE model (needed to initialize R0)
  
  # Set model
  pde <- new(Laplacian_2D_Order1, mesh_data_locations)
  # Set zero forcing term
  quadrature_nodes <- pde$get_quadrature_nodes()
  f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
  pde$set_forcing_term(as.matrix(f))
  # Define and init model
  model <- new(FPCA_Laplacian_2D_GeoStatNodes, pde)
  # Fix lambdas
  lambdas <- 10^seq(-4.8, -3.8, by = 0.01)
  model$set_lambdas(lambdas)
  model$init_regularization()
  
  
  # Eigenfunctions normalization
  R0 = model$R0()
  for(i in 1:3){
    L2normf <- sqrt(as.numeric(t(load_ex[[i]]) %*% R0 %*% load_ex[[i]]))
    load_ex[[i]] <- load_ex[[i]] / L2normf
  }

  
  # Orthogonality check
  # sqrt(abs(as.numeric(t(f[[1]]) %*% R0 %*% f[[2]])))
  # sqrt(abs(as.numeric(t(f[[1]]) %*% R0 %*% f[[3]])))
  # sqrt(abs(as.numeric(t(f[[2]]) %*% R0 %*% f[[3]])))
  
  
  truedatarange <- max(c(load_ex[[1]], load_ex[[2]], load_ex[[3]])) - 
                   min(c(load_ex[[1]], load_ex[[2]], load_ex[[3]]))
  
  ## Field ----
  ## ||||||||||
  
  # Standard deviation on scores and errors
  sd.scores <- c(0.5, 0.3, 0.2)
  sd.error  <- 0.1
  
  # Random scores
  score_ex <- list()
  for(i in 1:3){
    score_ex[[i]] <- rnorm(n = N, sd = sd.scores[i] * truedatarange)
  }
  
  # Clean data
  data_clean <- matrix(score_ex[[1]]) %*% t(matrix(load_ex[[1]])) +
                matrix(score_ex[[2]]) %*% t(matrix(load_ex[[2]])) +
                matrix(score_ex[[3]]) %*% t(matrix(load_ex[[3]]))
  data_centered_clean <- data_clean -
                         matrix(apply(data_clean, 2, mean),
                                ncol = ncol(data_clean),
                                nrow = nrow(data_clean),
                                byrow = TRUE)
  
  # Noise
  noise <- rnorm(n = N * n.locations, sd = sd.error * truedatarange)
  
  # Data
  data <- data_clean + noise
  data_centered <- data - 
                   matrix(apply(data, 2, mean),
                          ncol = ncol(data),
                          nrow = nrow(data),
                          byrow = TRUE)
  
  
  # Plots
  # plot_field(nodes, locations, data_centered[1,], range(data_centered), "Data at locations, noisy", TRUE)
  
  save(K, S, N,
       nodes, locations, mesh_data, R0,
       data, data_centered_clean, data_clean, data_centered,
       score_ex, load_ex,
       file = data.path)
}

save(generate_data, file = "scripts/functions/generate_data.RData")