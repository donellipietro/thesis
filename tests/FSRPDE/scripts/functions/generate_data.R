generate_data <- function(n.nodes, n.locations, N, data.path) {
  
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
  
  S <- n.locations^2
  x <- seq(0, 1, length.out = n.locations)
  y <- x
  locations <- as.matrix(expand.grid(x, y))
  
  
  ## Field ----
  ## ||||||||||
  
  data_clean = NULL
  noise = NULL
  
  for(n in 1:N){ 
    
    func_evaluation_locations <- generate_X(locations, 2)
    
    noise <- rbind(noise, stats::rnorm(nrow(locations), mean = 0, sd = 0.2))
    data_clean <- rbind(data_clean, func_evaluation_locations)
  }
  
  # Rename cols
  rownames(noise) <- 1:N
  rownames(data_clean) <- 1:N
  
  # Add noise to X-data
  data = data_clean + noise
  
  # Plots
  # plot_field(nodes, locations, data[1,], range(data), "Data at locations, noisy", TRUE)
  
  save(K, S, N,
       nodes, locations, mesh_data,
       data, data_clean,
       file = data.path)
}

save(generate_data, file = "scripts/functions/generate_data.RData")