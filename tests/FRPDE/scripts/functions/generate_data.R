generate_data <- function(N, mesh.path, data.path, generate_X, options = data.frame(), mesh_finer.path = "", mesh.area_refine = 0) {
  
  ## Options ----
  ## ||||||||||||
  
  if("X.index" %in% colnames(options)){
    X.index <- options$X.index
  } else{
    X.index <- 1
  }
  
  if("NSR" %in% colnames(options)){
    NSR.X <- as.numeric(options$NSR)
  } else {
    NSR.X <- (1/3)^2
  }
  
  
  ## Domain ----
  ## |||||||||||
  
  # Mesh
  mesh_import <- import_fdaPDE_mesh(mesh.path, FALSE)
  mesh <- mesh_import$mesh
  if(mesh_import$check){
    # cat("Mesh imported correctly\n")
  } else{ 
    # cat("Error in importing the mesh\n")
  }
  nodes <- mesh$nodes
  K <- dim(nodes)[1]
  
  # Prepare list data structure
  mesh_data <- list(
    "nodes"    = mesh$nodes,
    "edges"    = mesh$edges,
    "elements" = mesh$triangles,
    "neigh"    = mesh$neighbors,
    "boundary" = mesh$nodesmarkers
  )
  
  # Locations
  mesh_finer <- mesh
  if(mesh_finer.path != ""){
    mesh_import <- import_fdaPDE_mesh(mesh_finer.path, FALSE)
    mesh_finer <- mesh_import$mesh
    if(mesh_import$check){
      # cat("Mesh imported correctly\n")
    } else{
      # cat("Error in importing the mesh\n")
    }
  }
  if(mesh.area_refine != 0) {
    mesh_finer <- refine.mesh.2D(mesh_finer, minimum_angle = 30, maximum_area = mesh.area_refine, delaunay = TRUE)
  } else {
    mesh_finer <- refine.by.splitting.mesh.2D(mesh_finer)
  }
  locations <- mesh_finer$nodes
  S <- dim(locations)[1]
  
  
  ## Field ----
  ## ||||||||||
  
  # Data & noise
  set.seed(0)
  func_evaluation_locations <- generate_X(locations, X.index)
  data_clean <- matrix(func_evaluation_locations, nrow = N, ncol = S, byrow = TRUE)
  
  #Noise
  # NSR.X = sigma_noise^2/Var(X)
  # => sigma_noise^2 = NSR.X*Var(X)
  sigma_noise <- sqrt(NSR.X*var(func_evaluation_locations))
  # cat(paste("Sigma used:", sigma_noise))
  set.seed(round(NSR.X*1000))
  noise <- matrix(rnorm(S*N, mean = 0, sd = sigma_noise), nrow = N, ncol = S, byrow = TRUE)
    
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