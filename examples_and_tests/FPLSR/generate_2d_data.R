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

save(generate_2d_data, file = "generate_2d_data.RData")