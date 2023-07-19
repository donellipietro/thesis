

rm(list = ls())
graphics.off()
def.par = par()


# ||||||||||||||
# Libraries ----
# ||||||||||||||

library(fdaPDE2)
library(pracma)
library(plot3D)
library(tictoc)



# ||||||||||||||
# Functions ----
# ||||||||||||||

generate_X <- function(points, X.index){
  
  if(X.index == 1){
    a1 = stats::rnorm(1, mean = 1, sd = 0.2)
    a2 = stats::rnorm(1, mean = 1, sd = 0.2)
    z = a1*cos(2*pi*points[,1]) + a2*cos(2*pi*points[,2]) + 1
  }
  if(X.index == 2){
    a1 = stats::rnorm(1, mean = 1, sd = 0.1)
    a2 = stats::rnorm(1, mean = 1, sd = 0.1)
    a3 = stats::rnorm(1, mean = 1, sd = 0.1)
    f <- function(x, y, z = 1){
      coe <- function(x,y){
        1/2*sin(5*pi*x)*exp(-x^2)+1
      }
      a1*sin(2*pi*(coe(y,1)*x*cos(z-2)-a2*y*sin(z-2)))*cos(2*pi*(coe(y,1)*x*cos(z-2+pi/2)+a3*coe(x,1)*y*sin((z-2)*pi/2)))
    }
    # Exact solution (pointwise at nodes)
    z = f(points[,1], points[,2])
  }
  
  return(z)
  
}

plot_field <- function(nodes, locations, data, data.range, title, UNIQUE = TRUE, ONLY_KEY = FALSE) {
  
  if(!ONLY_KEY){
    cols <- hcl.colors(100, "YlOrRd", rev = FALSE)
    
    if(UNIQUE)
      par(mfrow = c(1,1), mar = c(3, 0, 2, 0), oma = c(3,1,1,1))
    
    xlim = range(nodes[,1])
    ylim = range(nodes[,2])
    
    xlim[1] = floor(xlim[1])
    ylim[1] = floor(ylim[1])
    
    xlim[2] = ceiling(xlim[2])
    ylim[2] = ceiling(ylim[2])
    
    data.range[1] = floor(data.range[1])
    data.range[1] = floor(data.range[1])
    
    plot(NA, xlim = xlim, ylim = ylim, asp = 1,
         main = title, bty = "n", xaxt = "n", yaxt = "n",
         ylab = "", xlab = "")
    points(nodes, cex = 0.8, pch = 20, col = "grey") 
    scatter2D(x = locations[, 1], y = locations[, 2], colvar = data, col = cols,
              pch = 16, colkey = FALSE, add = TRUE)
    
    if(UNIQUE)
      colkey(clim = data.range, col = cols, side = 1, add = TRUE, 
             width = 0.5, length = 0.6,
             dist = -0.1)
  }
  else{
    colkey(clim = data.range, col = cols, side = 2, add = FALSE, 
           width = 0.5, length = 0.5,
           dist = 0)
  }
  
}


# |||||||||
# Data ----
# |||||||||


# Domain
# ||||||

# Set square domain
n <- 31
x.2D <- seq(0, 1, length.out = n)
y.2D <- x.2D
nodes <- expand.grid(x.2D, y.2D)
mesh.2D <- fdaPDE::create.mesh.2D(nodes)

# Prepare list data structure
mesh_data <- list(
  "nodes"    = mesh.2D$nodes,
  "edges"    = mesh.2D$edges,
  "elements" = mesh.2D$triangles,
  "neigh"    = mesh.2D$neighbors,
  "boundary" = mesh.2D$nodesmarkers
)


# Data locations
n <- 60
x.2D <- seq(0, 1, length.out = n)
y.2D <- x.2D
locations.2D <- expand.grid(x.2D, y.2D)
mesh.2D <- fdaPDE::create.mesh.2D(locations.2D)
locations <- mesh.2D$nodes
n_spatial_locations <- nrow(locations)
n_locations = n_spatial_locations


# Field
# |||||

N <- 1200
X_clean_locations = NULL
noise_locations = NULL

set.seed(0)

for(n in 1:N){
  
  func_evaluation_locations <- numeric(nrow(locations))
  func_evaluation_locations <- generate_X(locations, 2)
  
  noise_locations <- rbind(noise_locations, stats::rnorm(nrow(locations), mean = 0, sd = 0.2))
  X_clean_locations <- rbind(X_clean_locations, func_evaluation_locations)
}

# Rename cols
rownames(noise_locations) <- 1:N
rownames(X_clean_locations) <- 1:N

# Add noise to X-data
X_locations = X_clean_locations + noise_locations

# Plots
plot_field(nodes, locations, X_locations[1,], range(X_locations), "Data at locations, noisy", TRUE)

save(X_locations, X_clean_locations, nodes, locations, mesh_data, n, N, file = "data_1200.RData")







# |||||||||||||||
# fPCA model ----
# |||||||||||||||

rm(list = ls())
graphics.off()
def.par = par()

library(fdaPDE2)
library(pracma)
library(plot3D)
library(tictoc)


load("data_1200.RData")

## set model
pde <- new(Laplacian_2D_Order1, mesh_data)

## set zero forcing term
quadrature_nodes <- pde$get_quadrature_nodes()
f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
pde$set_forcing_term(as.matrix(f))

lambda_s <- 10^-3

N_vect <- c(75, 100, 150, 200, 300, 400, 600, 800, 1200)
times <- matrix(0, nrow = 3, ncol = length(N_vect))
errors <- matrix(0, nrow = 3, ncol = length(N_vect))

exec <- 1

for(i in 1:length(N_vect)){
  
  cat(paste("\n\nExecution ", exec, ", N = ", N_vect[i], "\n\n", sep = ""))
  
  ## define and init model_FPCA
  rm(model_FPCA)
  model_FPCA <- new(FPCA_Laplacian_2D_GeoStatLocations, pde)
  model_FPCA$set_locations(locations)
  model_FPCA$set_lambdas(lambda_s)
  model_FPCA$init_regularization()
  
  ## define and init model_FPCA_CS
  rm(model_FPCA_CS)
  model_FPCA_CS <- new(FPCA_CS_Laplacian_2D_GeoStatLocations, pde)
  model_FPCA_CS$set_locations(locations)
  model_FPCA_CS$set_lambda_s(lambda_s)
  model_FPCA_CS$init_regularization()
  
  ## define and init model_FPCA_CS + Mass lumping
  rm(model_FPCA_CS_ML)
  model_FPCA_CS_ML <- new(FPCA_CS_Laplacian_2D_GeoStatLocations, pde)
  model_FPCA_CS_ML$set_locations(locations)
  model_FPCA_CS_ML$set_lambda_s(lambda_s)
  model_FPCA_CS_ML$set_mass_lumping(1)
  model_FPCA_CS_ML$init_regularization()
  
  # Sample data
  set.seed(i)
  sample <- sample.int(N, N_vect[i], replace =  TRUE)
  data <- X_locations[sample,]
  data_ex <- X_clean_locations[sample,]
  
  # Set data
  model_FPCA$set_observations(data)
  model_FPCA_CS$set_observations(data)
  model_FPCA_CS_ML$set_observations(data)
  
  # output_file <- paste("execution", exec, ".txt")
  # sink(output_file)
  
  # fPCA
  cat("// fPCA //\n")
  tic()
  model_FPCA$solve()
  elapsed <- toc()
  times[1, i] <- as.numeric(elapsed$toc - elapsed$tic)
  errors[1, i] <- norm(data_ex - model_FPCA$scores() %*% t(model_FPCA$loadings()), "2")/(N_vect[i]*n*n)
  
  # fPCA_CS
  cat("// fPCA - Closed solution //\n")
  tic()
  model_FPCA_CS$solve()
  elapsed <- toc()
  times[2, i] <- as.numeric(elapsed$toc - elapsed$tic)
  errors[2, i] <- norm(data_ex - model_FPCA_CS$scores() %*% t(model_FPCA_CS$loadings()), "2")/(N_vect[i]*n*n)
  
  # fPCA_CS + Mass Lumping
  cat("// fPCA - Closed solution - Mass Lumping //\n")
  tic()
  model_FPCA_CS_ML$solve()
  elapsed <- toc()
  times[3, i] <- as.numeric(elapsed$toc - elapsed$tic)
  errors[3, i] <- norm(data_ex - model_FPCA_CS_ML$scores() %*% t(model_FPCA_CS_ML$loadings()), "2")/(N_vect[i]*n*n)
  
  # close(output_file)
 
}

exec <- exec+1

# load("results_2.RData")

par(mfrow = c(1, 2), mar = def.par$mar, omi = def.par$omi)
matplot(log(N_vect), t(times),
        type = "o", col = c("red", "blue", "green"),
        lty = c(1,1,1), pch = c(1,2,3),
        main = "Execution time",
        ylab = "Execution time",
        xlab = "log(# statistical units)")
grid()
legend("topleft", c("fPCA", "fPCA CS", "fPCA_CS + ML"),
       col = c("red", "blue", "green"), lty = c(1, 1, 1), pch = c(1,2,3))
matplot(log(N_vect), t(errors),
        type = "o", col = c("red", "blue", "green"),
        lty = c(1,1,1), pch = c(1,2,3),
        main = "MSE", 
        ylab = "MSE",
        xlab = "log(# statistical units)")
grid()


save(N_vect, times, errors, X_locations, X_clean_locations, file = "results_2.RData")

