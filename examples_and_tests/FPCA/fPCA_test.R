

rm(list = ls())
graphics.off()


# ||||||||||||||
# Libraries ----
# ||||||||||||||

library(fdaPDE2)
library(pracma)
library(plot3D)


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

N <- 50
X_clean_locations = NULL
noise_locations = NULL

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


# |||||||||||||||
# fPCA model ----
# |||||||||||||||

cat("Prima pde")

## set model
pde <- new(Laplacian_2D_Order1, mesh_data)

cat("Dopo pde")

## set zero forcing term
quadrature_nodes <- pde$get_quadrature_nodes()
f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
pde$set_forcing_term(as.matrix(f))

## define and init model
model <- new(FPCA_Laplacian_2D_GeoStatLocations, pde)
model$set_locations(locations)
lambda_s <- 10^seq(-4.0, -3.0, by = 0.1)
model$set_lambdas(lambda_s)
model$init_regularization()


## extract principal components
model$set_observations(X_locations)
model$solve()

load1 <- model$loadings()[, 1]
load2 <- model$loadings()[, 2]
load3 <- model$loadings()[, 3]
scor1 <- model$scores()[, 1]
scor2 <- model$scores()[, 2]
scor3 <- model$scores()[, 3]

plot_field(nodes, locations, load1, range(load1), "load1", TRUE)
plot_field(nodes, locations, load2, range(load1), "load2", TRUE)
plot_field(nodes, locations, load3, range(load1), "load3", TRUE)


