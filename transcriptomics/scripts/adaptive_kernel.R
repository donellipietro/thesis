# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Function to generate data for FPLSR testing %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()
def.par = par()

# ||||||||||||||
# Libraries ----
# ||||||||||||||

library(plot3D)
library(matlib)
library(car)


# ||||||||||||||
# Functions ----
# ||||||||||||||

# Data generation and visualization
load("functions/plot_field.RData")

# Algorithm
adaptive_kernel <- function(X, V, x_ref, v_ref, N.kn_loc, N.kn_val, max_iter, toll = 1e-6) {
  
  # Initialization
  iter = 0
  convergence <- FALSE
  cov = matrix(c(1,0,0,1), nrow = 2, ncol = 2) # identity matrix
  
  # Difference between each spot and the reference one
  # d[i, ] = x[i, ] - x_ref, for all the spots i = 1, ..., N
  d <- X
  d[,1] <- d[,1] - x_ref[1]
  d[,2] <- d[,2] - x_ref[2]
  
  # L2 distance for the gene expression
  # dist.val[i] = l2norm(v[i,] - v_ref), for all the spots i = 1, ..., N
  dist.val <- c()
  for(i in 1:dim(V)[2]){
    dist.val[i] <- sqrt(sum((V[,i]-v_ref)^2))
  }
  # plot_field(X, X, dist.val, range, "Norm of the values distance")
  # points(interesting_points.loc, col = "green",  pch = 16)
  
  # Iterative algorithm
  while(iter < max_iter & !convergence){
    
    # Inverse of the covariance matrix
    inv_cov <- inv(cov) 
    
    # Distance computation
    # dist.loc[i] = t(d[i,]) %*% inv_cov %*% d[i,]
    # Rmk: if cov = identity matrix 
    #      => inv_cov = identity matrix
    #      => dist.loc[i] = d[,1]^2 + d[,2]^2 = euclidean_distance(x[i,] - x_ref)
    dist.loc <- inv_cov[1,1]*d[,1]^2 + 2*inv_cov[1,2]*d[,1]*d[,2] + inv_cov[2,2]*d[,2]^2
    
    # Identify the clostes N.kn_loc points to the reference spot in terms of location
    order.loc <- order(dist.loc)
    closest_points.loc <- X[order.loc[1:N.kn_loc],] # location of the closest points
    closest_points.val <- dist.val[order.loc[1:N.kn_loc]] # (distance) value of the closest points
    range = range(dist.val)
    plot_field(X, X, dist.val, range, "Spatial dependency")
    # points(closest_points.loc, col = "green",  pch = 16)
    points(x_ref[1], x_ref[2], col = "blue", pch = 16)
    
    # Among the spots just selected using the spatial criteria I identify the
    # clostes N.kn_val spots to the reference spot in terms of values
    order.values <- order(closest_points.val)
    interesting_points.loc <- closest_points.loc[order.values[1:N.kn_val], ]
    # plot_field(X, X, V, range, "Spatial & Value dependency")
    # points(interesting_points.loc, col = "green",  pch = 16)
    # points(x_ref[1], x_ref[2], col = "blue", pch = 16)
    
    # Save the previous covariance matrix
    cov_prev <- cov
    
    # Compute the covariace matrix of the cloud of interesting points
    cov <- cov(interesting_points.loc)
    
    # Convergence check
    if(norm(cov-cov_prev) < toll)
      convergence <- TRUE
    
    iter = iter + 1
    
  }
  
  if(convergence == FALSE)
    cov = matrix(c(1,0,0,1), nrow = 2, ncol = 2)
  
  # Kernel computation
  inv_cov <- inv(cov)
  dist.loc <- inv_cov[1,1]*d[,1]^2 + 2*inv_cov[1,2]*d[,1]*d[,2] + inv_cov[2,2]*d[,2]^2
  ker <- exp(-1/2*dist.loc)
  range.ker <- range(ker)
  # plot_field(X, X, ker, range.ker, "Kenel")
  # ellipse(x_ref, shape=cov, 2, col = 'blue', lty = 2, center.pch = 16)
  
  return(list(ker = ker, cov = cov, interesting_points.loc = interesting_points.loc)) # cov is needed only for the plot
  # return(ker)
}

# |||||||||
# Data ----
# |||||||||

X <- read.csv("../data/DLPFC/DLPFC_xy_coords.csv")  # Coords matrix:   Nx2 matrix (N: number of spots)
V <- read.csv("../data/DLPFC/DLPFC_gene_by_spot_mat_normalized_zero_mean_1stddev.csv", header = TRUE)[2:3640]
Clustering <- t(as.vector(read.csv("../data/DLPFC/DLPFC_Ground_Truth_Labels_Integers.csv", header = FALSE)))     # Features matrix: MxN matrix (M: number of genes) (in this case M = 1)

# |||||||||||||||||||
# Algorithm test ----
# |||||||||||||||||||


# Reference spot (blue spot)
ref_point.index <- round(runif(1, 1, dim(X)[1]))
# ref_point.index <- 2556+1# 2279
x_ref <- as.numeric(X[ref_point.index, ])
v_ref <- V[,ref_point.index]

# Graphic things (no need to translate them in python)
# graphics.off()
# range <- range(Clustering)
# plot_field(X, X, Clustering, range, "Data & reference spot")
# points(x_ref[1], x_ref[2], col = "blue", pch = 16)

# Hyperparameters
N.kn_loc <- 7
N.kn_val <- 3

# Stopping criteria
toll <- 1e-6
max_iter <- 10

# Execution
result <- adaptive_kernel(X, V, x_ref, v_ref, N.kn_loc, N.kn_val, max_iter, toll)
ker <- result$ker
cov <- result$cov
interesting_points.loc <- result$interesting_points.loc

# Final plot
# ||||||||||

# PC computation
ev_dev <- eigen(cov)
evec <- ev_dev$vectors
eval <- ev_dev$values
m1 = evec[2,1]/evec[1,1]
m2 = evec[2,2]/evec[1,2]

# Plot

# par(mfrow = c(1,2), mar = def.par$mar, omi = def.par$omi)

# plot_field(X, X, V, range, "Data", FALSE)
# points(x_ref[1], x_ref[2], col = "blue", pch = 16)

range = range(Clustering)
plot_field(X, X, Clustering, range, "Selected area", TRUE)
points(interesting_points.loc, col = "green",  pch = 16)
points(x_ref[1], x_ref[2], col = "blue", pch = 16)
# abline(a = - m1*x_ref[1] + x_ref[2], b = m1)
# abline(a = - m2*x_ref[1] + x_ref[2], b = m2)
ellipse(x_ref, shape=cov, 2, col = 'blue', lty = 2, center.pch = 16)

# range.ker <- range(ker)
# plot_field(X, X, ker, range.ker, "Kenel", FALSE)
# ellipse(x_ref, shape=cov, 2, col = 'blue', lty = 2, center.pch = 16)

#plot.new()


