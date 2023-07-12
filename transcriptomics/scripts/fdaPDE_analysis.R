# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% fdaPDE anlysis on DLPFC data %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()
def.par = par()


# ||||||||||||||
# Libraries ----
# ||||||||||||||

library(pracma)
library(plot3D)

library(fdaPDE2)


# ||||||||||||||
# Functions ----
# ||||||||||||||

load("functions/plot_field.RData")


# |||||||||
# Data ----
# |||||||||

load("../data/data_fdaPDE.RData")


# Domain
# ||||||

nodes = mesh$nodes
mesh_data <- list(
  "nodes"    = mesh$nodes,
  "edges"    = mesh$edges,
  "elements" = mesh$triangles,
  "neigh"    = mesh$neighbors,
  "boundary" = mesh$nodesmarkers
)


# Field
# |||||

# Locations 
locations = as.matrix(X_new)
n_spatial_locations <- nrow(locations)
n_locations = n_spatial_locations

# Field
data = as.matrix(V_new)
plot_field(nodes, locations, data[1,], range(data), "Gene 1", TRUE)


# |||||||||||||||
# fPCA model ----
# |||||||||||||||

# Set model
pde <- new(Laplacian_2D_Order1, mesh_data)
quadrature_nodes <- pde$get_quadrature_nodes()
f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
pde$set_forcing_term(as.matrix(f))

# Define and init model
model <- new(FPCA_Laplacian_2D_GeoStatLocations, pde)
model$set_locations(locations)
lambda_s <- 10^-1 # 10^seq(-4.0, -3.0, by = 0.2)
model$set_lambdas(lambda_s)
model$init_regularization()


## extract principal components
model$set_observations(data)
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


