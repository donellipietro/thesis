# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Test fPLSR Rcpp wrapper %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()
def.par = par()

cat("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
cat("\n%% Test fPLSR Rcpp wrapper %%")
cat("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n")


# ||||||||||||||
# Libraries ----
# ||||||||||||||

library(fdaPDE)
library(fdaPDE2)
library(pracma)
library(plot3D)


# ||||||||||||||
# Functions ----
# ||||||||||||||

load("../../utils/functions/tests_unit_square.RData")
load("../../utils/functions/plot_field.RData")
load("../../utils/functions/import_fdaPDE_mesh.RData")
load("../../utils/functions/generate_2d_data.RData")


# ||||||||||||||||||||
# Test parameters ----
# ||||||||||||||||||||

# Mesh
number.locations_per_element <- 1

# Mesh
# "c_shaped",  "unit_square", "unit_square_coarse", "unit_square_medium"
mesh.name <- "unit_square"
mesh_finer.name <- c("unit_square")
mesh.area_refine <- c(0)
directory.test_functions_vec <- c("scripts/functions/tests_unit_square.RData") 

# Noise
STRATEGY = 1
# NSR.X <- (1/3)^2
NSR.Y <- (1/3)^2

# Number of statistical units
N <- 50

# Options
X.index <- 1
data.name <- "test"

# Responses
L <- 1;
B.indexes <- 1



# |||||||||
# Data ----
# |||||||||

data <- generate_2D_data(mesh.name, mesh_finer.name, mesh.area_refine,    # mesh
                         number.locations_per_element,                    # locations
                         N,                                               # samples
                         generate_X, X.index, generate_B, B.indexes,      # tests
                         STRATEGY, NSR.Y,                                 # noise
                         FALSE, FALSE)                                    # plot

# Response
Y_clean <- data[["Y_clean"]]
Y <- data[["Y"]]

# Coefficient function
B <- data[["B"]]

# At nodes
X_clean_nodes <- data[["X_clean_nodes"]]
X_nodes <- data[["X_nodes"]]

# At locations (1 node for each element)
locations <- data[["locations"]]
X_clean_locations <- data[["X_clean_locations"]]
X_locations <- data[["X_locations"]]

# Mesh
FEM_basis <- data[["basisobj"]]
mesh <- data[["mesh"]]

# Prepare list data structure
mesh_data <- list(
  "nodes"    = mesh$nodes,
  "edges"    = mesh$edges,
  "elements" = mesh$triangles,
  "neigh"    = mesh$neighbors,
  "boundary" = mesh$nodesmarkers
)


# |||||||||
# Test ----
# |||||||||

cat("\n# |||||||||")
cat("\n# Test ----")
cat("\n# |||||||||\n\n")

# Set model
pde <- new(Laplacian_2D_Order1, mesh_data)

# Set zero forcing term
quadrature_nodes <- pde$get_quadrature_nodes()
f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
pde$set_forcing_term(as.matrix(f))

# Define and init model
lambda_s <- 10
model <- new(FPLSR_Laplacian_2D_GeoStatLocations, pde)
model$set_verbose(TRUE)
model$set_locations(locations)
model$set_lambda_s(lambda_s)

# Set observations
model$set_data(Y, X_locations)

# Solve
model$solve()

# Results
Y_hat <- model$fitted()
X_hat <- model$reconstructed_field()
plot_field(mesh$nodes, mesh$nodes, data$X_clean_nodes[1,], range(data$X_clean_nodes[1,]), "Clean data (locations)", TRUE)
plot_field(mesh$nodes, mesh$nodes, data$X_nodes[1,], range(data$X_nodes[1,]), "Noisy data (locations)", TRUE)
plot_field(mesh$nodes, mesh$nodes, X_hat[1,], range(X_hat[1,]), "Denoised data (locations)", TRUE)

cat("\n\n\n")
