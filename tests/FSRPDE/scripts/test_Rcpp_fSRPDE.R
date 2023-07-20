# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Test fSRPDE Rcpp wrapper %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()
def.par = par()

cat("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
cat("\n%% Test fSRPDE Rcpp wrapper %%")
cat("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")


# ||||||||||||||
# Libraries ----
# ||||||||||||||

library(fdaPDE2)
library(pracma)
library(plot3D)


# ||||||||||||||
# Functions ----
# ||||||||||||||

load("../../utils/tests_unit_square.RData")
load("../../utils/plot_field.RData")
load("scripts/functions/generate_data.RData")

# ||||||||||||||||||||
# Test parameters ----
# ||||||||||||||||||||

# Number of nodes (on the side)
n.nodes <- 31

# Number of locations (on the side)
n.locations <- 60

# Number of statistical units
N <- 120


# |||||||||
# Data ----
# |||||||||

data.name <- paste("data_", N, ".RData", sep = "")
data.path <- paste("data/", data.name, sep = "")

if(!file.exists(data.path)){
  cat("\n# ||||||||||||||||||||")
  cat("\n# Generating data ----")
  cat("\n# ||||||||||||||||||||\n")
  generate_data(n.nodes, n.locations, N, data.path)
}

cat("\n# |||||||||||||||||")
cat("\n# Loading data ----")
cat("\n# |||||||||||||||||\n")
load(data.path)


# |||||||||
# Test ----
# |||||||||

cat("\n# |||||||||")
cat("\n# Test ----")
cat("\n# |||||||||\n")

# Set model
pde <- new(Laplacian_2D_Order1, mesh_data)

# Set zero forcing term
quadrature_nodes <- pde$get_quadrature_nodes()
f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
pde$set_forcing_term(as.matrix(f))

# Define and init model
lambda_s <- 1e-4
model <- new(FSRPDE_Laplacian_2D_GeoStatLocations, pde)
model$set_locations(locations)
model$set_lambda_s(lambda_s)


# Set observations
model$set_observations(data)

# Solve
model$solve()

# Results
f <- model$f()
fitted <- model$fitted()
plot_field(nodes, locations, data[1,], range(data[1,]), "Denoised data (locations)", TRUE)
plot_field(nodes, locations, fitted, range(fitted), "Denoised data (locations)", TRUE)
plot_field(nodes, nodes, f, range(f), "Denoised data (nodes)", TRUE)

cat("\n\n\n")