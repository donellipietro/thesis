# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Test fSRPDE Rcpp wrapper %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
# graphics.off()
def.par = par()

cat("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
cat("\n%% Test fSRPDE Rcpp wrapper %%")
cat("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n")


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
load("scripts/functions/generate_data.RData")


# ||||||||||||||||||||
# Test parameters ----
# ||||||||||||||||||||

# Mesh
mesh.path <- "../../fdaPDE/test/data/mesh/unit_square_medium/"

# Number of statistical units
N <- 120

# To force the generation of the field
FORCE_GENERATION = TRUE


# |||||||||
# Data ----
# |||||||||

data.name <- paste("data_", X.index, "_", N, ".RData", sep = "")
data.path <- paste("data/", data.name, sep = "")

if(!file.exists(data.path) | FORCE_GENERATION){
  cat("\n# ||||||||||||||||||||")
  cat("\n# Generating data ----")
  cat("\n# ||||||||||||||||||||\n")
  generate_data(N, mesh.path, data.path, 
                generate_X)
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