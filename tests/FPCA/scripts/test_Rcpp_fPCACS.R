# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Test fPCA_CS Rcpp wrapper %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()
def.par = par()

cat("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
cat("\n%% Test fPCA_CS Rcpp wrapper %%")
cat("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")


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

# Set model
pde <- new(Laplacian_2D_Order1, mesh_data)

# Set zero forcing term
quadrature_nodes <- pde$get_quadrature_nodes()
f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
pde$set_forcing_term(as.matrix(f))

# Define and init model
lambda_s <- 10^-3
model <- new(FPCA_CS_Laplacian_2D_GeoStatLocations, pde)
model$set_locations(locations)
model$set_lambda_s(lambda_s)
model$init_regularization()

# Set observations
model$set_observations(data)

# Solve
model$solve()

# Results
load1 <- model$loadings()[, 1]
load2 <- model$loadings()[, 2]
load3 <- model$loadings()[, 3]
scor1 <- model$scores()[, 1]
scor2 <- model$scores()[, 2]
scor3 <- model$scores()[, 3]

plot_field(nodes, locations, data[1,], range(data[1,]), "Data (noisy)", TRUE)
plot_field(nodes, locations, load1, range(load1), bquote(paste('1'^.('st'),' loading function')), TRUE)
plot_field(nodes, locations, load2, range(load1), bquote(paste('2'^.('nd'),' loading function')), TRUE)
plot_field(nodes, locations, load3, range(load1), bquote(paste('3'^.('rd'),' loading function')), TRUE)

cat("\n\n\n")
