# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Crack Detection Project: Analysis %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()
def.par <- par()

cat("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
cat("\n%% Crack Detection Project: Analysis %%")
cat("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")


# ||||||||||||||
# Libraries ----
# ||||||||||||||

library(fdaPDE)
library(fdaPDE2)


# ||||||||||||||
# Functions ----
# ||||||||||||||


# |||||||||||||||
# Parameters ----
# |||||||||||||||

test.name <- "Test4"

experiment.mesh.names <- c("alpha10", "alpha20")
experiments.pressure.names <- paste(c(40, 80), "psi", sep = "")

# Select the tesh
select.mesh <- 1 # 1: muro, 2: no muro

# Select the pressure
select.pressure <- 2

# Select the smoothing parameter
lambda_s <- 10^-3

# |||||||||
# Data ----
# |||||||||

cat("\n")
cat("\n# |||||||||||||||||||||||")
cat("\n# Importing raw data ----")
cat("\n# |||||||||||||||||||||||\n\n")

cat(paste("- Mesh: ", experiment.mesh.names[select.mesh], "\n"))
cat(paste("- Pressure: ", experiments.pressure.names[select.pressure], "\n"))
cat("\n")


# Domain
# ||||||

load(paste("data/mesh/mesh_", experiment.mesh.names[select.mesh], ".RData", sep = ""))

nodes <- mesh$nodes
mesh_data <- list(
  "nodes"    = mesh$nodes,
  "edges"    = mesh$edges,
  "elements" = mesh$triangles,
  "neigh"    = mesh$neighbors,
  "boundary" = mesh$nodesmarkers
)


# Locations
# |||||||||

locations <- sensors[, 1:3]
S <- dim(locations)[1]


# Field
# |||||

name <- paste(experiment.mesh.names[select.mesh], experiments.pressure.names[select.pressure], sep = "_")
data <- as.matrix(read.csv(paste("data/strains/", test.name, "/X_", name, ".csv", sep = "")))
data <- data[, -1]
N <- dim(data)[1]
K <- dim(data)[2]


# |||||||||||||||||
# FRPDE model ----
# |||||||||||||||||

# Set model
pde <- new(Laplacian_Surface_Order1, mesh_data)

# Set zero forcing term
quadrature_nodes <- pde$get_quadrature_nodes()
f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
pde$set_forcing_term(as.matrix(f))

# Define and init model
model <- new(FRPDE_Laplacian_Surface_GeoStatLocations, pde)
model$set_locations(locations)
model$set_lambda_s(lambda_s)

# Set observations
model$set_observations(data)

# Fit the model
model$solve()

# Results
f <- model$f()
fitted <- model$fitted()

# 3D plot
FEMbasis <- create.FEM.basis(mesh)
plot(FEM(as.numeric(f), FEMbasis))
points3d(sensors, col = "purple", pch = ".", cex = 20)

# Smoothing effect
par(mfrow = c(2, 1))
plot(NA,
  xlim = c(1, S), ylim = range(data),
  ylab = "Strain", xlab = "Sensor index", type = "l",
  main = name
)
grid()
for (i in 1:N) {
  points(as.numeric(data[i, ]), type = "l", col = rainbow(N)[i])
}
points(fitted[1, ], type = "l", col = "black", lwd = 4)

plot(fitted[1, ],
  ylab = "Strain", xlab = "Sensor index", type = "l",
  main = "Smoothed data", col = "blue", lwd = 2
)
grid()
