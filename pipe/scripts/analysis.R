# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% fdaPDE anlysis on pipe data %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

experiment.type.names <- c("alpha10","alpha20")
experiments.names <- paste(c(40, 80), "psi", sep = '')

select.type <- 1 # 1: muro, 2: no muro
select.pressure <- 2

# Domain
# ||||||

load(paste("../data/fdaPDE_data/locations/mesh_",select.type,".RData", sep = ""))

nodes = mesh.unrolled$nodes
mesh_data <- list(
  "nodes"    = mesh.unrolled$nodes,
  "edges"    = mesh.unrolled$edges,
  "elements" = mesh.unrolled$triangles,
  "neigh"    = mesh.unrolled$neighbors,
  "boundary" = mesh.unrolled$nodesmarkers
)


# Field
# |||||

# Locations 
locations <- sensors.unrolled[,1:2]
n_spatial_locations <- nrow(locations)
n_locations = n_spatial_locations

# Field
name <- paste(experiment.type.names[select.type], experiments.names[select.pressure], sep = "_")
data <- as.matrix(read.csv(paste("../data/fdaPDE_data/data/X_unrolled_",name,".csv", sep = "")))
data <- data[, -1]
N <- dim(data)[1]
K <- dim(data)[2]
# plot_field(nodes, locations, data[1,], range(data), "Gene 1", TRUE)



# |||||||||||||||||
# fSRPDE model ----
# |||||||||||||||||

## set model
pde <- new(Laplacian_2D_Order1, mesh_data)

## set zero forcing term
quadrature_nodes <- pde$get_quadrature_nodes()
f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
pde$set_forcing_term(as.matrix(f))

## define and init model
model <- new(FSRPDE_Laplacian_2D_GeoStatLocations, pde)
model$set_locations(locations)
lambda_s <- 10^1.1
model$set_lambda_s(lambda_s)


## extract principal components
model$set_observations(data)
model$solve()

f <- model$f()
fitted <- model$fitted()

# plot_field(nodes, nodes, f, range(f), "Denoised data (nodes)", TRUE)
# plot_field(nodes, locations, fitted, range(fitted), "Denoised data (locations)", TRUE)





FEMbasis <- create.FEM.basis(mesh)
f.real <- f[extra==0]


plot(FEM(as.numeric(f.real), FEMbasis))
points3d(sensors, col = "purple", pch = ".", cex = 20)



par(mfrow = c(2,1))
real <- sensors.unrolled[, "extra"]==0
plot(NA, xlim = c(1,sum(real)), ylim = range(data),
     ylab = "Strain", xlab = "Sensor index", type = "l",
     main = name)
grid()
for(i in 1:10){
  points(as.numeric(data[i,real]), type = "l") #, col = rainbow(10)[1])
  # for(j in 2:10){
  #   points(as.numeric(distance), as.numeric(strains[[experiments.names[i]]][j,]), type = "l", col = rainbow(10)[i])
  # }
  # abline(v = c(1, 495, 581, 581+949)*0.01, col = "red", lwd = 2, lty = 2)
}
points(fitted[real], type = "l", col = "blue", lwd = 2)
grid()

plot(fitted[real], xlim = c(1,sum(real)),
     ylab = "Strain", xlab = "Sensor index", type = "l",
     main = "Smoothed data", col = "blue", lwd = 2)
grid()
