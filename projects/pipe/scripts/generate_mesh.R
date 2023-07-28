# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Crack Detection Project: Generate Mesh %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()
def.par = par()

cat("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
cat("\n%% Crack Detection Project: Generate Mesh %%")
cat("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&%%%%%%%%%%%%\n")


# ||||||||||||||
# Libraries ----
# ||||||||||||||

library(fdaPDE)


# ||||||||||||||
# Functions ----
# ||||||||||||||

load("scripts/functions/generate_mesh.RData")
load("scripts/functions/generate_sensors.RData")


# |||||||||||||||
# Parameters ----
# |||||||||||||||

# Fiber parameters
step = 0.01
mesh.alpha_vect <- c(10, 20)*(2*pi/360)
mesh.windings_vect = c(39,19) 
mesh.type.number <- length(mesh.alpha_vect)
mesh.names_vect <- c("alpha10", "alpha20")

# Pipe parameters
R <- 38.1*1e-3

# Mesh parameters
theta <- 2*pi/40
ds <- 0.01

# Plots
PLOT = FALSE


# |||||||||
# Mesh ----
# |||||||||

cat("\n")
cat("\n# ||||||||||||||||||||")
cat("\n# Generating mesh ----")
cat("\n# ||||||||||||||||||||\n\n")

dir <- paste("data/mesh/", sep = '')
if (!file.exists(dir)){
  dir.create(dir)
}

for(i in 1:mesh.type.number){
  
  # Mesh
  mesh <- generate_mesh(R, mesh.alpha_vect[i], mesh.windings_vect[i], theta, ds)
  
  # Sensors
  sensors <- projection.points.2.5D(mesh, generate_sensors(R, mesh.alpha_vect[i], mesh.windings_vect[i], step))
  
  # Plot
  if(PLOT){
    plot(mesh)
    points3d(sensors, color = "purple", size = 10)
  }
  
  # Mesh data
  boundary <- as.numeric(mesh$nodesmarkers)
  edges <- mesh$edges
  elements <- mesh$triangles
  neigh <- mesh$neighbors
  points <- format(mesh$nodes, nsmall = 5)
  
  # dir <-paste(paste("../data/fdaPDE_data/mesh/pipe", mesh.names_vect[i], sep = '_'), "/", sep = "")
  # if (!file.exists(dir)){
  #   dir.create(dir)
  # }
  # 
  # # Export mesh
  # options(digits=17)
  # write.csv(boundary, paste(dir, "boundary.csv", sep = ''))
  # write.csv(edges, paste(dir, "edges.csv", sep = ''))
  # write.csv(elements, paste(dir, "elements.csv", sep = ''))
  # write.csv(neigh, paste(dir, "neigh.csv", sep = ''))
  # write.csv(points, paste(dir, "points.csv", sep = ''))
  # 
  # dir <- "../data/fdaPDE_data/locations/"
  # if (!file.exists(dir)){
  #   dir.create(dir)
  # }
  # 
  # # Export locations
  # name <- paste("locations", mesh.names[i], sep = '_')
  # write.csv(sensors, paste(dir, name, ".csv", sep = ''))
  
  save(mesh, sensors, file = paste(dir, "mesh_", mesh.names_vect[i], ".RData", sep = ""))
  
}


