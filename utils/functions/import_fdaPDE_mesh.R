# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Function to import fdaPDE meshes %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# rm(list = ls())
# graphics.off()
# 
# 
# # ||||||||||||||
# # Libraries ----
# # ||||||||||||||
# 
# library(fdaPDE)
# 
# 
# # ||||||||||||||
# # Functions ----
# # ||||||||||||||
# 
# load("scripts/functions/import_fdaPDE_mesh.RData")
# 
# 
# # ||||||||||||||||
# # Import mesh ----
# # ||||||||||||||||
# 
# meshes.directory <- "../../fdaPDE/test/data/mesh/"
# mesh.name <- "c_shaped" # # "c_shaped", "unit_square", "unit_square_coarse", "unit_square_medium"
# 
# mesh.directory <- paste(meshes.directory, mesh.name, "/", sep = "")
# 
# mesh_import <- import_fdaPDE_mesh(mesh.directory, TRUE)
# mesh <- mesh_import$mesh
# mesh_import$check
# 
# refined_mesh <- refine.mesh.2D(mesh, minimum_angle = 30, maximum_area = 0.007, delaunay = TRUE)
# plot(refined_mesh, asp = 1)


# |||||||||||||
# Function ----
# |||||||||||||

import_fdaPDE_mesh <- function(mesh.directory, PLOT) {

  # Import data
  edges = read.csv(paste(mesh.directory, "edges.csv", sep = ""), header = TRUE)[,2:3]
  points = read.csv(paste(mesh.directory, "points.csv", sep = ""), header = TRUE)[,2:3]
  elements = read.csv(paste(mesh.directory, "elements.csv", sep = ""), header = TRUE)[,2:4]
  boundary = read.csv(paste(mesh.directory, "boundary.csv", sep = ""), header = TRUE)[,2]     # unused
  neigh = read.csv(paste(mesh.directory, "neigh.csv", sep = ""), header = TRUE)[,2:4]         # unused
  
  segments <- edges[boundary[edges[,1]] & boundary[edges[,2]],]

  # Create mesh data
  mesh = create.mesh.2D(nodes = points, segments = segments, triangles = elements)
  
  # Check results
  error <- 0
  error <- error + norm(as.matrix(edges - mesh$edges), "I")
  error <- error + norm(as.matrix(points - mesh$nodes) , "I")
  error <- error + norm(as.matrix(elements - mesh$triangles) , "I")
  error <- error + norm(as.matrix(boundary - mesh$nodesmarkers) , "I")
  error <- error + norm(as.matrix(neigh - mesh$neighbors), "I")
  
  # Plot
  if(PLOT)
    plot(mesh, asp = 1, main = mesh.name)
  
  return(list(mesh = mesh, check = (error == 0), error = error))
}

save(import_fdaPDE_mesh, file = "utils/functions/import_fdaPDE_mesh.RData")


