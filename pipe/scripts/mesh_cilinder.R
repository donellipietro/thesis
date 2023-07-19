rm(list = ls())
graphics.off()


# ||||||||||||||
# Libraries ----
# ||||||||||||||

library(fdaPDE)

# ||||||||||||||
# Functions ----
# ||||||||||||||

generate_mesh <- function(R, alpha, windings, theta, ds) {
  
  # Length of the pipe
  pitch <- 2*pi*R*tan(alpha)
  L <- round(pitch*windings, 2)
  
  # Helper elements
  N.angles <- 2*pi/theta
  N.circles <- L/ds + 1
  angles <- seq(0, 2*pi, by = theta)[1:N.angles]
  angles_ids <- seq(1, N.angles, by = 1)
  circle <- data.frame(x = R*cos(angles), y = R*sin(angles))
  
  # Nodes
  nodes <- NULL
  for(i in seq(0, N.circles-1, by = 1)){
    nodes <- rbind(nodes, cbind(circle, z = i*ds))
  }
  
  # Triangles
  triangles <- NULL
  for(i in seq(0, N.circles-1, by = 1)){
    angles_ids.this <- angles_ids + i*N.angles
    angles_ids.prev <- angles_ids.this - N.angles
    angles_ids.next <- angles_ids.this + N.angles
    if(i!=0)
      triangles <- rbind(triangles, cbind(angles_ids.this, 
                                          c(angles_ids.this[2:N.angles], angles_ids.this[1]),
                                          c(angles_ids.prev[2:N.angles], angles_ids.prev[1])))
    if(i!=N.circles-1)
      triangles <- rbind(triangles, cbind(angles_ids.this, 
                                          angles_ids.next, 
                                          c(angles_ids.this[2:N.angles], angles_ids.this[1])))
    
  }
  
  # Mesh
  mesh = create.mesh.2.5D(nodes = nodes, triangles = triangles)
  return(mesh)
}

generate_sensors <- function(R, alpha, windings, step) {
  w = 2*pi
  pitch <- 2*pi*R*tan(alpha)
  conv = step/sqrt(R^2*w^2+pitch^2)
  
  t = seq(0, windings, length = 1000)
  K = 0:(windings/conv)
  
  helic <- data.frame(x = R*cos(w*t),
                      y = R*sin(w*t),
                      z = pitch*t)
  sensors <- data.frame(x = R*cos(w*conv*K),
                        y = R*sin(w*conv*K),
                        z = pitch*conv*K)
  
  return(sensors)
}

# |||||||||||||||
# Parameters ----
# |||||||||||||||

# Fiber parameters
step = 0.01
experiment.alpha <- c(10, 20)*(2*pi/360)
experiment.windings = c(19,39) 
experiment.type.number <- length(experiment.alpha)
experiment.names <- c("1", "2")

# Pipe parameters
R <- 38.1*1e-3

# Mesh parameters
theta <- 2*pi/40
ds <- 0.01


# ||||||||||||||||||
# Generate data ----
# ||||||||||||||||||

dir <- paste("../data/fdaPDE_data/", sep = '')
if (!file.exists(dir)){
  dir.create(dir)
}

dir <- paste("../data/fdaPDE_data/mesh/", sep = '')
if (!file.exists(dir)){
  dir.create(dir)
}

for(i in 1:experiment.type.number){
  mesh <- generate_mesh(R, experiment.alpha[i], experiment.windings[i], theta, ds)
  sensors <- generate_sensors(R, experiment.alpha[i], experiment.windings[i], step)
  
  # Plot
  plot(mesh)
  points3d(sensors, color = "purple", size = 10)
  
  # Mesh data
  boundary <- as.numeric(mesh$nodesmarkers)
  edges <- mesh$faces
  elements <- mesh$tetrahedrons
  neigh <- mesh$neighbors
  points <- format(mesh$nodes, nsmall = 5)
  
  dir <-paste(paste("../data/fdaPDE_data/mesh/pipe", experiment.names[i], sep = '_'), "/", sep = "")
  if (!file.exists(dir)){
    dir.create(dir)
  }
  
  # Export mesh
  options(digits=17)
  write.csv(boundary, paste(dir, "boundary.csv", sep = ''))
  write.csv(edges, paste(dir, "edges.csv", sep = ''))
  write.csv(elements, paste(dir, "elements.csv", sep = ''))
  write.csv(neigh, paste(dir, "neigh.csv", sep = ''))
  write.csv(points, paste(dir, "points.csv", sep = ''))
  
  dir <- "../data/fdaPDE_data/locations/"
  if (!file.exists(dir)){
    dir.create(dir)
  }
  
  # Export locatins
  name <- paste("locations", experiment.names[i], sep = '_')
  write.csv(sensors, paste(dir, name, ".csv", sep = ''))
  
}


