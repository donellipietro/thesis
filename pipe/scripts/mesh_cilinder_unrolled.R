rm(list = ls())
graphics.off()


# ||||||||||||||
# Libraries ----
# ||||||||||||||

library(fdaPDE)

# ||||||||||||||
# Functions ----
# ||||||||||||||

cartesian_to_polar <- function(x, y) {
  r <- sqrt(x^2 + y^2)
  theta <- atan2(y, x)
  return(list(r = r, theta = theta))
}

generate_mesh <- function(R, alpha, windings, theta, ds) {
  
  # Length of the pipe
  pitch <- 2*pi*R*tan(alpha)
  number.elements <- ceiling(pitch*windings/ds)
  L <- number.elements*ds
  
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


unrolling <- function(points, overlap) {
  
  x <- points[,1]
  y <- points[,2]
  z <- points[,3]
  polar <- cartesian_to_polar(x, y)
  r <- polar$r 
  theta <- polar$theta
  
  points.unrolled <- NULL
  
  for(i in 1:dim(points)[1]){
    theta.point <- c(rev(seq(theta[i], overlap[1], by = -2*pi)[-1]),
                     seq(theta[i], overlap[2], by = 2*pi))
    s.point <- r[i] * theta.point
    
    extra <- rep(1, length(s.point))
    extra[which(s.point == r[i]*theta[i])] = 0
    
    points.unrolled <- rbind(points.unrolled,
                             cbind(s = s.point, z = z[i], id = i, extra = extra))
  }
  
  return(points.unrolled)
}

unroll_mesh <- function(mesh, overlap){
  
  nodes <- mesh$nodes
  unrol <- unrolling(nodes, overlap)
  nodes.unrolled <- unrol[,1:2]
  mesh.unrolled <- create.mesh.2D(nodes.unrolled)
  
  return(list(mesh.unrolled = mesh.unrolled, extra = unrol[,4]))
  
}

# |||||||||||||||
# Parameters ----
# |||||||||||||||

# Fiber parameters
step = 0.01
experiment.alpha <- c(10, 20)*(2*pi/360)
experiment.windings = c(39, 19) 
experiment.type.number <- length(experiment.alpha)
experiment.names <- c("1", "2")

# Pipe parameters
R <- 38.1*1e-3

# Mesh parameters
theta <- 2*pi/30
ds <- 0.03


# ||||||||||||||||||
# Generate data ----
# ||||||||||||||||||

overlap = c(-pi/3-pi, pi + pi/3)


dir <- paste("../data/fdaPDE_data/", sep = '')
if (!file.exists(dir)){
  dir.create(dir)
}

dir <- paste("../data/fdaPDE_data/mesh/", sep = '')
if (!file.exists(dir)){
  dir.create(dir)
}

for(i in 1:experiment.type.number){
  
  # Cilindric mesh
  mesh <- generate_mesh(R, experiment.alpha[i], experiment.windings[i], theta, ds)
  sensors <- generate_sensors(R, experiment.alpha[i], experiment.windings[i], step)
  
  # Unrolling
  gen <- unroll_mesh(mesh, overlap)
  mesh.unrolled <- gen$mesh.unrolled
  extra <-  gen$extra
  sensors.unrolled <- unrolling(sensors, overlap)
  
  mesh.range <- range(mesh.unrolled$nodes[,1])
  outside <- which(sensors.unrolled[,1] < mesh.range[1] | sensors.unrolled[,1] > mesh.range[2])
  
  plot(mesh.unrolled, asp = 1)
  points(sensors.unrolled[,1:2], col = "purple", pch = ".", cex = 4)
  points(sensors.unrolled[outside,1:2], col = "red", pch = ".", cex = 4)
  points(sensors.unrolled[sensors.unrolled[,"extra"]==1,1:2], col = "blue", pch = ".", cex = 4)
  abline(v = c(R*pi, -R*pi), col = "red", lty = 2)
  abline(v = R*(overlap), col = "green", lty = 2)
  
  if(length(outside) > 0)
    sensors.unrolled <- sensors.unrolled[-outside, ]
  
  # Mesh data
  boundary <- as.numeric(mesh.unrolled$nodesmarkers)
  edges <- mesh.unrolled$edges
  elements <- mesh.unrolled$triangles
  neigh <- mesh.unrolled$neighbors
  points <- format(mesh.unrolled$nodes, nsmall = 5)
  
  dir <-paste(paste("../data/fdaPDE_data/mesh/pipe_unrolled", experiment.names[i], sep = '_'), "/", sep = "")
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
  name <- paste("locations_unrolled", experiment.names[i], sep = '_')
  write.csv(sensors.unrolled, paste(dir, name, ".csv", sep = ''))
  
  # Export mesh
  name <- paste("mesh", experiment.names[i], sep = '_')
  save(mesh, mesh.unrolled, sensors, sensors.unrolled, extra, file = paste(dir, name, ".RData", sep = ''))
  
}


