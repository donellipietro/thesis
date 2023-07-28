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

save(generate_mesh, file = "scripts/functions/generate_mesh.RData")