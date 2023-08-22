plot_mesh <- function(mesh, plot_settings) {
  
  # Nodes
  nodes <- data.frame(mesh$nodes)
  colnames(nodes) <- c("x_coord", "y_coord")
  # Triangles
  triangles <- data.frame(mesh$triangles)
  colnames(triangles) <- c("x1", "x2", "x3")
  triangles$id <- 1:dim(triangles)[1]
  triangles <- reshape(triangles,
                       idvar = "id",
                       varying = c("x1", "x2", "x3"),
                       v.names = "index",
                       timevar = "order",
                       times = c(1:3),
                       direction = "long")
  triangles_cordinates <- nodes[triangles$index,]
  triangles = cbind(triangles, triangles_cordinates)
  # Boundary segments
  boundaries <- data.frame(mesh$segments[mesh$segmentsmarkers == 1,])
  boundaries$id <- 1:dim(boundaries)[1]
  boundaries <- reshape(boundaries,
                        idvar = "id",
                        varying = c("X1", "X2"),
                        v.names = "index",
                        timevar = "order",
                        times = c(1:2),
                        direction = "long")
  boundaries_cordinates <- nodes[boundaries$index,]
  boundaries = cbind(boundaries, boundaries_cordinates)
  
  # Mesh
  plot <- plot_settings +
    geom_polygon(aes(x_coord, y_coord, group = id), data = triangles, color = "black", fill = "transparent") +
    geom_polygon(aes(x_coord, y_coord, group = id), data = boundaries, color = "red", fill = "transparent")
  
  return(plot)
  
}

save(plot_mesh, file = "scripts/functions/plot_mesh.RData")