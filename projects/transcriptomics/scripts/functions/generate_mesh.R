generate_mesh <- function(locations, concavity, length_threshold, simplify_tol, maximum_area) {
  
  indexes <- locations$ID
  locations <- locations[, c("x_coord", "y_coord")]
  
  # Boundaries
  polygon <- concaveman(as.matrix(locations), concavity = concavity, length_threshold = length_threshold)
  polygon <- data.frame(polygon)
  colnames(polygon) <- c("x_coord", "y_coord")
  
  # Plot boundaries
  # polygon %>%
  #   ggplot() +
  #   geom_point(aes(x_coord, y_coord), locations, color = "grey") +
  #   geom_point(aes(x_coord, y_coord), polygon, color = "blue") +
  #   geom_polygon(aes(x_coord, y_coord), color = "black", fill = "transparent")
  
  # Simplify the polygon by removing closely spaced vertices
  polygon <- Polygon(coords = polygon[,1:2])
  polygon <- Polygons(list(polygon), ID = "my_polygon")
  spatial_polygons <- SpatialPolygons(list(polygon))
  simplified_polygon <- gSimplify(spatial_polygons, tol = simplify_tol)@polygons[[1]]@Polygons[[1]]@coords
  simplified_polygon <- data.frame(simplified_polygon)
  colnames(simplified_polygon) <- c("x_coord", "y_coord")
  simplified_polygon <- simplified_polygon %>% distinct()
  
  # Plot simplified polygon
  # simplified_polygon %>%
  #   ggplot() +
  #   geom_point(aes(x_coord, y_coord), locations, color = "grey") +
  #   geom_point(aes(x_coord, y_coord), simplified_polygon, color = "blue") +
  #   geom_polygon(aes(x_coord, y_coord), color = "black", fill = "transparent")
  
  # Mesh
  np <- dim(simplified_polygon)[1]
  segments <- data.frame(V1 = 1:np, V2 = c(2:np, 1))
  mesh <- create.mesh.2D(simplified_polygon, segments = segments)
  mesh <- refine.mesh.2D(mesh, minimum_angle = 30, maximum_area = maximum_area)
  # plot(mesh,  main = "Mesh", asp = 1)
  # points(locations, pch = 20, col = "blue", cex = 1)
  
  # Remove external locations
  polygon <- Polygon(coords = simplified_polygon)
  polygons <- Polygons(list(polygon), ID = "my_polygon")
  spatial_polygons <- SpatialPolygons(list(polygons))
  spatial_points <- SpatialPoints(locations)
  is_inside <- gCovers(spatial_polygons, spatial_points, byid = TRUE)
  
  return(list(mesh = mesh, final_locations_indexes = indexes[is_inside]))
  
}


save(generate_mesh, file = "scripts/functions/generate_mesh.RData")