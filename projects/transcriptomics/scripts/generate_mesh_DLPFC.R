# # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # %% Mesh & data generation for fdaPDE anlysis on DLPFC data %%
# # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# rm(list = ls())
# graphics.off()
# def.par = par()
# 
# 
# # ||||||||||||||
# # Libraries ----
# # ||||||||||||||
# 
# library(ggplot2)
# theme_set(theme_bw())
# 
# 
# library(plot3D)
# library(fdaPDE)
# 
# 
# library(concaveman)
# library(tidyverse)
# 
# 
# # ||||||||||||||
# # Functions ----
# # ||||||||||||||
# 
# # Data generation and visualization
# load("../../utils/functions/plot_field.RData")
# 
# 
# # |||||||||
# # Data ----
# # |||||||||
# 
# data.directory <- "data/SpatialPCA_example_data/"
# list.files(data.directory)
# list.files(paste(data.directory, "DLPFC/", sep = ""))
# load(paste(data.directory, "DLPFC/LIBD_sample9.RData", sep = ""))
# location <- xy_coords
# xy_coords <- data.frame(location)
# colnames(xy_coords) <- c("x_coord", "y_coord")

# X <-  xy_coords
# V <- data
# Clustering <- as.numeric(Layer_sub)

generate_mesh_DLPFC <- function(X, V, Clustering, concavity, length_threshold, simplify_tol, maximum_area){

  # ||||||||||||||||||||||
  # Identify boundary ----
  # ||||||||||||||||||||||
  
  max_iter <- 150
  
  iter <- 0
  convergence <- FALSE
  
  inside <- c()
  boundary <- c(1105)
  
  while(iter < max_iter & !convergence){
    
    boundary_prev <- boundary
    
    for(el in boundary){
      
      x_ref <- as.numeric(X[el,])
      
      d <- X
      d[, 1] <- d[, 1] - x_ref[1]
      d[, 2] <- d[, 2] - x_ref[2]
      
      dist <- d[, 1]^2 + d[, 2]^2
      
      neigh_el <- which(dist < 50 & dist > 0)
      
      neigh_el.number <- length(neigh_el)
      
      if(neigh_el.number == 6){
        boundary <- boundary[boundary != el]
        inside <- union(inside, el)
      }
      
      boundary <- union(boundary, neigh_el)
      
    }
    
    boundary <- boundary[!(boundary %in% inside)]
    
    iter <- iter + 1
    
    if(setequal(boundary, boundary_prev)){
      convergence <- TRUE
    }
  
  }
  
  # ggplot() +
  # geom_point(aes(x_coord, y_coord), X, color = "black") +
  # geom_point(aes(x_coord, y_coord), X[inside,], color = "green") +
  # geom_point(aes(x_coord, y_coord), X[boundary,], color = "red")
  
  
  # ||||||||||||||||||||
  # Remove outliers ----
  # ||||||||||||||||||||
  
  type <- rep(0, dim(X)[1])
  type[inside] <- 1
  type[boundary] <- 2
  
  X_new <- X[type != 0,]
  V_new <- V[, type != 0]
  Clustering_new <- Clustering[type != 0]
  type <- type[type != 0]
  inside <- which(type == 1)
  boundary <- which(type == 2)
  
  # ggplot() +
  #   geom_point(aes(x_coord, y_coord), X_new, color = "black") +
  #   geom_point(aes(x_coord, y_coord), X_new[inside,], color = "green") +
  #   geom_point(aes(x_coord, y_coord), X_new[boundary,], color = "red")
  # 
  # 
  # plot_field(X_new, X_new, Clustering_new, c(1,8), "Ground truth", TRUE)
  
  # ||||||||||||||||||||
  # Refine boundary ----
  # ||||||||||||||||||||
  
  R <- 3
  N <- 4
  D <- 20
  
  out_boundary <- data.frame(matrix(ncol = dim(X_new)[2], nrow = 0))
  colnames(out_boundary) <- colnames(X_new)
  
  for(el in boundary){
    
    x_ref <- as.numeric(X_new[el,])
    
    d <- X_new
    d[, 1] <- d[, 1] - x_ref[1]
    d[, 2] <- d[, 2] - x_ref[2]
    
    dist <- d[, 1]^2 + d[, 2]^2
    
    neigh_el <- which(dist < 200 & dist > 0)
    
    angles <- seq(0, 2*pi*(N-1)/N, by = 2*pi/N)
    new_points <- R*data.frame(x_coord = cos(angles), y_coord = sin(angles))
    new_points[,1] <- new_points[,1] + x_ref[1]
    new_points[,2] <- new_points[,2] + x_ref[2]
    
    new_points_ok <- data.frame(matrix(ncol = dim(X)[2], nrow = 0))
    colnames(new_points_ok) <- colnames(X_new)
    
    for(i in 1:dim(new_points)[1]){
      
      x_ref <- as.numeric(new_points[i,])
      
      d <- rbind(X_new[neigh_el, 1:2],
                 out_boundary)
      d[, 1] <- d[, 1] - x_ref[1]
      d[, 2] <- d[, 2] - x_ref[2]
      
      dist <- d[, 1]^2 + d[, 2]^2
      
      if(min(dist) > D){
        new_points_ok <- rbind(new_points_ok, new_points[i,])
      }
      
    }
    
    out_boundary <- rbind(out_boundary, new_points_ok)
    
  }
  
  
  # ggplot() +
  #   geom_point(aes(x_coord, y_coord), X_new, color = "black") +
  #   geom_point(aes(x_coord, y_coord), X_new[inside,], color = "green") +
  #   geom_point(aes(x_coord, y_coord), X_new[boundary,], color = "red") +
  #   geom_point(aes(x_coord, y_coord), out_boundary, color = "blue")
  
  
  
  # ||||||||||||
  # Polygon ----
  # ||||||||||||
  
  polygons <- concaveman(as.matrix(out_boundary), concavity = concavity, length_threshold = length_threshold)
  polygons <- data.frame(polygons)
  colnames(polygons) <- colnames(X_new[,1:2])
  
  polygon <- Polygon(coords = polygons[,1:2])
  polygon <- Polygons(list(polygon), ID = "my_polygon")
  spatial_polygons <- SpatialPolygons(list(polygon))
  simplified_polygon <- gSimplify(spatial_polygons, tol = simplify_tol)@polygons[[1]]@Polygons[[1]]@coords
  simplified_polygon <- data.frame(simplified_polygon)
  colnames(simplified_polygon) <- c("x_coord", "y_coord")
  simplified_polygon <- simplified_polygon %>% distinct()
  
  
  # simplified_polygon %>%
  #   ggplot() +
  #   geom_polygon(aes(x_coord, y_coord), color = "black", fill = "transparent") +
  #   geom_point(aes(x_coord, y_coord), X_new, color = "grey") +
  #   geom_point(aes(x_coord, y_coord), X_new[inside,], color = "green") +
  #   geom_point(aes(x_coord, y_coord), X_new[boundary,], color = "red") +
  #   geom_point(aes(x_coord, y_coord), simplified_polygon, color = "blue") 
  
  
  
  polygons_new <- simplified_polygon %>% distinct()
  
  
  # ||||||||||||||||||||
  # Holes Detection ----
  # ||||||||||||||||||||
  
  
  residual <- anti_join(out_boundary, polygons)
  
  holes <- data.frame(matrix(ncol = dim(X_new)[2], nrow = 0))
  colnames(holes) <- colnames(X_new)
  
  for(i in 1:dim(residual)[1]){
    
    x_ref <- as.numeric(residual[i,])
    
    d <- polygons[,]
    d[, 1] <- d[, 1] - x_ref[1]
    d[, 2] <- d[, 2] - x_ref[2]
    
    dist <- d[, 1]^2 + d[, 2]^2
    
    if(min(dist) > 40){
      holes <- rbind(holes,
                     residual[i,])
    }
    
  }
  
  
  # polygons %>%
  #   ggplot() +
  #   geom_point(aes(x_coord, y_coord), holes, color = "orange")
  
    
  dist <- dist(holes, method='euclidean')
  hier <- hclust(dist, "average")
  
  # plot(hier)
  # rect.hclust(hier, k=4)
  
  cluster = cutree(hier, k=4)
  
  # polygons %>%
  #   ggplot() +
  #   geom_point(aes(x_coord, y_coord), holes, color = rainbow(4)[cluster])
  # 
  
  clustered_holes <- list()
  
  for(i in 1:4){
    clustered_holes[[i]] <- holes[cluster==i,]
    
    clustered_holes[[i]] <- concaveman(as.matrix(clustered_holes[[i]]), concavity = 2, length_threshold = 1)
    clustered_holes[[i]] <- data.frame(clustered_holes[[i]])
    colnames(clustered_holes[[i]]) <- colnames(X_new[,1:2])
    
    clustered_holes[[i]] <- clustered_holes[[i]] %>% distinct()
    
  }
  
  polygons <- polygons_new
  
  # polygons %>%
  #   ggplot() +
  #   geom_polygon(aes(x_coord, y_coord), color = "black", fill = "transparent") +
  #   geom_point(aes(x_coord, y_coord), X_new, color = "grey") +
  #   geom_point(aes(x_coord, y_coord), X_new[inside,], color = "green") +
  #   geom_point(aes(x_coord, y_coord), X_new[boundary,], color = "red") +
  #   geom_point(aes(x_coord, y_coord), polygons, color = "blue") +
  #   geom_polygon(aes(x_coord, y_coord), data = clustered_holes[[1]],   color = "black", fill = "transparent")+
  #   geom_point(aes(x_coord, y_coord), data = clustered_holes[[1]],   color = "orange")+
  #   geom_polygon(aes(x_coord, y_coord), data = clustered_holes[[3]],   color = "black", fill = "transparent") +
  #   geom_point(aes(x_coord, y_coord), data = clustered_holes[[3]],   color = "orange")
  
  
  # |||||||||
  # Mesh ----
  # |||||||||
  
  np <- dim(polygons)[1]
  nc1 <- dim(clustered_holes[[1]])[1]
  nc2 <- dim(clustered_holes[[3]])[1]
  
  nodes <-rbind(polygons, 
                clustered_holes[[1]],
                clustered_holes[[3]]) #,
                # X_new)
  segments <- rbind(data.frame(V1 = 1:np, V2 = c(2:np, 1)),
                    data.frame(V1 = (np+1):(np+nc1), V2 = c((np+2):(np+nc1), (np+1))), 
                    data.frame(V1 = (np+nc1+1):(np+nc1+nc2), V2 = c((np+nc1+2):(np+nc1+nc2), (np+nc1+1))))
                    
  holes <- rbind(colMeans(clustered_holes[[1]]), 
                colMeans(clustered_holes[[3]]))
  
  mesh <- create.mesh.2D(nodes, segments = segments, holes = holes)
  mesh <- refine.mesh.2D(mesh, minimum_angle = 30, maximum_area = maximum_area)
  # plot(mesh,  main = "Mesh", asp = 1)
  # points(X_new, pch = 20, col = "blue", cex = 1)
  
  return(list(mesh = mesh,
              data = V_new,
              xy_coords = X_new,
              target = Clustering_new))

}

save(generate_mesh_DLPFC, file = "scripts/functions/generate_mesh_DLPFC.RData")
