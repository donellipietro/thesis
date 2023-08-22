# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Spatial Transcriptomics Project: Generate mesh %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()

cat("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
cat("\n%% Spatial Transcriptomics Project: Generate mesh %%")
cat("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&%%%%%%%%%%%%%%%%%%%%\n")


# ||||||||||||||
# Libraries ----
# ||||||||||||||

library(plot3D)
library(gridExtra)
library(ggplot2)
theme_set(theme_bw())

library(SpatialPCA)

# Outliers
library(dplyr)
library(FNN)

# Boundary
library(concaveman)
library(tidyverse)
library(sp)
library(rgeos)

# Mesh
library(fdaPDE)


# ||||||||||||||
# Functions ----
# ||||||||||||||

load("scripts/functions/plot_mesh.RData")
load("scripts/functions/generate_mesh.RData")


# |||||||||||||||
# Parameters ----
# |||||||||||||||

data.directory <- "data/SpatialPCA_example_data/"

plot_settings <-  ggplot() + 
                  coord_fixed() +
                  xlab(NULL) + ylab(NULL) +
                  theme(plot.title = element_text(hjust = 0.5, size = 20), 
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        panel.border = element_blank(),
                        axis.text.x=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks=element_blank())

images.directory <- "images/mesh_generation/"
if (!file.exists(images.directory)){
  dir.create(images.directory)
}

processed_data.directory <- "data/processed/"
if (!file.exists(processed_data.directory)){
  dir.create(processed_data.directory)
}


# ||||||||||
# DLPFC ----
# ||||||||||

cat("\n")
cat("\n# ||||||||||")
cat("\n# DLPFC ----")
cat("\n# ||||||||||\n\n")

rm(list = ls()[! ls() %in% c("plot_field", "plot_mesh", "plot_settings", "images.directory", "data.directory", "processed_data.directory", "generate_mesh")])

# Load data
cat("- Loading raw data \n")
list.files(data.directory)
list.files(paste(data.directory, "DLPFC/", sep = ""))
load(paste(data.directory, "DLPFC/LIBD_sample9.RData", sep = ""))
xy_coords$ID = 1:dim(xy_coords)[1]

# Na identification
cat("- Na identification \n")
na.indexes <- is.na(Layer_sub)

# "Islands" identification
cat("- Islands identification \n")
k <- 5
knn_result <- get.knn(xy_coords[,1:2], k)
xy_coords$avg_distance <- rowMeans(knn_result$nn.dist)
threshold <- 1.6 * mean(xy_coords$avg_distance)
islands.indexes <- xy_coords$avg_distance > threshold

# Mesh
cat("- Meshing \n")
indexes <- !na.indexes & !islands.indexes
locations <- xy_coords[indexes, c("x_coord", "y_coord", "ID")]
mesh_gen <- generate_mesh(locations,
                          concavity = 1,
                          length_threshold = 15,
                          simplify_tol = 0.1,
                          maximum_area = 80)
mesh <- mesh_gen$mesh
final_locations_indexes <- mesh_gen$final_locations_indexes

cat("- Plots \n")

# Plot mesh
p1 <- plot_mesh(mesh, plot_settings) +
      ggtitle("Mesh")
p2 <- p1 +
      geom_point(aes(x_coord, y_coord), xy_coords[final_locations_indexes, ], color = "blue") +
      ggtitle("Mesh and final locations")
plot <- grid.arrange(p1, p2, ncol=2)
ggsave(paste(images.directory, "DLPFC_mesh", ".jpg", sep = ""),
       plot = plot, width = 16, height = 8, dpi = 300)

# Plot locations
p1 <- plot_settings +
      geom_point(aes(x_coord, y_coord), xy_coords, color = "purple") + 
      geom_point(aes(x_coord, y_coord), xy_coords[final_locations_indexes, ], color = "grey") +
      geom_point(aes(x_coord, y_coord), xy_coords[islands.indexes, ], color = "red") +
      geom_point(aes(x_coord, y_coord), xy_coords[na.indexes, ], color = "orange") + 
      ggtitle("Original locations")
p2 <- plot_settings +
      geom_point(aes(x_coord, y_coord), xy_coords[final_locations_indexes,], color = "grey") +
      ggtitle("Final locations")
plot <- grid.arrange(p1, p2, ncol=2)
ggsave(paste(images.directory, "DLPFC", ".jpg", sep = ""),
       plot = plot, width = 16, height = 8, dpi = 300)


# Ground truth
target <- as.numeric(Layer_sub)
locations <- xy_coords[final_locations_indexes, c("x_coord", "y_coord")]
# plot_field(locations, locations, target, range(target),
#            "DLPFC - LIBD_sample9", TRUE)

# Export data
cat("- Exporting data \n")
save(mesh, locations,
     xy_coords, final_locations_indexes,
     target, 
     file = paste(processed_data.directory, "DLPFC_mesh", ".RData", sep = ""))


# |||||||||||||||||||||||||
# Slide-seq Cerebellum ----
# |||||||||||||||||||||||||

cat("\n")
cat("\n# |||||||||||||||||||||||||")
cat("\n# Slide-seq Cerebellum ----")
cat("\n# |||||||||||||||||||||||||\n\n")

rm(list = ls()[! ls() %in% c("plot_field", "plot_mesh", "plot_settings", "images.directory", "data.directory", "processed_data.directory", "generate_mesh")])

# Load data
cat("- Loading raw data \n")
list.files(data.directory)
list.files(paste(data.directory, "SlideseqCerebellum/", sep = ""))
load(paste(data.directory, "SlideseqCerebellum/slideseq.rds", sep = ""))
xy_coords <- data.frame(location)
colnames(xy_coords) <- c("x_coord", "y_coord")
xy_coords$ID = 1:dim(xy_coords)[1]

# "Islands" identification
cat("- Islands identification \n")
k <- 50
knn_result <- get.knn(xy_coords[,1:2], k)
xy_coords$avg_distance <- rowMeans(knn_result$nn.dist)
threshold <- 2 * mean(xy_coords$avg_distance)
islands.indexes <- xy_coords$avg_distance > threshold

# Mesh
cat("- Meshing \n")
indexes <- !islands.indexes
locations <- xy_coords[indexes, c("x_coord", "y_coord", "ID")]
mesh_gen <- generate_mesh(locations,
                          concavity = 2,
                          length_threshold = 200,
                          simplify_tol = 17,
                          maximum_area = 6000)
mesh <- mesh_gen$mesh
final_locations_indexes <- mesh_gen$final_locations_indexes

cat("- Plots \n")

# Plot mesh
p1 <- plot_mesh(mesh, plot_settings) +
  ggtitle("Mesh")
p2 <- p1 +
  geom_point(aes(x_coord, y_coord), xy_coords[final_locations_indexes, ], color = "blue") +
  ggtitle("Mesh and final locations")
plot <- grid.arrange(p1, p2, ncol=2)
ggsave(paste(images.directory, "SlideseqCerebellum_mesh", ".jpg", sep = ""),
       plot = plot, width = 16, height = 8, dpi = 300)

# Plot locations
p1 <- plot_settings +
      geom_point(aes(x_coord, y_coord), xy_coords, color = "purple") + 
      geom_point(aes(x_coord, y_coord), xy_coords[final_locations_indexes, ], color = "grey") +
      geom_point(aes(x_coord, y_coord), xy_coords[islands.indexes, ], color = "red") +
      ggtitle("Original locations")
p2 <- plot_settings +
      geom_point(aes(x_coord, y_coord), xy_coords[final_locations_indexes,], color = "grey") +
      ggtitle("Final locations")
plot <- grid.arrange(p1, p2, ncol=2)
ggsave(paste(images.directory, "SlideseqCerebellum", ".jpg", sep = ""),
       plot = plot, width = 16, height = 8, dpi = 300)

# Export data
cat("- Exporting data \n")
save(mesh, locations,
     xy_coords, final_locations_indexes,
     file = paste(processed_data.directory, "SlideseqCerebellum_mesh", ".RData", sep = ""))


# ||||||||||||||||||||||||||
# SlideseqV2Hippocampus ----
# ||||||||||||||||||||||||||

rm(list = ls()[! ls() %in% c("plot_field", "plot_mesh", "plot_settings", "images.directory", "data.directory", "processed_data.directory", "generate_mesh")])

# Load data
cat("- Loading raw data \n")
list.files(data.directory)
list.files(paste(data.directory, "SlideseqV2Hippocampus/", sep = ""))
load(paste(data.directory, "SlideseqV2Hippocampus/Puck_200115_08_count_location.RData", sep = ""))
xy_coords <- data.frame(location)
colnames(xy_coords) <- c("x_coord", "y_coord")
xy_coords$ID = 1:dim(xy_coords)[1]

# "Islands" identification
cat("- Islands identification \n")
k <- 50
knn_result <- get.knn(xy_coords[,1:2], k)
xy_coords$avg_distance <- rowMeans(knn_result$nn.dist)
threshold <- 2 * mean(xy_coords$avg_distance)
islands.indexes <- xy_coords$avg_distance > threshold

# Mesh
cat("- Meshing \n")
indexes <- !islands.indexes
locations <- xy_coords[indexes, c("x_coord", "y_coord", "ID")]
mesh_gen <- generate_mesh(locations,
                          concavity = 2,
                          length_threshold = 200,
                          simplify_tol = 20,
                          maximum_area = 4000)
mesh <- mesh_gen$mesh
final_locations_indexes <- mesh_gen$final_locations_indexes

cat("- Plots \n")

# Plot mesh
p1 <- plot_mesh(mesh, plot_settings) +
      ggtitle("Mesh")
p2 <- p1 +
      geom_point(aes(x_coord, y_coord), xy_coords[final_locations_indexes, ], color = "blue") +
      ggtitle("Mesh and final locations")
plot <- grid.arrange(p1, p2, ncol=2)
ggsave(paste(images.directory, "SlideseqV2Hippocampus_mesh", ".jpg", sep = ""),
       plot = plot, width = 16, height = 8, dpi = 300)

# Plot locations
p1 <- plot_settings +
      geom_point(aes(x_coord, y_coord), xy_coords, color = "purple") + 
      geom_point(aes(x_coord, y_coord), xy_coords[final_locations_indexes, ], color = "grey") +
      geom_point(aes(x_coord, y_coord), xy_coords[islands.indexes, ], color = "red") +
      ggtitle("Original locations")
p2 <- plot_settings +
      geom_point(aes(x_coord, y_coord), xy_coords[final_locations_indexes,], color = "grey") +
      ggtitle("Final locations")
plot <- grid.arrange(p1, p2, ncol=2)
ggsave(paste(images.directory, "SlideseqV2Hippocampus", ".jpg", sep = ""),
       plot = plot, width = 16, height = 8, dpi = 300)

# Export data
cat("- Exporting data \n")
save(mesh, locations,
     xy_coords, final_locations_indexes,
     file = paste(processed_data.directory, "SlideseqV2Hippocampus_mesh", ".RData", sep = ""))


# ||||||||||||||||
# BreastTumor ----
# ||||||||||||||||

rm(list = ls()[! ls() %in% c("plot_field", "plot_mesh", "plot_settings", "images.directory", "data.directory", "processed_data.directory", "generate_mesh")])

# Load data
cat("- Loading raw data \n")
list.files(data.directory)
list.files(paste(data.directory, "BreastTumor/", sep = ""))
load(paste(data.directory, "BreastTumor/Tumor_data.RData", sep = ""))
xy_coords <- data.frame(location)
colnames(xy_coords) <- c("x_coord", "y_coord")
xy_coords$ID = 1:dim(xy_coords)[1]

# location matrix: n x 2, count matrix: g x n.
# here n is spot number, g is gene number.
# here the column names of sp_count and rownames of location should be matched
ST = CreateSpatialPCAObject(counts=rawcount, location=location,
                            project = "SpatialPCA",gene.type="spatial",
                            sparkversion="spark", gene.number=3000, 
                            customGenelist=NULL, min.loctions = 20, min.features=20)
data <- ST@normalized_expr

# "Islands" identification
cat("- Islands identification \n")
k <- 50
knn_result <- get.knn(xy_coords[,1:2], k)
xy_coords$avg_distance <- rowMeans(knn_result$nn.dist)
threshold <- 2 * mean(xy_coords$avg_distance)
islands.indexes <- xy_coords$avg_distance > threshold

# Mesh
cat("- Meshing \n")
indexes <- !islands.indexes
locations <- xy_coords[indexes, c("x_coord", "y_coord", "ID")]
mesh_gen <- generate_mesh(locations,
                          concavity = 1,
                          length_threshold = 0,
                          simplify_tol = 0,
                          maximum_area = 1)
mesh <- mesh_gen$mesh
final_locations_indexes <- mesh_gen$final_locations_indexes

cat("- Plots \n")

# Plot mesh
p1 <- plot_mesh(mesh, plot_settings) +
      ggtitle("Mesh")
p2 <- p1 +
      geom_point(aes(x_coord, y_coord), xy_coords[final_locations_indexes, ], color = "blue") +
      ggtitle("Mesh and final locations")
plot <- grid.arrange(p1, p2, ncol=2)
ggsave(paste(images.directory, "BreastTumor_mesh", ".jpg", sep = ""),
       plot = plot, width = 16, height = 8, dpi = 300)

# Plot locations
p1 <- plot_settings +
      geom_point(aes(x_coord, y_coord), xy_coords, color = "purple") + 
      geom_point(aes(x_coord, y_coord), xy_coords[final_locations_indexes, ], color = "grey") +
      geom_point(aes(x_coord, y_coord), xy_coords[islands.indexes, ], color = "red") +
      ggtitle("Original locations")
p2 <- plot_settings +
      geom_point(aes(x_coord, y_coord), xy_coords[final_locations_indexes,], color = "grey") +
      ggtitle("Final locations")
plot <- grid.arrange(p1, p2, ncol=2)
ggsave(paste(images.directory, "BreastTumor", ".jpg", sep = ""),
       plot = plot, width = 16, height = 8, dpi = 300)

# Export data
cat("- Exporting data \n")
save(data, file = paste(processed_data.directory, "BreastTumor_data", ".RData", sep = ""))
save(mesh, locations,
     xy_coords, final_locations_indexes,
     file = paste(processed_data.directory, "BreastTumor_mesh", ".RData", sep = ""))


# |||||||||||
# Vizgen ----
# |||||||||||

rm(list = ls()[! ls() %in% c("plot_field", "plot_mesh", "plot_settings", "images.directory", "data.directory", "processed_data.directory", "generate_mesh")])

# Load data
cat("- Loading raw data \n")
list.files(data.directory)
list.files(paste(data.directory, "Vizgen/", sep = ""))
load(paste(data.directory, "Vizgen/Vizgen_Merfish_count_location.RData", sep = ""))
xy_coords <- data.frame(location)
colnames(xy_coords) <- c("x_coord", "y_coord")

# Data pre-processing
Vizgen <- CreateSpatialPCAObject(counts = raw_matrix, location = location,
                                 project = "SpatialPCA", 
                                 gene.type = "spatial", sparkversion = "sparkx",
                                 gene.number = 3000, customGenelist = NULL,
                                 min.loctions = 50, min.features = 85) 
data <- Vizgen@normalized_expr
cellnames <- colnames(data)
xy_coords <- xy_coords[cellnames, ]
xy_coords$ID = 1:dim(xy_coords)[1]

# "Islands" identification
cat("- Islands identification \n")
k <- 10
knn_result <- get.knn(xy_coords[,1:2], k)
xy_coords$avg_distance <- rowMeans(knn_result$nn.dist)
threshold <- 1.7 * mean(xy_coords$avg_distance)
islands.indexes <- xy_coords$avg_distance > threshold

plot_settings +
  geom_point(aes(x_coord, y_coord), xy_coords, color = "grey") +
  geom_point(aes(x_coord, y_coord), xy_coords[islands.indexes, ], color = "red") +
  ggtitle("Original locations")

# Mesh
cat("- Meshing \n")
indexes <- !islands.indexes
locations <- xy_coords[indexes, c("x_coord", "y_coord", "ID")]
mesh_gen <- generate_mesh(locations,
                          concavity = 3,
                          length_threshold = 400,
                          simplify_tol = 20,
                          maximum_area = 5000)
mesh <- mesh_gen$mesh
final_locations_indexes <- mesh_gen$final_locations_indexes

cat("- Plots \n")

# Plot mesh
p1 <- plot_mesh(mesh, plot_settings) +
      ggtitle("Mesh")
p2 <- p1 +
      geom_point(aes(x_coord, y_coord), xy_coords[final_locations_indexes, ], color = "blue") +
      ggtitle("Mesh and final locations")
plot <- grid.arrange(p1, p2, ncol=2)
ggsave(paste(images.directory, "Vizgen_mesh", ".jpg", sep = ""),
       plot = plot, width = 16, height = 8, dpi = 300)

# Plot locations
p1 <- plot_settings +
      geom_point(aes(x_coord, y_coord), xy_coords, color = "purple") + 
      geom_point(aes(x_coord, y_coord), xy_coords[final_locations_indexes, ], color = "grey") +
      geom_point(aes(x_coord, y_coord), xy_coords[islands.indexes, ], color = "red") +
      ggtitle("Original locations")
p2 <- plot_settings +
      geom_point(aes(x_coord, y_coord), xy_coords[final_locations_indexes,], color = "grey") +
      ggtitle("Final locations")
plot <- grid.arrange(p1, p2, ncol=2)
ggsave(paste(images.directory, "Vizgen", ".jpg", sep = ""),
       plot = plot, width = 16, height = 8, dpi = 300)

# Export data
cat("- Exporting data \n")
locations <- xy_coords[final_locations_indexes, c("x_coord", "y_coord")]
save(data,
     file = paste(processed_data.directory, "Vizgen_data", ".RData", sep = ""))
save(mesh, locations,
     xy_coords, final_locations_indexes,
     file = paste(processed_data.directory, "Vizgen_mesh", ".RData", sep = ""))
