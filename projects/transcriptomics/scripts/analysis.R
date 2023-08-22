# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Spatial Transcriptomics Project: Hyper-parameters tuning %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()
def.par = par()

cat("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
cat("\n%% Spatial Transcriptomics Project: Hyper-parameters tuning %%")
cat("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")


# ||||||||||||||
# Libraries ----
# ||||||||||||||

library(pracma)
library(plot3D)

library(fdaPDE2)
library(tictoc)

library(SpatialPCA)
library(ggplot2)
library(mclust)


# ||||||||||||||
# Functions ----
# ||||||||||||||

load("../../utils/functions/plot_field.RData")

test.name <- "Vizgen"
lambda <- 1e-2
npc <- 20

images_test.directory <- paste("images/", test.name, "/", sep = "")

# Create images directory
if (!file.exists(images_test.directory)){
  dir.create(images_test.directory)
}

images_test.directory <- paste(images_test.directory, "npc", npc, "_", "l", lambda, "/", sep = "")

# Create images directory
if (!file.exists(images_test.directory)){
  dir.create(images_test.directory)
}


# |||||||||
# Data ----
# |||||||||

cat("\n# |||||||||||||||||||")
cat("\n# Importing data ----")
cat("\n# |||||||||||||||||||\n\n")

load(paste("data/processed/", test.name, "_mesh.RData", sep = ""))
load(paste("data/processed/", test.name, "_data.RData", sep = ""))

# Domain
# ||||||

nodes = mesh$nodes
mesh_data <- list(
  "nodes"    = mesh$nodes,
  "edges"    = mesh$edges,
  "elements" = mesh$triangles,
  "neigh"    = mesh$neighbors,
  "boundary" = mesh$nodesmarkers
)

# Locations
# |||||||||

locations <- as.matrix(locations)
n_locations <- nrow(locations)

# Field
# |||||

data = as.matrix(data[, final_locations_indexes])
N <- dim(data)[1]
K <- dim(data)[2]
# plot_field(nodes, locations, data[170,], range(data), "Gene 1", TRUE)


# |||||||||||||
# Analysis ----
# |||||||||||||

cat("\n# |||||||||||||")
cat("\n# Analysis ----")
cat("\n# |||||||||||||\n\n")


# fPCA model
# ||||||||||

cat("\n# fPCA ----")
cat("\n# |||||||||\n\n")

dim(data)
dim(locations)

# Set model
pde <- new(Laplacian_2D_Order1, mesh_data)
quadrature_nodes <- pde$get_quadrature_nodes()
f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
pde$set_forcing_term(as.matrix(f))

# Define and init model
model <- new(FPCA_CS_Laplacian_2D_GeoStatLocations, pde)
model$set_locations(locations)
model$set_npc(npc)
model$set_lambda_s(lambda)
model$set_mass_lumping(TRUE)
model$init_regularization()

# Set data
model$set_observations(data)

# Fit the model
tic()
model$solve()
elapsed <- toc()
time <- elapsed$toc - elapsed$tic

# Results
load <- matrix(0, nrow = npc, ncol = K)
scores <- matrix(0, nrow = N, ncol = npc)
weight <- c()
for(i in 1:npc){
  scores[, i] <- model$scores()[, i]
  weight[i] <- norm(scores[, i], "2")
  load[i, ] <- weight[i]*model$loadings()[, i]
  
  jpeg(paste(images_test.directory, "load", i, ".jpg", sep = ""), width = 1500, height = 1500, units = "px", quality = 100, pointsize = 37)
  plot_field(nodes, locations, load[i, ], range(load[i, ]), paste("load", i, sep = ""), TRUE)
  dev.off()
}



cat("\n# Clustering ----")
cat("\n# ||||||||||||||||\n\n")







clusterlabel= louvain_clustering(clusternum=30,latent_dat=load, knearest=200 )


# set color
D3=c("#1F77B4", "#FF7F0E", "#2CA02C" ,"#D62728", "#9467BD" ,"#8C564B", "#E377C2",
     "#7F7F7F", "#BCBD22", "#17BECF", "#AEC7E8" ,"#FFBB78" ,"#98DF8A", "#FF9896",
     "#C5B0D5" ,"#C49C94", "#F7B6D2", "#C7C7C7" ,"#DBDB8D" ,"#9EDAE5", "#393B79",
     "#637939", "#8C6D31", "#843C39", "#7B4173" ,"#5254A3" ,"#8CA252", "#BD9E39",
     "#AD494A", "#A55194", "#6B6ECF", "#B5CF6B", "#E7BA52" ,"#D6616B", "#CE6DBD",
     "#9C9EDE", "#CEDB9C" ,"#E7CB94", "#E7969C", "#DE9ED6" ,"#3182BD", "#E6550D",
     "#31A354", "#756BB1" ,"#636363", "#6BAED6" ,"#FD8D3C" ,"#74C476", "#9E9AC8",
     "#969696", "#9ECAE1" ,"#FDAE6B", "#A1D99B" ,"#BCBDDC" ,"#BDBDBD", "#C6DBEF",
     "#FDD0A2" ,"#C7E9C0" ,"#DADAEB", "#D9D9D9")

# plot_field(locations, locations, as.numeric(clusterlabel), range(as.numeric(clusterlabel)),
#            "prova", TRUE)

plot <- plot_cluster(legend="none",location=locations,
             clusterlabel,
             pointsize=1,text_size=20,
             title_in=paste0("fPCA cluster, lambda = ",lambda, " npc = ", npc, sep = ""),
             color_in=D3)


ggsave(paste(images_test.directory, "clustering.jpg", sep = ""),
       plot = plot, width = 16, height = 16, dpi = 300)

# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
#   
# clusternum <- 7
# knearest_vect <- seq(25, 225, by = 25)
# 
# clusterlabel <- list()
# clusterlabel_refine <- list()
#   
# 
# for(i in 1:length(knearest_vect)){
#   
#   cat(paste("- knn =", knearest_vect[i], "\n"))
#   
#   
#   clusterlabel[[i]] = walktrap_clustering(clusternum = clusternum,
#                                           latent_dat = load[1:npc, ],
#                                           knearest = knearest_vect[i]) 
#   
#   clusterlabel_refine[[i]] = refine_cluster_10x(clusterlabels=clusterlabel[[i]],
#                                                 locations, shape="hexagon")
#   
#   
# }
# 
# # Clustering
# jpeg(paste(images_test.directory, "clustering_npc", npc, ".jpg", sep = ""),
#      width = 4500, height = 4500, units = "px", quality = 100, pointsize = 70,
#      type = "windows")
# par(mfrow = c(3, 3), mar = c(0,0,2,0), oma = c(0,0,0,0))
# for(i in 1:length(knearest_vect)){
#   plot_field(locations, locations, as.numeric(clusterlabel[[i]]), c(1,clusternum),
#              paste("knn =",knearest_vect[i]), FALSE)
# }
# dev.off()
# 
# # Clustering refined
# jpeg(paste(images_test.directory, "clustering_refined_npc", npc, ".jpg", sep = ""),
#      width = 4500, height = 4500, units = "px", quality = 100, pointsize = 70,
#      type = "windows")
# par(mfrow = c(3, 3), mar = c(0,0,2,0), oma = c(0,0,0,0))
# for(i in 1:length(knearest_vect)){
#   plot_field(locations, locations, as.numeric(clusterlabel_refine[[i]]), c(1,clusternum),
#              paste("knn =",knearest_vect[i]), FALSE)
# }
# dev.off()
# 
