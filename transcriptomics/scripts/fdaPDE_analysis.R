# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% fdaPDE anlysis on DLPFC data %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()
def.par = par()


# ||||||||||||||
# Libraries ----
# ||||||||||||||

library(pracma)
library(plot3D)

library(fdaPDE2)
library(tictoc)


# ||||||||||||||
# Functions ----
# ||||||||||||||

load("functions/plot_field.RData")


# |||||||||
# Data ----
# |||||||||

load("../data/data_fdaPDE.RData")


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


# Field
# |||||

# Locations 
locations = as.matrix(X_new)
n_spatial_locations <- nrow(locations)
n_locations = n_spatial_locations

# Field
N <- dim(V_new)[1]
K <- dim(V_new)[2]
data = as.matrix(V_new)# [sample.int(N, 2000, replace =  TRUE),]
plot_field(nodes, locations, data[1,], range(data), "Gene 1", TRUE)


# |||||||||||||||
# fPCA model ----
# |||||||||||||||

# Set model
pde <- new(Laplacian_2D_Order1, mesh_data)
quadrature_nodes <- pde$get_quadrature_nodes()
f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
pde$set_forcing_term(as.matrix(f))

# Define and init model
model <- new(FPCA_CS_Laplacian_2D_GeoStatLocations, pde)
model$set_locations(locations)
npc <- 20
model$set_npc(npc)
lambda_s <- 100 # 10^seq(-4.0, -3.0, by = 0.2)
model$set_lambda_s(lambda_s)
model$set_mass_lumping(1)
model$init_regularization()


## extract principal components
model$set_observations(data)
tic()
model$solve()
elapsed <- toc()

load <- matrix(0, nrow = npc, ncol = K)
for(i in 1:npc){
  load[i, ] <- model$loadings()[, i]
  jpeg(paste("../images/laod", i, ".jpg", sep = ""), width = 1500, height = 1500, units = "px", quality = 100, pointsize = 37)
  plot_field(nodes, locations, load[i, ], range(load[i, ]), paste("load", i, sep = ""), TRUE)
  dev.off()
}


save(load, file = "../results/results_3000.RData")






load("../results/results_3000.RData")

dist <- dist(t(load[c(1, 2, 7, 11, 15) , ]), method='canberra')
hier <- hclust(dist, "average")

plot(hier)
rect.hclust(hier, k=10)



cluster = cutree(hier, k=10)
plot_field(nodes, locations, cluster, range(cluster), "Clustering", TRUE)




if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
if (!requireNamespace("spatialLIBD", quietly = TRUE)) {BiocManager::install("spatialLIBD")}

library(SPARK)
library(Seurat)
library(peakRAM)
library(SpatialPCA)
library(ggplot2)

# library(devtools)
# install_github("shangll123/SpatialPCA")

clusternum = 7
knn = 70
SPCA_clusterlabels <- walktrap_clustering(clusternum = clusternum, latent_dat=load, knearest = knn)
SPCA_clusterlabels_refined <- refine_cluster_10x(clusterlabels = SPCA_clusterlabels, location=DLPFC_SPCA@location, shape="hexagon")


