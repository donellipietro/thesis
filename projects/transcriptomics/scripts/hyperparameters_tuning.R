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
library(gridExtra)
library(mclust)


# ||||||||||||||
# Functions ----
# ||||||||||||||

load("../../utils/functions/plot_field.RData")


# |||||||||
# Data ----
# |||||||||

cat("\n# |||||||||||||||||||")
cat("\n# Importing data ----")
cat("\n# |||||||||||||||||||\n\n")

test.name <- "SlideseqCerebellum" # "DLPFC", "SlideseqV2Hippocampus", "BreastTumor", "SlideseqCerebellum", Vizgen"
TARGET <- FALSE
REFINE <- FALSE

load(paste("data/processed/", test.name, "_mesh.RData", sep = ""))
load(paste("data/processed/", test.name, "_data.RData", sep = ""))
#data <- read.csv(paste("data/processed/", test.name, "_data.csv", sep = ""), header = TRUE)[-1]

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

# Ground truth
# ||||||||||||

if(TARGET)
  target <- target


# |||||||||||||||||||||||||||
# Model hyper-parameters ----
# |||||||||||||||||||||||||||

hyperparameters_file <- paste("data/processed/", test.name, "_hyperparameters.RData", sep = "")
if (!file.exists(hyperparameters_file)){
  
  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    
  # fPCA
  lambda_vect <- 10^seq(-4, 0, by = 0.5)
  npc_vect <- c(5, 10, 15, 20)
  npc_max <- max(npc_vect)
  
  #  Clustering
  clustering.algorithm.name <- "Louvain" # "Walktrap", "Louvain"
  clusternum <- 8
  knearest_vect <- seq(100, 220, by = 10)
  
  # Plots
  pointsize <- 1.5
  
  save(lambda_vect, npc_vect, npc_max,
       clustering.algorithm.name, clusternum, knearest_vect,
       pointsize,
       file = hyperparameters_file)
  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

} else {
  load(hyperparameters_file)
}

# Names
lambda_vect_names <- paste("l", formatC(lambda_vect, format = "e", digits = 1), sep = "")
npc_vect_names <- paste("npc", npc_vect, sep = "")
knearest_vect_names <- paste("knn", knearest_vect, sep = "")

# Room for results
time <- rep(0, length(lambda_vect_names))
names(time) <- lambda_vect_names
ARI <- matrix(0, nrow = length(lambda_vect), length(npc_vect))
rownames(ARI) <- lambda_vect_names
colnames(ARI) <- npc_vect_names
ARI_refined <- matrix(0, nrow = length(lambda_vect), length(npc_vect))
rownames(ARI_refined) <- lambda_vect_names
colnames(ARI_refined) <- npc_vect_names
best_knn <- matrix(0, nrow = length(lambda_vect), length(npc_vect))
rownames(best_knn) <- lambda_vect_names
colnames(best_knn) <- npc_vect_names
best_knn_refined <- matrix(0, nrow = length(lambda_vect), length(npc_vect))
rownames(best_knn_refined) <- lambda_vect_names
colnames(best_knn_refined) <- npc_vect_names

cbp = c("#1F77B4", "#FF7F0E", "#2CA02C" ,"#D62728", "#9467BD" ,"#8C564B", "#E377C2",
        "#7F7F7F", "#BCBD22", "#17BECF", "#AEC7E8" ,"#FFBB78" ,"#98DF8A", "#FF9896",
        "#C5B0D5" ,"#C49C94", "#F7B6D2", "#C7C7C7" ,"#DBDB8D" ,"#9EDAE5", "#393B79",
        "#637939", "#8C6D31", "#843C39", "#7B4173" ,"#5254A3" ,"#8CA252", "#BD9E39",
        "#AD494A", "#A55194", "#6B6ECF", "#B5CF6B", "#E7BA52" ,"#D6616B", "#CE6DBD",
        "#9C9EDE", "#CEDB9C" ,"#E7CB94", "#E7969C", "#DE9ED6" ,"#3182BD", "#E6550D",
        "#31A354", "#756BB1" ,"#636363", "#6BAED6" ,"#FD8D3C" ,"#74C476", "#9E9AC8",
        "#969696", "#9ECAE1" ,"#FDAE6B", "#A1D99B" ,"#BCBDDC" ,"#BDBDBD", "#C6DBEF",
        "#FDD0A2" ,"#C7E9C0" ,"#DADAEB", "#D9D9D9")


# |||||||||||||
# Analysis ----
# |||||||||||||

cat("\n# |||||||||||||")
cat("\n# Analysis ----")
cat("\n# |||||||||||||\n\n")

# Create images directory
images.directory <- paste("images/", test.name, "/", sep = "")
if (!file.exists(images.directory)){
  dir.create(images.directory)
}
images.directory <- paste(images.directory, "hyperparameters_tuning", "/", sep = "")
if (!file.exists(images.directory)){
  dir.create(images.directory)
}

principal_components <- list()
clustering <- list()

for(l in 1:length(lambda_vect)){
  
  lambda <- lambda_vect[l]
  clustering[[lambda_vect_names[l]]] <- list()
  
  name <- paste("l", formatC(lambda, format = "e", digits = 1), sep = "")
  images_test.directory <- paste(images.directory, name, "/", sep = "")
  
  # Create images directory
  if (!file.exists(images_test.directory)){
    dir.create(images_test.directory)
  }
  
  cat(paste("\n\nTest ", l, ": lambda = ", name, "\n\n", sep = ""))
  
  # |||||||||||||||
  # fPCA model ----
  # |||||||||||||||
  
  cat("\n# fPCA ----")
  cat("\n# |||||||||\n\n")
  
  # Set model
  pde <- new(Laplacian_2D_Order1, mesh_data)
  quadrature_nodes <- pde$get_quadrature_nodes()
  f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
  pde$set_forcing_term(as.matrix(f))
  
  # Define and init model
  model <- new(FPCA_CS_Laplacian_2D_GeoStatLocations, pde)
  model$set_locations(locations)
  model$set_npc(npc_max)
  model$set_lambda_s(lambda)
  model$set_mass_lumping(TRUE)
  model$set_coefficients_position(2)
  model$init_regularization()
  
  # Set data
  model$set_observations(data)
  
  # Fit the model
  tic()
  model$solve()
  elapsed <- toc()
  
  time[l] <- elapsed$toc - elapsed$tic
  
  # Results
  load <- matrix(0, nrow = npc_max, ncol = K)
  for(i in 1:npc_max){
    load[i, ] <- temp <-  model$loadings()[, i]
    
    jpeg(paste(images_test.directory, "load", i, ".jpg", sep = ""), width = 1200, height = 1200, units = "px", quality = 100, pointsize = 37)
    plot_field(nodes, locations, load[i, ], range(load[i, ]), paste("load", i, sep = ""), TRUE)
    dev.off()
  }
  
  principal_components[[name]] <- load
  
  cat("\n")
  
  
  cat("\n# Clustering ----")
  cat("\n# ||||||||||||||||\n\n")
  
  
  for(h in 1:length(npc_vect)){
    
    npc <- npc_vect[h]
    
    cat(paste("Number of components: ", npc, "\n", sep = ""))
    
    clusterlabel <- list()
    clusterlabel_refine <- list()
    
    
    for(i in 1:length(knearest_vect)){
      
      cat(paste("- knn =", knearest_vect[i], "\n"))
      
      if(clustering.algorithm.name == "Walktrap"){ 
        clusterlabel[[knearest_vect_names[i]]] = walktrap_clustering(clusternum = clusternum,
                                                latent_dat = load[1:npc, ],
                                                knearest = knearest_vect[i])
        if(REFINE){
          clusterlabel_refine[[knearest_vect_names[i]]] = refine_cluster_10x(clusterlabels=clusterlabel[[knearest_vect_names[i]]],
                                                        locations, shape="hexagon")
        }
      }
      
      if(clustering.algorithm.name == "Louvain"){
        clusterlabel[[knearest_vect_names[i]]] = louvain_clustering(clusternum = clusternum,
                                               latent_dat = load[1:npc, ],
                                               knearest = knearest_vect[i])
        if(REFINE){
          clusterlabel_refine[[knearest_vect_names[i]]] = refine_cluster_10x(clusterlabels=clusterlabel[[knearest_vect_names[i]]],
                                                        locations, shape="hexagon")
        }
      }
      
      
    }
    
    clustering[[lambda_vect_names[l]]][[npc_vect_names[h]]] <- list()
    
    # Clustering
    clustering[[lambda_vect_names[l]]][[npc_vect_names[h]]]$normal <- clusterlabel
    ari <- rep(0, length(knearest_vect))
    plots <- list()
    for(i in 1:length(knearest_vect)){
      title <- paste("knn =", knearest_vect[i])
      if(TARGET){
        ari[i] <- adjustedRandIndex(target, clusterlabel[[knearest_vect_names[i]]])
        title <- paste(title, "-", "ARI =", round(ari[i], digit = 3))
      }
      plots[[i]] <- plot_cluster(legend="none",location=locations,
                                 clusterlabel[[knearest_vect_names[i]]],
                                 pointsize=pointsize,text_size=20,
                                 title_in=title,
                                 color_in=cbp)

    }
    n <- length(plots)
    nCol <- floor(sqrt(n))
    plot <- do.call("grid.arrange", c(plots, ncol=nCol))
    ggsave(paste(images_test.directory, "clustering_npc", npc, ".jpg", sep = ""),
           plot = plot, width = 16, height = 18, dpi = 300)
    
    if(TARGET){
      ARI[l, h] <- max(ari)
      best_knn[l, h] <- knearest_vect[which(ari == max(ari))]
    }
    
    # Clustering refined
    if(REFINE){
      clustering[[lambda_vect_names[l]]][[npc_vect_names[h]]]$refined <- clusterlabel_refine
      ari_refined <- rep(0, length(knearest_vect))
      for(i in 1:length(knearest_vect)){
        title <- paste("knn =",knearest_vect[i])
        if(TARGET){
          ari_refined[i] <- adjustedRandIndex(target, clusterlabel_refine[[knearest_vect_names[i]]])[1]
          title <- paste(title, "-", "ARI =", round(ari_refined[i], digit = 3))
        }
        plots[[i]] <- plot_cluster(legend="none",location=locations,
                                   clusterlabel_refine[[knearest_vect_names[i]]],
                                   pointsize=pointsize,text_size=20,
                                   title_in=title,
                                   color_in=cbp)
  
      }
      n <- length(plots)
      nCol <- floor(sqrt(n))
      plot <- do.call("grid.arrange", c(plots, ncol=nCol))
      ggsave(paste(images_test.directory, "clustering_refined_npc", npc, ".jpg", sep = ""),
             plot = plot, width = 16, height = 18, dpi = 300)
    
      if(TARGET){
        ARI_refined[l, h] <- max(ari_refined)
        best_knn_refined[l, h] <- knearest_vect[which(ari_refined == max(ari_refined))][1]
      }
    }
    
    cat("\n")
    
  }
  
  rm(model)

}

# cat("\nTime:\n")
# cat(time)
# 
# cat("\nARI:\n")
# cat(ARI)
# 
# cat("\nBest KNN:\n")
# cat(best_knn)
# 
# cat("\nARI:\n")
# cat()
# 
# cat("\nBest KNN:\n")
# cat(best_knn_refined)

results.directory <- paste("results/", "hyperparameters_tuning", "/", sep = "")
# Create results directory
if (!file.exists(results.directory)){
  dir.create(results.directory)
}

save(lambda_vect, lambda_vect_names, npc_vect, npc_vect_names,
     knearest_vect, knearest_vect_names, clusternum, clustering.algorithm.name,
     principal_components,
     clustering,
     target,
     locations,
     time,
     ARI, best_knn,
     ARI_refined, best_knn_refined,
     file = paste(results.directory, "results_", test.name, "_", format(Sys.time(), "_%Y%m%d_%H%M%S"), ".RData", sep = ""))

