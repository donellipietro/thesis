# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Test fSRPDE: noise suppression %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()
def.par = par()

cat("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
cat("\n%% Test fSRPDE: noise suppression %%")
cat("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n")


# ||||||||||||||
# Libraries ----
# ||||||||||||||

library(fdaPDE)
library(fdaPDE2)
library(pracma)
library(plot3D)
library(RColorBrewer)


# ||||||||||||||
# Functions ----
# ||||||||||||||

load("../../utils/functions/plot_field.RData")
load("../../utils/functions/import_fdaPDE_mesh.RData")
load("../../utils/functions/generate_data_combinations.RData")
load("scripts/functions/generate_data.RData")


# ||||||||||||||||||||
# Test parameters ----
# ||||||||||||||||||||

# Data directory
mesh.directory <- "../../fdaPDE/test/data/mesh/"
data.directory <- "data/noise/"
results.directory <- "results/noise/"
images.directory <- "images/noise/"

# Mesh
mesh.name_vect <- c("unit_square_medium", "c_shaped")
mesh.options <- list("unit_square_medium" = list("mesh_finer.path" = "",
                                                 "mesh.area_refine" = 0,
                                                 "generate_X" = "../../utils/functions/tests_unit_square.RData"),
                     "c_shaped" = list("mesh_finer.path" = "",
                                       "mesh.area_refine" = 0.004,
                                       "generate_X" = "../../utils/functions/tests_c_shaped.RData"))

# Number of statistical units
N <- 50
N.batches <- 10
N.generated <- N*N.batches

# Type of field
X.index_vect <- c(1,2)

# Noise level
NSR_vect <- c(1/4, 1/3, 1/2, 1)^2

# To force the generation of the field
FORCE_GENERATION = FALSE


# |||||||||
# Data ----
# |||||||||

cat("\n")
cat("\n# ||||||||||||||||||||")
cat("\n# Generating data ----")
cat("\n# ||||||||||||||||||||\n")

# Create test directory
if (!file.exists(data.directory)){
  dir.create(data.directory)
}

# List of test options
test.options <- list(
  mesh.name = mesh.name_vect,
  X.index = X.index_vect,
  NSR = NSR_vect
)

# List of test options names
test.options_names <- list(
  mesh.name_vect = mesh.name_vect,
  X.index = paste("xi", c(1,2), sep = ""),
  NSR_vect = paste("nl", 1:length(NSR_vect), sep = "")
)

test.options <- generate_data_combinations(N, 
                                           test.options, test.options_names, 
                                           mesh.directory, mesh.options,
                                           data.directory,
                                           generate_X,
                                           FORCE_GENERATION)

cat("\n")

# |||||||||
# Test ----
# |||||||||

cat("\n")
cat("\n# |||||||||")
cat("\n# Test ----")
cat("\n# |||||||||\n")

# Create results directory
if (!file.exists(results.directory)){
  dir.create(results.directory)
}

lambda_s <- c(1e-4, 1e-6)

errors <- matrix(0, nrow = N.batches, ncol = dim(test.options)[1])
data.examples <- list()
data_clean.examples <- list()
nodes.examples <- list()
locations.examples <- list()
results <- list()

for(i in 1:dim(test.options)[1]){
  
  cat(paste("\n", test.options$test.name[i], sep= ""))
  
  # Load data
  load(test.options$data.path[i])
  
  # Save 1 sample for data visualization
  data.examples[[i]] <- data[1,]
  data_clean.examples[[i]] <- data_clean[1,]
  nodes.examples[[i]] <- nodes
  locations.examples[[i]] <- locations
  
  # Set model
  pde <- new(Laplacian_2D_Order1, mesh_data)
  
  # Set zero forcing term
  quadrature_nodes <- pde$get_quadrature_nodes()
  f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
  pde$set_forcing_term(as.matrix(f))
  
  # Define and init model
  
  model <- new(FSRPDE_Laplacian_2D_GeoStatLocations, pde)
  model$set_locations(locations)
  model$set_lambda_s(lambda_s[as.numeric(test.options$X.index[i])])
  
  # Generate a grouping variable for split
  group_var <- rep(1:N, each = N %/% N.batches, length.out =  N)
  
  # Split the vector into N subvectors without repetition
  set.seed(0)
  subvectors <- split(sample(1:N), group_var)
  
  for(j in 1:N.batches){
  
    # Set observations
    model$set_observations(data[subvectors[[j]],])
    
    # Solve
    model$solve()
    
    # Results
    errors[j,i] <- norm(data_clean[1,] - model$fitted())/S
  
  }
  
  results[[i]] <- model$fitted()
  
  
  rm(model)
  
}

save(X.index_vect, NSR_vect, test.options,
     errors,
     data.examples, data_clean.examples,
     nodes.examples, locations.examples,
     results,
     file = paste(results.directory,"results", format(Sys.time(), "_%Y%m%d_%H%M%S"), ".RData", sep = ""))

cat("\n")

# ||||||||||||
# Results ----
# ||||||||||||

cat("\n")
cat("\n# ||||||||||||||||||||||")
cat("\n# Exporting results ----")
cat("\n# ||||||||||||||||||||||\n")


# Create images directory
if (!file.exists(images.directory)){
  dir.create(images.directory)
}

# list.files("../results/")
load(paste(results.directory, tail(list.files(results.directory), n = 1), sep = ""))

# Choose a color palette from RColorBrewer
palette_name <- "Spectral"

# Generate N different colors from the chosen palette
colors <- brewer.pal(length(NSR_vect), palette_name)

# Clean vs Noisy vs Denoised
cat("\nComparison")
for(m in 1:length(mesh.name_vect)){
  for(i in 1:length(X.index_vect)){
    
    selected <- which(test.options$X.index == X.index_vect[i] & test.options$mesh.name == mesh.name_vect[m])
    
    for(j in 1:length(NSR_vect)){
      
      jpeg(paste(images.directory, "comparison_",test.options$test.name[selected[j]], ".jpg", sep = ""),
           width = 2000, height = 2000, units = "px", quality = 100, pointsize = 37,
           type = "windows")
      
      par(mfrow = c(2, 2), mar = c(1,1,1,1), oma = c(0,0,1,0))
      
      range <- range(data.examples[[selected[j]]])
      
      plot_field(nodes.examples[[selected[j]]],
                 locations.examples[[selected[j]]],
                 data_clean.examples[[selected[j]]],
                 range,
                 "Clean data", FALSE)
      
      plot_field(nodes.examples[[selected[j]]],
                 locations.examples[[selected[j]]],
                 data.examples[[selected[j]]],
                 range,
                 "Noisy data", FALSE)
      
      plot_field(nodes.examples[[selected[j]]],
                 locations.examples[[selected[j]]],
                 results[[selected[j]]],
                 range,
                 "Denoised data", FALSE)
      
      plot_field(nodes.examples[[selected[j]]],
                 locations.examples[[selected[j]]],
                 data_clean.examples[[selected[j]]] - results[[selected[j]]],
                 range,
                 "Error", FALSE)
      
      title(paste("NSR =", round(NSR_vect[j], digits = 2)), outer = TRUE, line = 2, cex.main = 2)
      
      dev.off()
      
    }
    
  }
}

# MSE boxplots
cat("\nMSE boxplots")
for(m in 1:length(mesh.name_vect)){
  for(i in 1:length(X.index_vect)){
    
    jpeg(paste(images.directory, "boxplot_", mesh.name_vect[m], "_xi", X.index_vect[i], ".jpg", sep = ""), width = 2000, height = 1500, units = "px", quality = 100, pointsize = 37)
    
    selected <- which(test.options$X.index == X.index_vect[i] & test.options$mesh.name == mesh.name_vect[m])
      
    par(mfrow = c(1,1), mar = def.par$mar, omi = def.par$omi)
    boxplot(errors[, selected],
            xlab = "Noise level", ylab = "MSE",
            col = colors,
            main = paste("Mesh:", mesh.name_vect[m], "| X.index =", X.index_vect[i]))
    grid()
    legend("topleft", paste("NSR =", round(NSR_vect, digits = 2)),
           col = colors, pch = 15)
    
    dev.off()
    
  }
}


cat("\n\n\n")
