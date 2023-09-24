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
data.directory <- "data/number_statistical_units/"
results.directory <- "results/number_statistical_units/"
images.directory <- "images/number_statistical_units/"

# Mesh
mesh.name_vect <- c("unit_square_medium")
mesh.options <- list("unit_square_medium" = list("mesh_finer.path" = "",
                                                 "mesh.area_refine" = 0,
                                                 "generate_X" = "../../utils/functions/tests_unit_square.RData"))

# Number of statistical units
N_vect <- c(10, 50, 100, 200, 400)
N.batches <- 10
N.generated <- max(N_vect*N.batches)

# Type of field
X.index_vect <- c(2)

# Noise level
NSR_vect <- c(1/4, 1/3, 1/2, 1, 5/4, 6/4)^2

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
  X.index = paste("xi", c(1), sep = ""),
  NSR_vect = paste("nl", 1:length(NSR_vect), sep = "")
)

test.options <- generate_data_combinations(N.generated, 
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

lambda_s <- c(1e-8, 1e-6, 1e-4)

errors <- list()
for(l in 1:length(lambda_s)){
  errors[[l]] <- list()
  for(i in 1:dim(test.options)[1])
    errors[[l]][[test.options_names$NSR_vect[i]]] <- matrix(0, nrow = N.batches, ncol = length(N_vect))
}

for(i in 1:dim(test.options)[1]){
  
  cat(paste("\n", test.options$test.name[i], sep= ""))
  
  # Load data
  load(test.options$data.path[i])
  
  # Save 1 sample for data visualization
  # data.examples[[i]] <- data[1,]
  # data_clean.examples[[i]] <- data_clean[1,]
  # nodes.examples[[i]] <- nodes
  # locations.examples[[i]] <- locations
  
  # Set model
  pde <- new(Laplacian_2D_Order1, mesh_data)
  
  # Set zero forcing term
  quadrature_nodes <- pde$get_quadrature_nodes()
  f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
  pde$set_forcing_term(as.matrix(f))
  
  for(l in 1:length(lambda_s)){
  
    # Define and init model
    model <- new(FRPDE_Laplacian_2D_GeoStatLocations, pde)
    model$set_locations(locations)
    model$set_lambda_s(lambda_s[l])
    
    for(n in 1:length(N_vect)){
      
      N <- N_vect[n]
    
      # Generate a grouping variable for split
      group_var <- rep(1:N, each = N, length.out =  N*N.batches)
      
      # Split the vector into N subvectors without repetition
      set.seed(0)
      subvectors <- split(sample(1:(N*N.batches)), group_var)
      
      for(j in 1:N.batches){
        
        # Set observations
        model$set_observations(data[subvectors[[j]],])
        
        # Solve
        model$solve()
        
        # Results
        errors[[l]][[test.options_names$NSR_vect[i]]][j,n] <- norm(data_clean[1,] - model$fitted())/sqrt(S)
        
      }
    
    }

    rm(model)
  
  }
  
}

save(N_vect, X.index_vect, NSR_vect, test.options,
     errors,
     # data.examples, data_clean.examples,
     # nodes.examples, locations.examples,
     # results,
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



merge_tests <- function(test1, test2, n_tests, n_test_options_1, n_test_options_2){
  n_batches <- min(dim(test1)[1], dim(test2)[1])
  n_test_options <- n_test_options_1 + n_test_options_2
  test_merged <- matrix(0, nrow = n_batches, ncol = 0)
  for(i in 1:n_tests){
    test_merged <- cbind(test_merged, test1[1:n_batches, n_test_options_1*(i-1) + 1:n_test_options_1])
    test_merged <- cbind(test_merged, test2[1:n_batches, n_test_options_2*(i-1) + 1:n_test_options_2])
  }
  
  return(test_merged)
}

plot_comparison <- function(comparison, mains, title, ylab, xticks, models, colors) {
  
  par(mfrow = c(length(lambda_s),1), mar = c(4,4,5,1), oma = c(0,0,3,0))
  
  for(i in 1:length(lambda_s)){
    
    n_test <- dim(comparison[[i]][[models[1]]])[2]
    merged <- comparison[[i]][[models[1]]]
    n_test_options_1 <- 1
    for(m in models[-1]){
      merged <- merge_tests(merged, comparison[[i]][[m]], n_test, n_test_options_1, 1)
      n_test_options_1 <- n_test_options_1 + 1
    }
    
    ylim = c(0, max(merged))
    t_temp <- NULL
    x_temp <- NULL
    n_models <- length(models)
    iter <- 1
    for(m in models){
      t_temp <- rbind(t_temp, 0:(n_test-1)*n_models+iter)
      x_temp <- rbind(x_temp, colMeans(comparison[[i]][[m]]))
      iter <- iter + 1
    }
    matplot(t(t_temp), 
            t(x_temp),
            main = mains[i],
            col = colors, type = "o", lwd = 2, pch = 19, lty = 1,
            xlab = "N", xaxt = "n",
            ylab = ylab, ylim = ylim)
    boxplot(merged, col = colors, add = TRUE, xaxt = "n")
    
    axis(1, at=(1:(n_test))*(n_models)-(n_models-1)/2, labels=xticks)
    grid(nx = NA, ny = NULL)
    abline(v = 0:(n_test-1)*n_models+0.5, col = "lightgrey", lty = 3)
    if(i == 1)
      legend("topright", paste("NSR =", round(NSR_vect, digits = 2)),
             col = colors, pch = 15)
    
  }
  
  title(title, outer = TRUE, line = 1, cex.main = 1.5)
  
}


# Pictures settings
width = 2000
height_per_row = 800
pointsize = 40


# Execution time
jpeg(paste(images.directory, "time", ".jpg", sep = ""),
     width = width, height = height_per_row*length(lambda_s), units = "px", quality = 100, pointsize = pointsize,
     type = "windows")
plot_comparison(errors,
                paste("lambda =", formatC(lambda_s, format = "e", digits = 1)),
                "Differential Regularization effect",
                "RMSE",
                N_vect, test.options_names$NSR_vect, colors)
dev.off()










