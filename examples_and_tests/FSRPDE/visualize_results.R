# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Utility to visualize results of FSRPDE testing %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()


# ||||||||||||||||||||||||||
# Environment variables ----
# ||||||||||||||||||||||||||

# Where to get the data
tests_dir = "../../fdaPDE/test/data/models/FSRPDE/"


# ||||||||||||||
# Functions ----
# ||||||||||||||

plot_results <- function(path, at_nodes, plot_locations) {
  
  par(mfrow = c(2,2))
  
  # Original data
  X_clean <- read.csv(paste(path,"X_clean.csv", sep = ''), header = TRUE)[,2]
  image(matrix(as.numeric(X_clean), 60, 60), main = "X_clean")
  
  # Sample of noisy data
  if(at_nodes){
    X <- read.csv(paste(path,"X.csv", sep = ''), header = TRUE)[,2:3601]
    image(matrix(as.numeric(X[1,]), 60, 60), main = "X sample")
  }
  
  # Estimated field
  f <- read.csv(paste(path,"results/","f.csv", sep = ''), header = FALSE)
  image(matrix(as.numeric(f[1,]), 60, 60), main = "f")
  
  # Residual
  image(matrix(as.numeric(f[1,]-X_clean), 60, 60), main = "error")
  
  # Locations
  if(plot_locations){
    locs <- read.csv(paste(path,"locations.csv", sep = ''), header = TRUE)[,2:3]
    plot(locs, pch = 8, col = "red", main = "Locations")
  }
  
  print(max(abs(X_clean-f)))
  print(min(abs(X_clean-f)))
}


# ||||||||||||||||||||||||||
# Results visualization ----
# ||||||||||||||||||||||||||

path_images <- "images/"
if (!file.exists(path_images)){
  dir.create(path_images)
}

## Test 1 ----
## |||||||||||

# domain:       unit square [0,1] x [0,1]
# sampling:     locations = nodes
# penalization: simple laplacian
# covariates:   no
# BC:           no
# order FE:     1

test_dir = paste(tests_dir, "2D_test1/", sep = '')
# jpeg(file=paste(path_images, "test1.jpg", sep = ''))
plot_results(test_dir, at_nodes = TRUE, plot_locations = FALSE)
# dev.off();

## Test 2 ----
## |||||||||||

# domain:       [0,1] x [0,1]
# sampling:     locations != nodes, #nodes = #locations
# penalization: simple laplacian
# covariates:   no
# BC:           no
# order FE:     1

test_dir = paste(tests_dir, "2D_test2/", sep = '') 
jpeg(file=paste(path_images, "test2.jpg", sep = ''))
plot_results(test_dir, at_nodes = FALSE, plot_locations = FALSE)
dev.off();

## Test 3 ----
## |||||||||||

# domain:       [0,1] x [0,1]
# sampling:     locations != nodes, #locations < #nodes (10%)
# penalization: simple laplacian
# covariates:   no
# BC:           no
# order FE:     1

test_dir = paste(tests_dir, "2D_test3/", sep = '') 
# jpeg(file=paste(path_images, "test3.jpg", sep = ''))
plot_results(test_dir, at_nodes = FALSE, plot_locations = TRUE)
# dev.off();


## Test 4 ----
## |||||||||||

# domain:       [0,1] x [0,1]
# sampling:     locations != nodes, #locations << #nodes (equispaced)
# penalization: simple laplacian
# covariates:   no
# BC:           no
# order FE:     1

test_dir = paste(tests_dir, "2D_test4/", sep = '') 
# jpeg(file=paste(path_images, "test4.jpg", sep = ''))
plot_results(test_dir, at_nodes = FALSE, plot_locations = TRUE)
# dev.off();

