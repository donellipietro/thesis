# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Utility to generate data for FSRPDE testing %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()

# ||||||||||||||||||||||||||
# Environment variables ----
# ||||||||||||||||||||||||||

# Where to save the data
tests_dir <- "../../fdaPDE/test/data/models/FSRPDE/" # "data"


# ||||||||||||||
# Functions ----
# ||||||||||||||

generate_2d_data <- function(x, y, num_samples = 100) {
  
  # Nodes:
  nodes <- expand.grid(x = x, y = y)
  
  # Subgrid
  x_sub <- seq(0.1, 0.9, length.out = 9)
  y_sub <- seq(0.1, 0.9, length.out = 9)
  nodes_sub <- expand.grid(x = x_sub, y = y_sub)
  
  # Generate X:
  X_clean_nodes = NULL
  X_nodes = NULL
  
  X_clean_locations = NULL
  X_locations = NULL
  
  X_clean_locations_sub = NULL
  X_locations_sub = NULL
  
  locations = NULL
  locations_sub = NULL
  
  set.seed(0)
  
  ds = 1/length(x)/2
  Sigma = (ds/3)^2*matrix(c(1,0,0,1), 2, 2)
  dx = MASS::mvrnorm(dim(nodes)[1], mu = c(0,0), Sigma = Sigma);
  locations = nodes + dx
  
  # manage points outside the domain
  locations = abs(locations) # bottom and left
  locations[locations > 1] = 1-(locations[locations > 1]-1)
  
  b = rep(1, num_samples)
  
  set.seed(0)
  
  a1 = 1
  a2 = 1
  
  func_evaluation_nodes = numeric(nrow(nodes))
  func_evaluation_locations = numeric(nrow(locations))
  
  for (i in 1:nrow(nodes)){
    
    X_clean_nodes[i] = a1* cos(2*pi*nodes[i,1]) +
      a2* cos(2*pi*nodes[i,2]) + 1
    X_clean_locations[i] = a1* cos(2*pi*locations[i,1]) +
      a2* cos(2*pi*locations[i,2]) + 1
    
  }
  
  for (i in 1:nrow(nodes_sub)){
    
    X_clean_locations_sub[i] = a1* cos(2*pi*nodes_sub[i,1]) +
      a2* cos(2*pi*nodes_sub[i,2]) + 1
    
  }
  
  for(ii in 1:num_samples){
    
    noise = stats::rnorm(nrow(nodes), mean = 0, sd = 0.5)
    noise_sub = stats::rnorm(nrow(nodes_sub), mean = 0, sd = 0.5)
    X_nodes = rbind(X_nodes, X_clean_nodes + noise)
    X_locations = rbind(X_locations, X_clean_locations + noise)
    X_locations_sub = rbind(X_locations_sub, X_clean_locations_sub + noise_sub)
    
  }
  
  return(list(X_clean_nodes = X_clean_nodes,
              X_nodes = X_nodes,
              locations = locations,
              X_clean_locations = X_clean_locations,
              X_locations = X_locations,
              locations_sub = nodes_sub,
              X_clean_locations_sub = X_clean_locations_sub,
              X_locations_sub = X_locations_sub,
              b = b))
  
}


# ||||||||||||||||||
# Generate data ----
# ||||||||||||||||||

x <- seq(0, 1, length.out = 60)
y <- seq(0, 1, length.out = 60)

set.seed(0)
myL <- generate_2d_data(x, y, 50)

# Covariates
b <- myL[["b"]]

# At nodes
X_nodes <- myL[["X_nodes"]]
X_clean_nodes <- myL[["X_clean_nodes"]]

# At locations (#locations == #nodes)
locations <- myL[["locations"]]
X_locations <- myL[["X_locations"]]
X_clean_locations <- myL[["X_clean_locations"]]

# At locations (#locations < #nodes)
picked <- sample.int(dim(locations)[1], dim(locations)[1]*0.1)
locations_less <- locations[picked, ]
X_locations_less <- X_locations[, picked]
X_clean_locations_less <- X_clean_locations[picked]

# At locations (#locations < #nodes, equispaced subgrid)
locations_sub <- myL[["locations_sub"]]
X_locations_sub <- myL[["X_locations_sub"]]
X_clean_locations_sub <- myL[["X_clean_locations_sub"]]


# ||||||||||||||||
# Export data ----
# ||||||||||||||||

if (!file.exists(tests_dir)){
  dir.create(tests_dir)
}
  
## Test 1 ----
## |||||||||||

test_dir <- paste(tests_dir, "/2D_test1", sep = '')
if (!file.exists(test_dir)){
  dir.create(test_dir)
}

write.csv(X_nodes, paste(test_dir, "/X.csv", sep = ''))
write.csv(X_clean_nodes, paste(test_dir, "/X_clean.csv", sep = ''))

write.csv(b, paste(test_dir, "/b.csv", sep = ''))

## Test 2 ----
## |||||||||||

test_dir <- paste(tests_dir, "/2D_test2", sep = '')
if (!file.exists(test_dir)){
  dir.create(test_dir)
}

write.csv(X_clean_nodes, paste(test_dir, "/X_clean.csv", sep = ''))

write.csv(locations, paste(test_dir, "/locations.csv", sep = ''))
write.csv(X_locations, paste(test_dir, "/X_locations.csv", sep = ''))
write.csv(X_clean_locations, paste(test_dir, "/X_clean_locations.csv", sep = ''))

write.csv(b, paste(test_dir, "/b.csv", sep = ''))

## Test 3 ----
## |||||||||||

test_dir <- paste(tests_dir, "/2D_test3", sep = '')
if (!file.exists(test_dir)){
  dir.create(test_dir)
}

write.csv(X_clean_nodes, paste(test_dir, "/X_clean.csv", sep = ''))

write.csv(locations_less, paste(test_dir, "/locations.csv", sep = ''))
write.csv(X_locations_less, paste(test_dir, "/X_locations.csv", sep = ''))
write.csv(X_clean_locations_less, paste(test_dir, "/X_clean_locations.csv", sep = ''))

write.csv(b, paste(test_dir, "/b.csv", sep = ''))

## Test 4 ----
## |||||||||||

test_dir <- paste(tests_dir, "/2D_test4", sep = '')
if (!file.exists(test_dir)){
  dir.create(test_dir)
}

write.csv(X_clean_nodes, paste(test_dir, "/X_clean.csv", sep = ''))

write.csv(locations_sub, paste(test_dir, "/locations.csv", sep = ''))
write.csv(X_locations_sub, paste(test_dir, "/X_locations.csv", sep = ''))
write.csv(X_clean_locations_sub, paste(test_dir, "/X_clean_locations.csv", sep = ''))

write.csv(b, paste(test_dir, "/b.csv", sep = ''))

