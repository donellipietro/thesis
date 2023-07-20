# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Test performance comparison between fPCA and fPCA_CS %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()
def.par = par()

cat("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
cat("\n%% Test performance comparison between fPCA and fPCA_CS %%")
cat("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")


# ||||||||||||||
# Libraries ----
# ||||||||||||||

library(fdaPDE2)
library(pracma)
library(plot3D)
library(tictoc)


# ||||||||||||||
# Functions ----
# ||||||||||||||

load("../../utils/tests_unit_square.RData")
load("../../utils/plot_field.RData")
load("scripts/functions/generate_data.RData")


# ||||||||||||||||||||
# Test parameters ----
# ||||||||||||||||||||

# Number of nodes (on the side)
n.nodes <- 31

# Number of locations (on the side)
n.locations <- 60

# Test options
N_vect <- c(50, 75, 100, 150, 200, 300, 400, 600, 800, 1200)

# Number of statistical units
N <- max(N_vect)


# |||||||||
# Data ----
# |||||||||

data.name <- paste("data_", N, ".RData", sep = "")
data.path <- paste("data/", data.name, sep = "")

if(!file.exists(data.path)){
  cat("\n# ||||||||||||||||||||")
  cat("\n# Generating data ----")
  cat("\n# ||||||||||||||||||||\n")
  generate_data(n.nodes, n.locations, N, data.path)
}

cat("\n# |||||||||||||||||")
cat("\n# Loading data ----")
cat("\n# |||||||||||||||||\n")
load(data.path)



# |||||||||
# Test ----
# |||||||||

cat("\n# |||||||||||||")
cat("\n# Analysis ----")
cat("\n# |||||||||||||\n")

# Set model
pde <- new(Laplacian_2D_Order1, mesh_data)

# set zero forcing term
quadrature_nodes <- pde$get_quadrature_nodes()
f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
pde$set_forcing_term(as.matrix(f))

# Model options
lambda_s <- 10^-3

# Room for results
times <- matrix(0, nrow = 3, ncol = length(N_vect))
errors <- matrix(0, nrow = 3, ncol = length(N_vect))

for(i in 1:length(N_vect)){
  
  cat(paste("\n\nExecution: N = ", N_vect[i], "\n\n", sep = ""))
  
  ## define and init model_FPCA
  model_FPCA <- new(FPCA_Laplacian_2D_GeoStatLocations, pde)
  model_FPCA$set_locations(locations)
  model_FPCA$set_lambdas(lambda_s)
  model_FPCA$init_regularization()
  
  ## define and init model_FPCA_CS
  model_FPCA_CS <- new(FPCA_CS_Laplacian_2D_GeoStatLocations, pde)
  model_FPCA_CS$set_locations(locations)
  model_FPCA_CS$set_lambda_s(lambda_s)
  model_FPCA_CS$init_regularization()
  
  ## define and init model_FPCA_CS + Mass lumping
  model_FPCA_CS_ML <- new(FPCA_CS_Laplacian_2D_GeoStatLocations, pde)
  model_FPCA_CS_ML$set_locations(locations)
  model_FPCA_CS_ML$set_lambda_s(lambda_s)
  model_FPCA_CS_ML$set_mass_lumping(1)
  model_FPCA_CS_ML$init_regularization()
  
  # Sample data
  set.seed(i)
  sample <- sample.int(N, N_vect[i], replace =  TRUE)
  
  # Set data
  model_FPCA$set_observations(data[sample,])
  model_FPCA_CS$set_observations(data[sample,])
  model_FPCA_CS_ML$set_observations(data[sample,])
  
  # output_file <- paste("execution", exec, ".txt")
  # sink(output_file)
  
  # fPCA
  cat("\nfPCA\n")
  tic()
  model_FPCA$solve()
  elapsed <- toc()
  times[1, i] <- as.numeric(elapsed$toc - elapsed$tic)
  errors[1, i] <- norm(data_clean[sample,] - model_FPCA$scores() %*% t(model_FPCA$loadings()), "2")/(N_vect[i]*S)
  
  # fPCA_CS
  cat("\nfPCA - Closed solution\n")
  tic()
  model_FPCA_CS$solve()
  elapsed <- toc()
  times[2, i] <- as.numeric(elapsed$toc - elapsed$tic)
  errors[2, i] <- norm(data_clean[sample,] - model_FPCA_CS$scores() %*% t(model_FPCA_CS$loadings()), "2")/(N_vect[i]*S)
  
  # fPCA_CS + Mass Lumping
  cat("\nfPCA - Closed solution - Mass Lumping\n")
  tic()
  model_FPCA_CS_ML$solve()
  elapsed <- toc()
  times[3, i] <- as.numeric(elapsed$toc - elapsed$tic)
  errors[3, i] <- norm(data_clean[sample,] - model_FPCA_CS_ML$scores() %*% t(model_FPCA_CS_ML$loadings()), "2")/(N_vect[i]*S)
 
  # Clean
  rm(model_FPCA)
  rm(model_FPCA_CS)
  rm(model_FPCA_CS_ML)
}


save(data.path,
     N, K, S,
     N_vect, times, errors,
     file = paste("results/results", format(Sys.time(), "_%Y%m%d_%H%M%S"), ".RData", sep = ""))


cat("\n\n\n")

# ||||||||||||
# Results ----
# ||||||||||||

# list.files("../results/")
load(paste("results/", tail(list.files("results/"), n = 1), sep = ""))

jpeg(paste("images/fPCA_vsfPCAcs.jpg", sep = ""), width = 2000, height = 1500, units = "px", quality = 100, pointsize = 37)
par(mfrow = c(1, 2), mar = def.par$mar, omi = def.par$omi)
matplot(log(N_vect), t(times),
        type = "o", col = c("red", "blue", "green"),
        lty = c(1,1,1), pch = c(1,2,3),
        lwd = c(2,2,2),
        main = "Execution time",
        ylab = "Execution time [Seconds]",
        xlab = "log(# of statistical units)")
grid()
legend("topleft", c("fPCA", "fPCA CS", "fPCA CS + ML"),
       col = c("red", "blue", "green"), lty = c(1, 1, 1), pch = c(1,2,3))
matplot(log(N_vect), t(errors),
        type = "o", col = c("red", "blue", "green"),
        lty = c(1,1,1), pch = c(1,2,3),
        lwd = c(2,2,2),
        main = "MSE", 
        ylab = "MSE",
        xlab = "log(# of statistical units)")
grid()
dev.off()

