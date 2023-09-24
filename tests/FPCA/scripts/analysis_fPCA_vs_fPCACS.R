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
library(RColorBrewer)


# ||||||||||||||
# Functions ----
# ||||||||||||||

load("../../utils/functions/tests_unit_square.RData")
load("../../utils/functions/plot_field.RData")
load("scripts/functions/generate_data.RData")
load("scripts/functions/generate_data_combinations.RData")
norm_vec <- function(x) sqrt(sum(x^2))


# ||||||||||||||||||||
# Test parameters ----
# ||||||||||||||||||||

# Data directory
data.directory <- "data/fPCA_vs_fPCACS/"
results.directory <- "results/fPCA_vs_fPCACS/"
images.directory <- "images/fPCA_vs_fPCACS/"


# Number of nodes (on the side)
n.nodes_vect <- c(10, 20, 30, 40, 50)

# Test options
N_vect <- c(50, 75, 100, 150, 200, 300, 400, 600, 800, 1200, 1600)

# Number of statistical units
N <- max(N_vect)
N.batches <- 10
N.generated <- N*N.batches

# To force the generation of the field
FORCE_GENERATION = TRUE


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
  n.nodes = n.nodes_vect
)

# List of test options names
test.options_names <- list(
  n.nodes = paste("nnodes", n.nodes_vect, sep = "")
)

test.options <- generate_data_combinations(N.generated, 
                                           test.options, test.options_names, 
                                           data.directory,
                                           FORCE_GENERATION)

cat("\n")


# |||||||||
# Test ----
# |||||||||

cat("\n# |||||||||||||")
cat("\n# Analysis ----")
cat("\n# |||||||||||||\n")


# Create results directory
if (!file.exists(results.directory)){
  dir.create(results.directory)
}

models <- c("fPCA", "fPCA_CS", "fPCA_CS_ML", "fPCA_CS_I", "fPCA_CS_ML_I")
n.nodes_vect_names <- paste("nv", n.nodes_vect, sep = "")

RMSE.data <- list()
RMSE.f1 <- list()
RMSE.f2 <- list()
RMSE.f3 <- list()
RMSE.s1 <- list()
RMSE.s2 <- list()
RMSE.s3 <- list()
PrincAngl <- list()
times <- list()
for(m in models){
  RMSE.data[[m]] <- list()
  RMSE.f1[[m]] <- list()
  RMSE.f2[[m]] <- list()
  RMSE.f3[[m]] <- list()
  RMSE.s1[[m]] <- list()
  RMSE.s2[[m]] <- list()
  RMSE.s3[[m]] <- list()
  PrincAngl[[m]] <- list()
  times[[m]] <- list()
  for(i in 1:length(n.nodes_vect)){
    RMSE.data[[m]][[n.nodes_vect_names[i]]] <- matrix(0, nrow = N.batches, ncol = length(N_vect))
    RMSE.f1[[m]][[n.nodes_vect_names[i]]] <- matrix(0, nrow = N.batches, ncol = length(N_vect))
    RMSE.f2[[m]][[n.nodes_vect_names[i]]] <- matrix(0, nrow = N.batches, ncol = length(N_vect))
    RMSE.f3[[m]][[n.nodes_vect_names[i]]] <- matrix(0, nrow = N.batches, ncol = length(N_vect))
    RMSE.s1[[m]][[n.nodes_vect_names[i]]] <- matrix(0, nrow = N.batches, ncol = length(N_vect))
    RMSE.s2[[m]][[n.nodes_vect_names[i]]] <- matrix(0, nrow = N.batches, ncol = length(N_vect))
    RMSE.s3[[m]][[n.nodes_vect_names[i]]] <- matrix(0, nrow = N.batches, ncol = length(N_vect))
    PrincAngl[[m]][[n.nodes_vect_names[i]]] <- matrix(0, nrow = N.batches, ncol = length(N_vect))
    times[[m]][[n.nodes_vect_names[i]]] <- matrix(0, nrow = N.batches, ncol = length(N_vect))
  }
}


for(i in 1:dim(test.options)[1]){
  
  cat(paste("\n\n\nNumber of nodes: ", test.options$test.name[i], "\n", sep= ""))
  
  # Load data
  load(test.options$data.path[i])
  
  # Set model
  pde <- new(Laplacian_2D_Order1, mesh_data)
  
  # Set zero forcing term
  quadrature_nodes <- pde$get_quadrature_nodes()
  f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
  pde$set_forcing_term(as.matrix(f))
  
  
  for(j in 1:length(N_vect)){
      
    N <- N_vect[j]
    
    cat(paste("\n- Number of statistical units ", N, ": ", sep= ""))
    
    model <- list()
    
    # Define and init model: fPCA
    model[["fPCA"]]<- new(FPCA_Laplacian_2D_GeoStatLocations, pde)
    model[["fPCA"]]$set_locations(locations)
    model[["fPCA"]]$set_lambdas(1e-4)
    
    # Define and init model: fPCA_CS
    model[["fPCA_CS"]] <- new(FPCA_CS_Laplacian_2D_GeoStatLocations, pde)
    model[["fPCA_CS"]]$set_locations(locations)
    model[["fPCA_CS"]]$set_lambda_s(1e-4)
    model[["fPCA_CS"]]$set_verbose(FALSE)
    model[["fPCA_CS"]]$set_mass_lumping(FALSE)
    model[["fPCA_CS"]]$set_iterative(FALSE)
    
    # Define and init model: fPCA_CS_ML
    model[["fPCA_CS_ML"]] <- new(FPCA_CS_Laplacian_2D_GeoStatLocations, pde)
    model[["fPCA_CS_ML"]]$set_locations(locations)
    model[["fPCA_CS_ML"]]$set_lambda_s(1e-4)
    model[["fPCA_CS_ML"]]$set_verbose(FALSE)
    model[["fPCA_CS_ML"]]$set_mass_lumping(TRUE)
    model[["fPCA_CS_ML"]]$set_iterative(FALSE)
    
    # Define and init model: fPCA_CS_I
    model[["fPCA_CS_I"]] <- new(FPCA_CS_Laplacian_2D_GeoStatLocations, pde)
    model[["fPCA_CS_I"]]$set_locations(locations)
    model[["fPCA_CS_I"]]$set_lambda_s(1e-4)
    model[["fPCA_CS_I"]]$set_verbose(FALSE)
    model[["fPCA_CS_I"]]$set_mass_lumping(FALSE)
    model[["fPCA_CS_I"]]$set_iterative(TRUE)
    
    # Define and init model: fPCA_CS_ML_I
    model[["fPCA_CS_ML_I"]] <- new(FPCA_CS_Laplacian_2D_GeoStatLocations, pde)
    model[["fPCA_CS_ML_I"]]$set_locations(locations)
    model[["fPCA_CS_ML_I"]]$set_lambda_s(1e-4)
    model[["fPCA_CS_ML_I"]]$set_verbose(FALSE)
    model[["fPCA_CS_ML_I"]]$set_mass_lumping(TRUE)
    model[["fPCA_CS_ML_I"]]$set_iterative(TRUE)
    
    for(m in models){
      
      cat(paste("\n  - Model ", m, ": \t\t", sep= ""))
    
      # Generate a grouping variable for split
      group_var <- rep(1:(N*N.batches), each = N, length.out = N*N.batches)
      
      # Split the vector into N subvectors without repetition
      set.seed(0)
      subvectors <- split(sample(1:(N*N.batches)), group_var)
      
      for(k in 1:N.batches){
        
        cat(paste("|", sep= ""))
        
        # Set observations
        selected <- subvectors[[k]]
        model[[m]]$set_observations(data_centered[selected,])
        
        # Solve
        tic(quiet = TRUE)
        model[[m]]$solve()
        elapsed <- toc(quiet = TRUE)
        
        # Results
        load <- list()
        score <- list()
        for(h in 1:3){
          load[[h]] <- model[[m]]$loadings()[, h]
          score[[h]] <- model[[m]]$scores()[, h]
          L2norm <- sqrt(as.numeric(t(load[[h]]) %*% R0 %*% load[[h]]))
          load[[h]] <- load[[h]] / L2norm
          score[[h]] <-  score[[h]] * L2norm
          if (norm_vec(load[[h]] + load_ex[[h]]) < norm_vec(load[[h]] - load_ex[[h]])) {
            load[[h]] = -load[[h]]
            score[[h]] = -score[[h]]
          }
        }
        
        # h <- 2
        # par(mfrow = c(1, 2))
        # plot_field(nodes, locations, load[[h]], range(load_ex[[h]]), "Estimated", FALSE)
        # plot_field(nodes, locations, load_ex[[h]], range(load_ex[[h]]), "Expected", FALSE)
        
        # Compute and store errors
        RMSE.f1[[m]][[n.nodes_vect_names[i]]][k,j] = norm_vec(load[[1]]-load_ex[[1]])/sqrt(4*n.nodes_vect[i]^2)
        RMSE.f2[[m]][[n.nodes_vect_names[i]]][k,j] = norm_vec(load[[2]]-load_ex[[2]])/sqrt(4*n.nodes_vect[i]^2)
        RMSE.f3[[m]][[n.nodes_vect_names[i]]][k,j] = norm_vec(load[[3]]-load_ex[[3]])/sqrt(4*n.nodes_vect[i]^2)
        RMSE.s1[[m]][[n.nodes_vect_names[i]]][k,j] = norm_vec(score[[1]]-score_ex[[1]][selected])/sqrt(N)
        RMSE.s2[[m]][[n.nodes_vect_names[i]]][k,j] = norm_vec(score[[2]]-score_ex[[2]][selected])/sqrt(N)
        RMSE.s3[[m]][[n.nodes_vect_names[i]]][k,j] = norm_vec(score[[3]]-score_ex[[3]][selected])/sqrt(N)
        
        # Reconstruction error
        data_projection <- (data - data_centered)[selected,] + # mean
                           (
                            matrix(score[[1]]) %*% t(load[[1]]) +
                            matrix(score[[2]]) %*% t(load[[2]]) +
                            matrix(score[[3]]) %*% t(load[[2]])
                            )
        RMSE.data[[m]][[n.nodes_vect_names[i]]][k,j] = sqrt(mean((data_projection - data_clean[selected,])^2))/(sqrt(4*n.nodes_vect[i]^2)*sqrt(N))
        PrincAngl[[m]][[n.nodes_vect_names[i]]][k,j] = subspace(cbind(load_ex[[1]], load_ex[[2]], load_ex[[3]]),
                                   cbind(load[[1]], load[[2]], load[[3]]))
        
        # Execution time
        times[[m]][[n.nodes_vect_names[i]]][k,j] <- elapsed$toc - elapsed$tic
        
      }
      
      cat(paste(" \t", mean(times[[m]][[n.nodes_vect_names[i]]][,j]), sep= ""))
      
    }
      
    rm(model)
    
  }
  
  rm(K, S, N,
     nodes, locations, mesh_data, R0,
     data, data_centered_clean, data_clean, data_centered,
     score_ex, load_ex)
  
}


save(data.directory, results.directory, images.directory,
     N_vect, n.nodes_vect, models,
     RMSE.data,
     RMSE.f1, RMSE.f2, RMSE.f3,
     RMSE.s1, RMSE.s2, RMSE.s3,
     PrincAngl, times,
     file = paste(results.directory, "results", format(Sys.time(), "_%Y%m%d_%H%M%S"), ".RData", sep = ""))

cat("\n\n\n")


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

list.files(results.directory)
load(paste(results.directory, tail(list.files(results.directory), n = 1), sep = ""))


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

plot_comparison <- function(comparison, title, mains, ylab, xticks, n.nodes_vect, models, colors, lines = FALSE) {
  
  par(mfrow = c(ceiling((length(n.nodes_vect)+1)/2),2), mar = c(4,4,5,1), oma = c(0,0,3,0))
  
  plot.new()
  legend("center", models, col = colors, pch = 16)
  
  for(n.nodes in n.nodes_vect){
    
    n_test <- dim(comparison[[models[1]]][[n.nodes]])[2]
    merged <- comparison[[models[1]]][[n.nodes]]
    n_test_options_1 <- 1
    for(m in models[-1]){
      merged <- merge_tests(merged, comparison[[m]][[n.nodes]], n_test, n_test_options_1, 1)
      n_test_options_1 <- n_test_options_1 + 1
    }
    fake = matrix(NaN, nrow = nrow(comparison[[m]][[n.nodes]]), ncol = ncol(comparison[[m]][[n.nodes]]))
    merged <- merge_tests(merged, fake, n_test, n_test_options_1, 1)
    
    ylim = c(0, max(merged, na.rm = TRUE))
    t_temp <- NULL
    x_temp <- NULL
    n_models <- length(models)
    iter <- 1
    for(m in models){
      t_temp <- rbind(t_temp, 0:(n_test-1)*(n_models+1)+iter)
      x_temp <- rbind(x_temp, colMeans(comparison[[m]][[n.nodes]]))
      iter <- iter + 1
    }
    if(lines){
      matplot(t(t_temp), 
              t(x_temp),
              main = mains[[n.nodes]],
              col = colors, type = "o", lwd = 2, pch = 19, lty = 1,
              xlab = "N", xaxt = "n", xlim = c(min(t_temp)-0.3, max(t_temp)+0.3),
              ylab = ylab, ylim = ylim)
    } else{
      matplot(NaN, NaN,
              main = mains[[n.nodes]],
              col = colors, type = "o", lwd = 2, pch = 19, lty = 1,
              xlab = "N", xaxt = "n", xlim = c(min(t_temp)-0.3, max(t_temp)+0.3),
              ylab = ylab, ylim = ylim)
    }
    
    boxplot(merged, col = c(colors, "white"), add = TRUE, xaxt = "n")
    
    axis(1, at=(1:(n_test))*(n_models+1)-(n_models)/2, labels=xticks)
    grid(nx = NA, ny = NULL)
    abline(v = 0:(n_test)*(n_models+1), col = "lightgrey", lty = 3)
    
  }
  
  title(title, outer = TRUE, line = 1, cex.main = 1.5)
  
}


## ||||||||||||
## Results ----
## ||||||||||||

models <- c("fPCA", "fPCA_CS", "fPCA_CS_ML", "fPCA_CS_I", "fPCA_CS_ML_I")
n.nodes_vect_names <- paste("nv", n.nodes_vect, sep = "")


# Choose a color palette from RColorBrewer
palette_name <- "Paired"

# Generate N different colors from the chosen palette
colors <- rev(brewer.pal(length(models), palette_name))

# Models names
models_names <- models
n.nodes_chosen <- (1:length(n.nodes_vect))[-2]

# Pictures settings
width = 2000
height_per_row = 650
pointsize = 40

# Execution time
jpeg(paste(images.directory, "time", ".jpg", sep = ""),
     width = width, height = height_per_row*length(n.nodes_vect)/2, units = "px", quality = 100, pointsize = pointsize,
     type = "windows")
plot_comparison(times,
                "Execution time",
                as.list(setNames(paste("K =", n.nodes_vect[n.nodes_chosen]^2), n.nodes_vect_names[n.nodes_chosen])),
                "Execution time [s]",
                N_vect, n.nodes_vect_names[n.nodes_chosen], models, colors, lines = TRUE)
dev.off()

# Pictures settings
width = 2000
height_per_row = 850
pointsize = 40

# Execution time (only new methods)
jpeg(paste(images.directory, "time_only_new", ".jpg", sep = ""),
     width = width, height = height_per_row*length(n.nodes_vect)/2, units = "px", quality = 100, pointsize = pointsize,
     type = "windows")
plot_comparison(times,
                "Execution time",
                as.list(setNames(paste("K =", n.nodes_vect[n.nodes_chosen]^2), n.nodes_vect_names[n.nodes_chosen])),
                "Execution time [s]",
                N_vect, n.nodes_vect_names[n.nodes_chosen], models[-1], colors[-1], lines = TRUE)
dev.off()


# Pictures settings
width = 2000
height_per_row = 800
pointsize = 40

# Time analysis in N
n_rows <- ceiling((length(n.nodes_vect[n.nodes_chosen])+1)/2)
jpeg(paste(images.directory, "time_analysis_N", ".jpg", sep = ""),
     width = width, height = height_per_row*n_rows, units = "px", quality = 100, pointsize = pointsize,
     type = "windows")
par(mfrow = c(n_rows,2), mar = c(4,4,5,1), oma = c(0,0,3,0))
plot.new()
legend("center", models, col = colors, pch = 16)
mains <- as.list(setNames(paste("K =", n.nodes_vect[n.nodes_chosen]^2), n.nodes_vect_names[n.nodes_chosen]))
for(n.nodes in n.nodes_vect_names[n.nodes_chosen]){
  mean_time <- list()
  for(m in models){
    mean_time[[m]] <- colMeans(times[[m]][[n.nodes]])
  }
  mean_time <- data.frame(mean_time)
  mean_time <- mean_time / mean_time[1,1]# as.data.frame(apply(mean_time, 1, function(x, factor) x / factor, as.numeric(mean_time[1,])))
  matplot(N_vect,
          mean_time, log="xy",
          type = "l", lty = 1, col = colors, lwd = 4,
          ylab = "log(time)", xlab = "log(N)",
          main = mains[[n.nodes]])
  matplot(N_vect,
          cbind((N_vect/N_vect[1])^(1/2),
                (N_vect/N_vect[1])^1,
                (N_vect/N_vect[1])^2), log="xy",
          type = "l", col = "grey", lty = c(2, 3, 4), lwd = 2,
          add = TRUE)
  grid()
  legend("topleft", legend = (c(bquote("N"^"1/2"), bquote("N"^"1"), bquote("N"^"2"))), 
         col = c(rep("grey",3)), lty = c(2,3,4), lwd = 3)
}
title("Time complexity analysis (in N)", outer = TRUE, line = 1, cex.main = 1.5)
dev.off()


# Time analysis in K
n_rows <- ceiling((length(models)+1)/2)
jpeg(paste(images.directory, "time_analysis_K", ".jpg", sep = ""),
     width = width, height = height_per_row*n_rows, units = "px", quality = 100, pointsize = pointsize,
     type = "windows")
par(mfrow = c(n_rows,2), mar = c(4,4,5,1), oma = c(0,0,3,0))
plot.new()
legend("center", models, col = colors, pch = 16)
mean_time <- list()
for(m in models){
  for(n.nodes in n.nodes_vect_names[n.nodes_chosen]){
    mean_time[[m]] <- rbind(mean_time[[m]], colMeans(times[[m]][[n.nodes]]))
  }
}
for(m in 1:length(models)){
  mean_time[[m]] <- mean_time[[m]]/mean_time[[m]][1,1]# as.data.frame(apply(mean_time[[m]], 1, function(x, factor) x / factor, as.numeric(mean_time[[m]][1,])))
  n_total_nodes <- n.nodes_vect[n.nodes_chosen]^2
  matplot(n_total_nodes,
          mean_time[[m]], 
          type = "l", lty = 1, col = colors[m], lwd = 4, log="xy",
          xlab = "log(K)", ylab = "log(time)", main = models_names[m])
  matplot(n_total_nodes,
          cbind((n_total_nodes/n_total_nodes[1])^1,
                (n_total_nodes/n_total_nodes[1])^2,
                (n_total_nodes/n_total_nodes[1])^3),
          type = "l", col = "grey", lty = c(2,3,4), lwd = 2,
          add = TRUE)
  grid()
  legend("topleft", legend = (c(bquote("K"^"1"), bquote("K"^"2"), bquote("K"^"3"))), 
         col = c(rep("grey",3)), lty = c(2,3,4), lwd = 3)
}
title("Time complexity analysis (in K)", outer = TRUE, line = 1, cex.main = 1.5)
dev.off()

# Time analysis in K (normalized)
n_rows <- ceiling((length(models)+1)/2)
jpeg(paste(images.directory, "time_analysis_K_norm", ".jpg", sep = ""),
     width = width, height = height_per_row*n_rows, units = "px", quality = 100, pointsize = pointsize,
     type = "windows")
par(mfrow = c(n_rows,2), mar = c(4,4,5,1), oma = c(0,0,3,0))
plot.new()
legend("center", models, col = colors, pch = 16)
mean_time <- list()
for(m in models){
  for(n.nodes in n.nodes_vect_names[n.nodes_chosen]){
    mean_time[[m]] <- rbind(mean_time[[m]], colMeans(times[[m]][[n.nodes]]))
  }
}
for(m in 1:length(models)){
  mean_time[[m]] <- t(as.data.frame(apply(mean_time[[m]], 1, function(x, factor) x / factor, as.numeric(mean_time[[m]][1,])))) # mean_time[[m]]/mean_time[[m]][1,1]#
  n_total_nodes <- n.nodes_vect[n.nodes_chosen]^2
  matplot(n_total_nodes,
          mean_time[[m]], 
          type = "l", lty = 1, col = colors[m], lwd = 4, log="xy",
          xlab = "log(K)", ylab = "log(time)", main = models_names[m])
  matplot(n_total_nodes,
          cbind((n_total_nodes/n_total_nodes[1])^1,
                (n_total_nodes/n_total_nodes[1])^2,
                (n_total_nodes/n_total_nodes[1])^3),
          type = "l", col = "grey", lty = c(2,3,4), lwd = 2,
          add = TRUE)
  grid()
  legend("topleft", legend = (c(bquote("K"^"1"), bquote("K"^"2"), bquote("K"^"3"))), 
         col = c(rep("grey",3)), lty = c(2,3,4), lwd = 3)
}
title("Time complexity analysis (in K)", outer = TRUE, line = 1, cex.main = 1.5)
dev.off()


# Pictures settings
width = 2000
height_per_row = 650
pointsize = 40

# Reconstruction error
jpeg(paste(images.directory, "reconstruction_error", ".jpg", sep = ""),
     width = width, height = height_per_row*length(n_total_nodes)/2, units = "px", quality = 100, pointsize = pointsize,
     type = "windows")
plot_comparison(RMSE.data,
                "Reconstruction error",
                as.list(setNames(paste("K =", n.nodes_vect[n.nodes_chosen]^2), n.nodes_vect_names[n.nodes_chosen])),
                "RMSE",
                N_vect, n.nodes_vect_names[n.nodes_chosen], models_names, colors)
dev.off()


# Pictures settings
width = 2000
height_per_row = 900
pointsize = 40

# Loading functions
jpeg(paste(images.directory, "load1", ".jpg", sep = ""),
     width = width, height = height_per_row*length(n_total_nodes)/2, units = "px", quality = 100, pointsize = pointsize,
     type = "windows")
plot_comparison(RMSE.f1,
                paste("First functional principal component error"),
                as.list(setNames(paste("K =", n.nodes_vect[n.nodes_chosen]^2), n.nodes_vect_names[n.nodes_chosen])),
                "RMSE",
                N_vect, n.nodes_vect_names[n.nodes_chosen], models_names, colors)
dev.off()
jpeg(paste(images.directory, "load2", ".jpg", sep = ""),
     width = width, height = height_per_row*length(n_total_nodes)/2, units = "px", quality = 100, pointsize = pointsize,
     type = "windows")
plot_comparison(RMSE.f2,
                paste("Second functional principal component error"),
                as.list(setNames(paste("K =", n.nodes_vect[n.nodes_chosen]^2), n.nodes_vect_names[n.nodes_chosen])),
                "RMSE",
                N_vect, n.nodes_vect_names[n.nodes_chosen], models_names, colors)
dev.off()
jpeg(paste(images.directory, "load3", ".jpg", sep = ""),
     width = width, height = height_per_row*length(n_total_nodes)/2, units = "px", quality = 100, pointsize = pointsize,
     type = "windows")
plot_comparison(RMSE.f3,
                paste("Third functional principal component error"),
                as.list(setNames(paste("K =", n.nodes_vect[n.nodes_chosen]^2), n.nodes_vect_names[n.nodes_chosen])),
                "RMSE",
                N_vect, n.nodes_vect_names[n.nodes_chosen], models_names, colors)
dev.off()


# Scores
jpeg(paste(images.directory, "score1", ".jpg", sep = ""),
     width = width, height = height_per_row*length(n_total_nodes)/2, units = "px", quality = 100, pointsize = pointsize,
     type = "windows")
plot_comparison(RMSE.s1,
                paste("First componet scores error"),
                as.list(setNames(paste("K =", n.nodes_vect[n.nodes_chosen]^2), n.nodes_vect_names[n.nodes_chosen])),
                "RMSE",
                N_vect, n.nodes_vect_names[n.nodes_chosen], models_names, colors)
dev.off()
jpeg(paste(images.directory, "score2", ".jpg", sep = ""),
     width = width, height = height_per_row*length(n_total_nodes)/2, units = "px", quality = 100, pointsize = pointsize,
     type = "windows")
plot_comparison(RMSE.s2,
                paste("Second componet scores error"),
                as.list(setNames(paste("K =", n.nodes_vect[n.nodes_chosen]^2), n.nodes_vect_names[n.nodes_chosen])),
                "RMSE",
                N_vect, n.nodes_vect_names[n.nodes_chosen], models_names, colors)
dev.off()
jpeg(paste(images.directory, "score3", ".jpg", sep = ""),
     width = width, height = height_per_row*length(n_total_nodes)/2, units = "px", quality = 100, pointsize = pointsize,
     type = "windows")
plot_comparison(RMSE.s3,
                paste("Third componet scores error"),
                as.list(setNames(paste("K =", n.nodes_vect[n.nodes_chosen]^2), n.nodes_vect_names[n.nodes_chosen])),
                "RMSE",
                N_vect, n.nodes_vect_names[n.nodes_chosen], models_names, colors)
dev.off()


# PrincAngl
# jpeg(paste(images.directory, "princangl", ".jpg", sep = ""),
#      width = width, height = height_per_row*length(n_total_nodes)/2, units = "px", quality = 100, pointsize = pointsize,
#      type = "windows")
# plot_comparison(PrincAngl,
#                 bquote(paste("PrincAngl")),
#                 as.list(setNames(paste("K =", n.nodes_vect[n.nodes_chosen]^2), n.nodes_vect_names[n.nodes_chosen])),
#                 "",
#                 N_vect, n.nodes_vect_names[n.nodes_chosen], models_names, colors)
# dev.off()


## |||||||||
## Data ----
## |||||||||


images.directory <- paste(images.directory, "data/", sep = "")
if (!file.exists(images.directory)){
  dir.create(images.directory)
}

pointsize <- 25

for(i in c(5)){ # dim(test.options)[1]
  
  cat(paste("\n\n\nNumber of nodes: ", test.options$test.name[i], "\n", sep= ""))
  
  # Load data
  load(test.options$data.path[i])
  
  jpeg(paste(images.directory, "loading_ex.jpg", sep = ""),
       width = 2000, height = 600, units = "px", quality = 100, pointsize = pointsize,
       type = "windows")
  par(mfrow = c(1,3))
  plot_field(nodes, locations, load_ex[[1]], range(load_ex[[1]]), title = paste("First exact functional principal component"), UNIQUE = FALSE)
  plot_field(nodes, locations, load_ex[[2]], range(load_ex[[2]]), title = paste("Second exact functional principal component"), UNIQUE = FALSE)
  plot_field(nodes, locations, load_ex[[3]], range(load_ex[[3]]), title = paste("Third exact functional principal component"), UNIQUE = FALSE)
  dev.off()
  
  # jpeg(paste(images.directory, "scores_ex.jpg", sep = ""),
  #      width = 2000, height = 700, units = "px", quality = 100, pointsize = pointsize,
  #      type = "windows")
  # boxplot(score_ex, main = "Exact scores distributions", xlab = "Components")
  # grid()
  # dev.off()
  
  jpeg(paste(images.directory, "data_examples.jpg", sep = ""),
       width = 2000, height = 1800, units = "px", quality = 100, pointsize = pointsize,
       type = "windows")
  par(mfrow = c(3,3))
  for(j in 1:9)
    plot_field(nodes, locations, data[j,], range(data[j,]), title = paste("Sample", j), UNIQUE = FALSE)
  dev.off()
}




