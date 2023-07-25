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
n.nodes_vect <- c(20, 30, 40, 50, 60)

# Test options
N_vect <- c(50, 75, 100, 150, 200, 300, 400, 600, 800, 1200)

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

models <- c("fPCA", "fPCA_CS", "fPCA_CS_ML")

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
    RMSE.data[[m]][[i]] <- matrix(0, nrow = N.batches, ncol = length(N_vect))
    RMSE.f1[[m]][[i]] <- matrix(0, nrow = N.batches, ncol = length(N_vect))
    RMSE.f2[[m]][[i]] <- matrix(0, nrow = N.batches, ncol = length(N_vect))
    RMSE.f3[[m]][[i]] <- matrix(0, nrow = N.batches, ncol = length(N_vect))
    RMSE.s1[[m]][[i]] <- matrix(0, nrow = N.batches, ncol = length(N_vect))
    RMSE.s2[[m]][[i]] <- matrix(0, nrow = N.batches, ncol = length(N_vect))
    RMSE.s3[[m]][[i]] <- matrix(0, nrow = N.batches, ncol = length(N_vect))
    PrincAngl[[m]][[i]] <- matrix(0, nrow = N.batches, ncol = length(N_vect))
    times[[m]][[i]] <- matrix(0, nrow = N.batches, ncol = length(N_vect))
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
    model[["fPCA_CS"]]$set_mass_lumping(FALSE)
    
    # Define and init model: fPCA_CS
    model[["fPCA_CS_ML"]] <- new(FPCA_CS_Laplacian_2D_GeoStatLocations, pde)
    model[["fPCA_CS_ML"]]$set_locations(locations)
    model[["fPCA_CS_ML"]]$set_lambda_s(1e-4)
    model[["fPCA_CS_ML"]]$set_mass_lumping(TRUE)
    
    for(m in models){
      
      cat(paste("\n- - Model ", m, ": \t\t", sep= ""))
    
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
        RMSE.f1[[m]][[i]][k,j] = norm_vec(load[[1]]-load_ex[[1]])/sqrt(4*n.nodes_vect[i]^2)
        RMSE.f2[[m]][[i]][k,j] = norm_vec(load[[2]]-load_ex[[2]])/sqrt(4*n.nodes_vect[i]^2)
        RMSE.f3[[m]][[i]][k,j] = norm_vec(load[[3]]-load_ex[[3]])/sqrt(4*n.nodes_vect[i]^2)
        RMSE.s1[[m]][[i]][k,j] = norm_vec(score[[1]]-score_ex[[1]][selected])/sqrt(N)
        RMSE.s2[[m]][[i]][k,j] = norm_vec(score[[2]]-score_ex[[2]][selected])/sqrt(N)
        RMSE.s3[[m]][[i]][k,j] = norm_vec(score[[3]]-score_ex[[3]][selected])/sqrt(N)
        
        # Reconstruction error
        data_projection <- (data - data_centered)[selected,] + # mean
                           (
                            matrix(score[[1]]) %*% t(load[[1]]) +
                            matrix(score[[2]]) %*% t(load[[2]]) +
                            matrix(score[[3]]) %*% t(load[[2]])
                            )
        RMSE.data[[m]][[i]][k,j] = sqrt(mean((data_projection - data_clean[selected,])^2))
        PrincAngl[[m]][[i]][k,j] = subspace(cbind(load_ex[[1]], load_ex[[2]], load_ex[[3]]),
                                   cbind(load[[1]], load[[2]], load[[3]]))
        
        # Execution time
        times[[m]][[i]][k,j] <- elapsed$toc - elapsed$tic
        
      }
      
      cat(paste(" \t", mean(times[[m]][[i]][,j]), sep= ""))
      
    }
      
    rm(model)
    
  }
  
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

plot_comparison <- function(comparison, title, mains, ylab, xticks, models, colors) {
  
  par(mfrow = c(length(n.nodes_vect),1), mar = c(4,4,5,1), oma = c(0,0,3,0))
  
  for(i in 1:length(n.nodes_vect)){
    
    n_test <- dim(comparison[["fPCA_CS"]][[i]])[2]
    merged <- merge_tests(comparison[["fPCA"]][[i]], comparison[["fPCA_CS"]][[i]], n_test, 1, 1)
    merged <- merge_tests(merged, comparison[["fPCA_CS_ML"]][[i]], n_test, 2, 1)
    
    ylim = c(0, max(merged))
    matplot(t(rbind(0:9*3+1, 0:9*3+2, 0:9*3+3)), 
            t(rbind(colMeans(comparison[["fPCA"]][[i]]),
                    colMeans(comparison[["fPCA_CS"]][[i]]),
                    colMeans(comparison[["fPCA_CS_ML"]][[i]]))),
            main = mains[i],
            col = colors, type = "o", lwd = 2, pch = 19, lty = 1,
            xlab = "# of statistical units", xaxt = "n",
            ylab = ylab, ylim = ylim)
    boxplot(merged, col = colors, add = TRUE, xaxt = "n")
    n_tests <- dim(merged)[2]/length(models)
    axis(1, at=(1:n_tests)*(length(models))-(length(models)-1)/2, labels=xticks)
    grid(nx = NA, ny = NULL)
    abline(v = 0:9*3+0.5, col = "lightgrey", lty = 3)
    # legend("topleft", models, col = colors, pch = 16)
    
  }
  
  title(title, outer = TRUE, line = 1, cex.main = 1.5)
  
}



# Choose a color palette from RColorBrewer
palette_name <- "Paired"

# Generate N different colors from the chosen palette
colors <- rev(brewer.pal(length(models), palette_name))

# Models names
models_names <- c("fPCA", "fPCA + CS", "fPCA + CS + ML")

# Pictures settings
width = 3000
height = 6000
pointsize = 80

# Execution time
jpeg(paste(images.directory, "time", ".jpg", sep = ""),
     width = width, height = height, units = "px", quality = 100, pointsize = pointsize,
     type = "windows")
plot_comparison(times,
                "Execution time",
                paste("# of nodes:", n.nodes_vect^2),
                "Execution time [s]",
                N_vect, models_names, colors)
dev.off()

# Time analysis in N
jpeg(paste(images.directory, "time_analysis_N", ".jpg", sep = ""),
     width = width, height = height*2/3, units = "px", quality = 100, pointsize = pointsize,
     type = "windows")
par(mfrow = c(ceiling(length(n.nodes_vect)/2),2), mar = c(4,4,5,1), oma = c(0,0,3,0))
mains <- paste("# of nodes:", n.nodes_vect^2)
for(i in 1:length(n.nodes_vect)){
  mean_time <- list()
  mean_time[["fPCA"]] <- colMeans(times[["fPCA"]][[i]])
  mean_time[["fPCA_CS"]] <- colMeans(times[["fPCA_CS"]][[i]])
  mean_time[["fPCA_CS_ML"]] <- colMeans(times[["fPCA_CS_ML"]][[i]])
  matplot(N_vect,
          cbind((N_vect/N_vect[1])^(1/2),
                (N_vect/N_vect[1])^1,
                (N_vect/N_vect[1])^2), log="xy",
          type = "l", col = "grey", lty = c(2,3,4), lwd = 2,
          ylab = "log(time)", xlab = "log(# of statistical units)",
          main = mains[i])
  mean_time <- data.frame(mean_time)
  mean_time <- as.data.frame(apply(mean_time, 1, function(x, factor) x / factor, as.numeric(mean_time[1,])))
  matplot(N_vect,
          t(mean_time), log="xy",
          type = "l", lty = 1, col = colors, lwd = 3,
          add = TRUE)
  grid()
}
plot.new()
legend("center", legend = (c(models, bquote("N"^"1/2"), bquote("N"^"1"), bquote("N"^"2"))), 
       col = c(colors, rep("grey",3)), lty = c(1,1,1,2,3,4), lwd = 3)
title("Computational cost analysis (in N)", outer = TRUE, line = 1, cex.main = 1.5)
dev.off()


# Time analysis in n.nodes
jpeg(paste(images.directory, "time_analysis_nnodes", ".jpg", sep = ""),
     width = width, height = height*3/5, units = "px", quality = 100, pointsize = pointsize,
     type = "windows")
par(mfrow = c(2,2), mar = c(4,4,5,1), oma = c(0,0,3,0))
mean_time <- list()
for(i in 1:length(n.nodes_vect)){
  mean_time[["fPCA"]] <- rbind(mean_time[["fPCA"]], colMeans(times[["fPCA"]][[i]]))
  mean_time[["fPCA_CS"]] <- rbind(mean_time[["fPCA_CS"]], colMeans(times[["fPCA_CS"]][[i]]))
  mean_time[["fPCA_CS_ML"]] <- rbind(mean_time[["fPCA_CS_ML"]], colMeans(times[["fPCA_CS_ML"]][[i]]))
}
for(m in 1:3){
  mean_time[[m]] <- as.data.frame(apply(mean_time[[m]], 1, function(x, factor) x / factor, as.numeric(mean_time[[m]][1,])))
  n_total_nodes <- n.nodes_vect^2
  matplot(n_total_nodes,
          cbind((n_total_nodes/n_total_nodes[1])^1,
                (n_total_nodes/n_total_nodes[1])^2,
                (n_total_nodes/n_total_nodes[1])^3), log="xy",
          type = "l", col = "grey", lty = c(2,3,4), lwd = 2,
          xlab = "log(# of nodes)", ylab = "log(time)", main = models_names[m])
  matplot(n_total_nodes,
          t(mean_time[[m]]), 
          type = "l", lty = 1, col = colors[m], lwd = 3,
          add = TRUE)
  grid()
}
plot.new()
legend("center", legend = (c(models, bquote("n.nodes"^"1"), bquote("n.nodes"^"2"), bquote("n.nodes"^"3"))), 
       col = c(colors, rep("grey",3)), lty = c(1,1,1,2,3,4), lwd = 3)
title("Computational cost analysis (in # of nodes)", outer = TRUE, line = 1, cex.main = 1.5)
dev.off()


# Reconstruction error
jpeg(paste(images.directory, "reconstruction_error", ".jpg", sep = ""),
     width = width, height = height, units = "px", quality = 100, pointsize = pointsize,
     type = "windows")
plot_comparison(RMSE.data,
                "Reconstruction error",
                paste("# of nodes:", n.nodes_vect^2),
                "MSE",
                N_vect, models_names, colors)
dev.off()


# Loading functions
jpeg(paste(images.directory, "load1", ".jpg", sep = ""),
     width = width, height = height, units = "px", quality = 100, pointsize = pointsize,
     type = "windows")
plot_comparison(RMSE.f1,
                bquote(paste("1"^"st"," loading function error")),
                paste("# of nodes:", n.nodes_vect^2),
                "MSE",
                N_vect, models_names, colors)
dev.off()
jpeg(paste(images.directory, "load2", ".jpg", sep = ""),
     width = width, height = height, units = "px", quality = 100, pointsize = pointsize,
     type = "windows")
plot_comparison(RMSE.f2,
                bquote(paste("2"^"nd"," loading function error")),
                paste("# of nodes:", n.nodes_vect^2),
                "MSE",
                N_vect, models_names, colors)
dev.off()
jpeg(paste(images.directory, "load3", ".jpg", sep = ""),
     width = width, height = height, units = "px", quality = 100, pointsize = pointsize,
     type = "windows")
plot_comparison(RMSE.f3,
                bquote(paste("3"^"rd"," loading function error")),
                paste("# of nodes:", n.nodes_vect^2),
                "MSE",
                N_vect, models_names, colors)
dev.off()


# Scores
jpeg(paste(images.directory, "score1", ".jpg", sep = ""),
     width = width, height = height, units = "px", quality = 100, pointsize = pointsize,
     type = "windows")
plot_comparison(RMSE.s1,
                bquote(paste("1"^"st"," score error")),
                paste("# of nodes:", n.nodes_vect^2),
                "MSE",
                N_vect, models_names, colors)
dev.off()
jpeg(paste(images.directory, "score2", ".jpg", sep = ""),
     width = width, height = height, units = "px", quality = 100, pointsize = pointsize,
     type = "windows")
plot_comparison(RMSE.s2,
                bquote(paste("2"^"nd"," score error")),
                paste("# of nodes:", n.nodes_vect^2),
                "MSE",
                N_vect, models_names, colors)
dev.off()
jpeg(paste(images.directory, "score3", ".jpg", sep = ""),
     width = width, height = height, units = "px", quality = 100, pointsize = pointsize,
     type = "windows")
plot_comparison(RMSE.s3,
                bquote(paste("3"^"rd"," score error")),
                paste("# of nodes:", n.nodes_vect^2),
                "MSE",
                N_vect, models_names, colors)
dev.off()


# PrincAngl
jpeg(paste(images.directory, "princangl", ".jpg", sep = ""),
     width = width, height = height, units = "px", quality = 100, pointsize = pointsize,
     type = "windows")
plot_comparison(PrincAngl,
                bquote(paste("PrincAngl")),
                paste("# of nodes:", n.nodes_vect^2),
                "",
                N_vect, models_names, colors)
dev.off()



  




