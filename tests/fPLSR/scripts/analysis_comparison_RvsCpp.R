# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Test fPLSR: noise R vs C++ %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()
def.par = par()

cat("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
cat("\n%% Test fPLSR: noise R vs C++ %%")
cat("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n")


# ||||||||||||||
# Libraries ----
# ||||||||||||||

library(fdaPDE)
library(fdaPDE2)
library(pracma)
library(plot3D)
library(tictoc)
library(RColorBrewer)
library(pR1FPLS)


# ||||||||||||||
# Functions ----
# ||||||||||||||

load("scripts/functions/generate_data_RvsCpp.RData")
load("scripts/functions/generate_data_combinations.RData")
load("scripts/functions/PLSR.RData")
load("../../utils/functions/import_fdaPDE_mesh.RData")
load("../../utils/functions/tests_unit_square.RData")
norm_vec <- function(x) sqrt(sum(x^2))


# ||||||||||||||||||||
# Test parameters ----
# ||||||||||||||||||||

# Data directory
mesh.directory <- "../../fdaPDE/test/data/mesh/"
data.directory <- "data/RvsCpp/"
results.directory <- "results/RvsCpp/"
images.directory <- "images/RvsCpp/"

# Mesh
n.nodes_vect <-  c(30, 40, 60)
n.nodes_vect_names <- paste("nnodes", n.nodes_vect, sep = "")

# Noise
L = 1
options.global = list(
  number.locations_per_element = 1,
  L = L,
  STRATEGY = 1,
  NSR.Y = (1/3)^2
)

# Number of statistical units
N_vect <- c(50, 75, 100, 150, 200, 300, 400)
N.batches <- 10
N.generated <- max(N_vect*N.batches)

# Data options
X.index_vect <- c(1) # 2
X.index_vect_names = paste("xi", X.index_vect, sep = "")

# Responses
B.index_vect <- c(1) #, 2, 3, 4, 5) 
B.index_vect_names = paste("bi", B.index_vect, sep = "")

# Options
FORCE_GENERATION = FALSE


# |||||||||
# Data ----
# |||||||||

cat("\n")
cat("\n# ||||||||||||||||||||")
cat("\n# Generating data ----")
cat("\n# ||||||||||||||||||||\n\n")

# Create test directory
if (!file.exists(data.directory)){
  dir.create(data.directory)
}

# List of test options
test.options <- list(
  n.nodes = n.nodes_vect,
  X.index = X.index_vect,
  B.index = B.index_vect
)

# List of test options names
test.options_names <- list(
  n.nodes = n.nodes_vect_names,
  X.index = X.index_vect_names,
  B.index = B.index_vect_names
)

test.options <- generate_data_combinations(N.generated, 
                                           test.options, test.options_names,
                                           options.global,
                                           data.directory,
                                           FORCE_GENERATION)

cat("\n")


# |||||||||
# Test ----
# |||||||||

# Penalty
lambda_s <- 10

cat("\n")
cat("\n# |||||||||")
cat("\n# Test ----")
cat("\n# |||||||||\n")

# Create results directory
if (!file.exists(results.directory)){
  dir.create(results.directory)
}

models_R <- c("R", "R m")
models_cpp <- c("C++ ff", "C++ m", "C++ s")
models <- c("Multivariate", models_R, models_cpp)
results_names <- c("Y_hat", "X_hat", "B_hat")

options.models_R <- list()
options.models_R[["R"]] <- list(lambda_s = 1e-5)
options.models_R[["R m"]] <- list(lambda_s = 1e-10)

options.models_cpp <- list()
options.models_cpp[["C++ ff"]] <- list(full_functional = TRUE, smoothing = FALSE, lambda_s = 1e-5)
options.models_cpp[["C++ m"]] <- list(full_functional = FALSE, smoothing = FALSE, lambda_s = 1e-16)
options.models_cpp[["C++ s"]] <- list(full_functional = FALSE, smoothing = TRUE, lambda_s = 1e-7)


times <- list()
for(m in models){
  times[[m]] <- list()
  for(i in 1:length(X.index_vect)){
    times[[m]][[X.index_vect_names[i]]] <- list()
    for(j in 1:length(B.index_vect)){
      times[[m]][[X.index_vect_names[i]]][[B.index_vect_names[j]]] <- list()
      for(k in 1:length(n.nodes_vect)){
        times[[m]][[X.index_vect_names[i]]][[B.index_vect_names[j]]][[n.nodes_vect_names[k]]] <- matrix(0, nrow = N.batches, ncol = length(N_vect))
      }
    }
  }
}

error <- list()
for(r in results_names){
  error[[r]] <- list()
  for(m in models){
    error[[r]][[m]] <- list()
    for(i in 1:length(X.index_vect)){
      error[[r]][[m]][[X.index_vect_names[i]]] <- list()
      for(j in 1:length(B.index_vect)){
        error[[r]][[m]][[X.index_vect_names[i]]][[B.index_vect_names[j]]] <- list()
        for(k in 1:length(n.nodes_vect)){
          error[[r]][[m]][[X.index_vect_names[i]]][[B.index_vect_names[j]]][[n.nodes_vect_names[k]]] <- matrix(0, nrow = N.batches, ncol = length(N_vect))
        }
      }
    }
  }
}

diff <- list()
for(r in results_names){
  diff[[r]] <- list()
  for(m in models[-c(1,2,3)]){
    diff[[r]][[m]] <- list()
    for(i in 1:length(X.index_vect)){
      diff[[r]][[m]][[X.index_vect_names[i]]] <- list()
      for(j in 1:length(B.index_vect)){
        diff[[r]][[m]][[X.index_vect_names[i]]][[B.index_vect_names[j]]] <- list()
        for(k in 1:length(n.nodes_vect)){
          diff[[r]][[m]][[X.index_vect_names[i]]][[B.index_vect_names[j]]][[n.nodes_vect_names[k]]] <- matrix(0, nrow = N.batches, ncol = length(N_vect))
        }
      }
    }
  }
}

diff_m <- list()
for(r in results_names){
  diff_m[[r]] <- list()
  for(m in models[-c(1)]){
    diff_m[[r]][[m]] <- list()
    for(i in 1:length(X.index_vect)){
      diff_m[[r]][[m]][[X.index_vect_names[i]]] <- list()
      for(j in 1:length(B.index_vect)){
        diff_m[[r]][[m]][[X.index_vect_names[i]]][[B.index_vect_names[j]]] <- list()
        for(k in 1:length(n.nodes_vect)){
          diff_m[[r]][[m]][[X.index_vect_names[i]]][[B.index_vect_names[j]]][[n.nodes_vect_names[k]]] <- matrix(0, nrow = N.batches, ncol = length(N_vect))
        }
      }
    }
  }
}


for(i in 1:dim(test.options)[1]){
  
  cat(paste("\n\n\nTest: ", test.options$test.name[i], "\n", sep= ""))
  
  # Load data
  load(test.options$data.path[i])
  
  n.nodes_name <- paste("nnodes", test.options$n.nodes[i], sep = "")
  X.index_name <- paste("xi", test.options$X.index[i], sep = "")
  B.index_name <- paste("bi", test.options$B.index[i], sep = "")
  
  # Data
  Y <- data$Y
  X <- data$X_nodes
  Y_clean <- data$Y_clean
  X_clean <- data$X_clean_nodes
  B_clean <- data$B
  
  # Mesh
  FEM_basis <- data$basisobj
  mesh <- data$mesh
  mesh_data <- list(
    "nodes"    = mesh$nodes,
    "edges"    = mesh$edges,
    "elements" = mesh$triangles,
    "neigh"    = mesh$neighbors,
    "boundary" = mesh$nodesmarkers
  )
  
  # Dimensions
  S <- length(B_clean)
  
  for(j in 1:length(N_vect)){
    
    N <- N_vect[j]
    cat(paste("\n- Number of statistical units ", N, ": ", sep= ""))
    
    # Generate a grouping variable for split
    group_var <- rep(1:N, each = N, length.out =  N*N.batches)
    
    # Split the vector into N subvectors without repetition
    set.seed(0)
    subvectors <- split(sample(1:(N*N.batches)), group_var)
    
    # Multivariate
    # ||||||||||||
    
    cat(paste("\n  - Multivariate implementation:   \t", sep= ""))
    
    reference_multivariate <- list()
    for(r in results_names)
      reference_multivariate[[r]] <- list()
    
    for(k in 1:N.batches){
      
      cat(paste("|", sep= ""))
      
      # Data
      Y_batch <- matrix(Y[subvectors[[k]]], nrow = N, ncol = L)
      X_batch <- X[subvectors[[k]],]
      Y_clean_batch <- matrix(Y_clean[subvectors[[k]]], nrow = N, ncol = L)
      X_clean_batch <- X_clean[subvectors[[k]],]
      
      
      # Solve
      tic(quiet = TRUE)
      results_PLSR <- PLSR(X_batch, Y_batch, A = 3, deflation_Y = FALSE)
      elapsed <- toc(quiet = TRUE)
      
      times[["Multivariate"]][[X.index_name]][[B.index_name]][[n.nodes_name]][k,j] <- elapsed$toc - elapsed$tic
      
      reference_multivariate[["Y_hat"]][[k]] <- results_PLSR$Y_hat
      reference_multivariate[["B_hat"]][[k]] <- results_PLSR$Beta
      reference_multivariate[["X_hat"]][[k]] <- results_PLSR$X_hat
      
      error[["X_hat"]][["Multivariate"]][[X.index_name]][[B.index_name]][[n.nodes_name]][k,j] <- norm_vec(X_clean_batch - reference_multivariate[["X_hat"]][[k]])/sqrt(S*N)
      error[["Y_hat"]][["Multivariate"]][[X.index_name]][[B.index_name]][[n.nodes_name]][k,j] <- norm_vec(Y_clean_batch - reference_multivariate[["Y_hat"]][[k]])/sqrt(N)
      error[["B_hat"]][["Multivariate"]][[X.index_name]][[B.index_name]][[n.nodes_name]][k,j] <- norm_vec(B_clean - reference_multivariate[["B_hat"]][[k]])/sqrt(S)
      
      
    }
    
    # fPLSR: R
    # ||||||||
    
    reference <- list()
    for(r in results_names)
      reference[[r]] <- list()
    
    for(m in models_R){
      
      cat(paste("\n  - R implementation (", m,"):   \t", sep= ""))
    
      for(k in 1:N.batches){
        
        cat(paste("|", sep= ""))
        
        # Data
        Y_batch <- matrix(Y[subvectors[[k]]], nrow = N, ncol = L)
        X_batch <- X[subvectors[[k]],]
        Y_clean_batch <- matrix(Y_clean[subvectors[[k]]], nrow = N, ncol = L)
        X_clean_batch <- X_clean[subvectors[[k]],]
  
        # Solve
        tic(quiet = TRUE)
        results_fpls_fem <- r1fpls_fem(X_batch, Y_batch, ncomp = 3, center = TRUE,
                                       basisobj = FEM_basis, penalty = options.models_R[[m]]$lambda_s,
                                       verbose = FALSE)
        elapsed <- toc(quiet = TRUE)
        
        times[[m]][[X.index_name]][[B.index_name]][[n.nodes_name]][k,j] <- elapsed$toc - elapsed$tic
        
        TT <- as.matrix(results_fpls_fem[["TT"]])
        C <- as.matrix(results_fpls_fem[["C"]])
        X_mean <- matrix(as.matrix(results_fpls_fem[["X_mean"]]), ncol = dim(C)[1], nrow = N, byrow = TRUE)
        
        if(m == "R"){
          reference[["Y_hat"]][[k]] <- results_fpls_fem$fitted.values
          reference[["B_hat"]][[k]] <- results_fpls_fem$coefficient_function
          reference[["X_hat"]][[k]] <- TT %*% t(C) + X_mean
        }
       
        
        diff_m[["X_hat"]][[m]][[X.index_name]][[B.index_name]][[n.nodes_name]][k,j] <- norm(TT %*% t(C) + X_mean - reference_multivariate[["X_hat"]][[k]], "I")
        diff_m[["Y_hat"]][[m]][[X.index_name]][[B.index_name]][[n.nodes_name]][k,j] <- norm(results_fpls_fem$fitted.values - reference_multivariate[["Y_hat"]][[k]], "I")
        diff_m[["B_hat"]][[m]][[X.index_name]][[B.index_name]][[n.nodes_name]][k,j] <- norm(results_fpls_fem$coefficient_function - reference_multivariate[["B_hat"]][[k]], "I")
        
        error[["X_hat"]][[m]][[X.index_name]][[B.index_name]][[n.nodes_name]][k,j] <- norm_vec(X_clean_batch - (TT %*% t(C) + X_mean))/sqrt(S*N)
        error[["Y_hat"]][[m]][[X.index_name]][[B.index_name]][[n.nodes_name]][k,j] <- norm_vec(Y_clean_batch - results_fpls_fem$fitted.values)/sqrt(N)
        error[["B_hat"]][[m]][[X.index_name]][[B.index_name]][[n.nodes_name]][k,j] <- norm_vec(B_clean - results_fpls_fem$coefficient_function)/sqrt(S)
        
      }
    
    }
    
    # fPLSR: C++
    # ||||||||||
    
    for(m in models_cpp){
    
      cat(paste("\n  - C++ implementation (",m,"): \t", sep= ""))
      
      for(k in 1:N.batches){
        
        cat(paste("|", sep= ""))
        
        # Data
        Y_batch <- matrix(Y[subvectors[[k]]], nrow = N, ncol = L)
        X_batch <- X[subvectors[[k]],]
        Y_clean_batch <- matrix(Y_clean[subvectors[[k]]], nrow = N, ncol = L)
        X_clean_batch <- X_clean[subvectors[[k]],]
        
        tic(quiet = TRUE)
        # Set model
        pde <- new(Laplacian_2D_Order1, mesh_data)
        
        # Set zero forcing term
        quadrature_nodes <- pde$get_quadrature_nodes()
        f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
        pde$set_forcing_term(as.matrix(f))
        
        # Define and init model
        model <- new(FPLSR_Laplacian_2D_GeoStatNodes, pde)
        model$set_verbose(FALSE)
        model$set_lambda_s(options.models_cpp[[m]]$lambda_s)
        model$set_full_functional(options.models_cpp[[m]]$full_functional)
        model$set_smoothing3(options.models_cpp[[m]]$smoothing, options.models_cpp[[m]]$smoothing, 1e-5, 1e-9)
        
        # Set observations
        model$set_data(Y_batch, X_batch)
        
        # Solve
        model$solve()
        elapsed <- toc(quiet = TRUE)
        
        times[[m]][[X.index_name]][[B.index_name]][[n.nodes_name]][k,j] <- elapsed$toc - elapsed$tic
        
        diff[["X_hat"]][[m]][[X.index_name]][[B.index_name]][[n.nodes_name]][k,j] <- norm(model$reconstructed_field() - reference[["X_hat"]][[k]], "I")
        diff[["Y_hat"]][[m]][[X.index_name]][[B.index_name]][[n.nodes_name]][k,j] <- norm(model$fitted() - reference[["Y_hat"]][[k]], "I")
        diff[["B_hat"]][[m]][[X.index_name]][[B.index_name]][[n.nodes_name]][k,j] <- norm(model$B() - reference[["B_hat"]][[k]], "I")
        
        diff_m[["X_hat"]][[m]][[X.index_name]][[B.index_name]][[n.nodes_name]][k,j] <- norm(model$reconstructed_field() - reference_multivariate[["X_hat"]][[k]], "I")
        diff_m[["Y_hat"]][[m]][[X.index_name]][[B.index_name]][[n.nodes_name]][k,j] <- norm(model$fitted() - reference_multivariate[["Y_hat"]][[k]], "I")
        diff_m[["B_hat"]][[m]][[X.index_name]][[B.index_name]][[n.nodes_name]][k,j] <- norm(model$B() - reference_multivariate[["B_hat"]][[k]], "I")
        
        error[["X_hat"]][[m]][[X.index_name]][[B.index_name]][[n.nodes_name]][k,j] <- norm_vec(X_clean_batch - model$reconstructed_field())/sqrt(S*N)
        error[["Y_hat"]][[m]][[X.index_name]][[B.index_name]][[n.nodes_name]][k,j] <- norm_vec(Y_clean_batch - model$fitted())/sqrt(N)
        error[["B_hat"]][[m]][[X.index_name]][[B.index_name]][[n.nodes_name]][k,j] <- norm_vec(B_clean - model$B())/sqrt(S)
        
      }
      
      rm(model)
    
    }
    
  }

}


save(data.directory, results.directory, images.directory,
     N_vect, n.nodes_vect, models, X.index_vect, B.index_vect,
     n.nodes_vect_names, X.index_vect_names, B.index_vect_names,
     times, error, diff, diff_m,
     file = paste(results.directory, "results", format(Sys.time(), "_%Y%m%d_%H%M%S"), ".RData", sep = ""))

cat("\n\n\n")


# ||||||||||||||||||||||||||
# Results visualization ----
# ||||||||||||||||||||||||||

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

plot_comparison <- function(comparison, title, mains, ylab, xticks, models, colors, ylim = 0) {
  
  ylim_flag = length(ylim)>1
  
  par(mfrow = c(ceiling(length(n.nodes_vect)/2),2), mar = c(4,4,5,1), oma = c(0,0,3,0))
  
  plot.new()
  legend("center", models, col = colors, pch = 16)
  
  for(i in 1:length(n.nodes_vect)){
    
    n_test <- dim(comparison[[models[1]]][[i]])[2]
    merged <- comparison[[models[1]]][[i]]
    n_test_options_1 <- 1
    for(m in models[-1]){
      merged <- merge_tests(merged, comparison[[m]][[i]], n_test, n_test_options_1, 1)
      n_test_options_1 <- n_test_options_1 + 1
    }
    fake <- matrix(NaN, ncol = ncol(comparison[[1]][[i]]), nrow = nrow(comparison[[1]][[i]]))
    merged <- merge_tests(merged, fake, n_test, n_test_options_1, 1)
    
    if(!ylim_flag){
      if(max(merged, na.rm = TRUE) > 0 & min(merged, na.rm = TRUE) > 0)
        ylim = c(0, max(merged, na.rm = TRUE))
      else
        ylim = c(min(merged, na.rm = TRUE), max(merged, na.rm = TRUE))
    }
    t_temp <- NULL
    x_temp <- NULL
    n_models <- length(models)
    iter <- 1
    for(m in models){
      t_temp <- rbind(t_temp, 0:(n_test-1)*(n_models+1)+iter)
      x_temp <- rbind(x_temp, colMeans(comparison[[m]][[i]]))
      iter <- iter + 1
    }
    matplot(t(t_temp), 
            t(x_temp),
            main = mains[i],
            col = colors, type = "o", lwd = 2, pch = 19, lty = 1,
            xlab = "N", xaxt = "n", xlim = c(min(t_temp)-0.3, max(t_temp)+0.3),
            ylab = ylab, ylim = ylim)
    boxplot(merged, col = c(colors, "white"), add = TRUE, xaxt = "n")
    
    axis(1, at=(1:(n_test))*(n_models+1)-(n_models)/2, labels=xticks)
    grid(nx = NA, ny = NULL)
    abline(v = 0:(n_test)*(n_models+1), col = "lightgrey", lty = 3)
    
  }
  
  title(title, outer = TRUE, line = 1, cex.main = 1.5)
  
}


cat("\n")
cat("\n# ||||||||||||||||||||||||||")
cat("\n# Results visualization ----")
cat("\n# ||||||||||||||||||||||||||\n")

if (!file.exists(images.directory)){
  dir.create(images.directory)
}

list.files(results.directory)
load(paste(results.directory, tail(list.files(results.directory), n = 1), sep = ""))


# Choose a color palette from RColorBrewer
palette_name <- "Paired"

# Generate N different colors from the chosen palette
colors <- c(rev(brewer.pal(3, "Reds"))[1],
            rev(brewer.pal(3, "Greens"))[1:2],
            rev(brewer.pal(3, "Blues")))


## Execution time ----
## |||||||||||||||||||

# Pictures settings
width = 2000
height_per_row = 650
pointsize = 40

selected <- c(2,4)
speedup <- list()
speedup[["R vs C++"]] <- list()
for(i in X.index_vect_names){
  for(j in B.index_vect_names){
    comparison <- list()
    for(m in models[selected]){
      comparison[[m]] <- times[[m]][[i]][[j]]
    }
    jpeg(paste(images.directory, "time_",i,"_",j, ".jpg", sep = ""),
         width = width, height = height_per_row*ceiling(length(n.nodes_vect)/2), units = "px", quality = 100, pointsize = pointsize,
         type = "windows")
    plot_comparison(comparison,
                    "Execution time",
                    paste("K =", n.nodes_vect^2),
                    "[s]",
                    N_vect, models[selected], colors[selected])
    dev.off()
    
    # Speed up
    for(k in n.nodes_vect_names){
      speedup[["R vs C++"]][[k]] <- (comparison[["R"]][[k]] - comparison[["C++ ff"]][[k]])/comparison[["R"]][[k]]*100
    }
    jpeg(paste(images.directory, "speedup_",i,"_",j, ".jpg", sep = ""),
         width = width, height = height_per_row*ceiling(length(n.nodes_vect)/2), units = "px", quality = 100, pointsize = pointsize,
         type = "windows")
    plot_comparison(speedup,
                    "Speedup",
                    paste("K =", n.nodes_vect^2),
                    "[%]",
                    N_vect, c("R vs C++"), colors[3], c(0,100))
    dev.off()
    
  }
}


## Difference with R implementation ----
## |||||||||||||||||||||||||||||||||||||


# Pictures settings
width = 2000
height_per_row = 800
pointsize = 40

colors_errors <- rev(brewer.pal(3, "Reds"))
for(i in X.index_vect_names){
  for(j in B.index_vect_names){
    comparison <- list()
    for(r in results_names)
      comparison[[r]] <- lapply(diff[[r]][["C++ ff"]][[i]][[j]], log10)

    jpeg(paste(images.directory, "difference_",i,"_",j, ".jpg", sep = ""),
         width = width, height = height_per_row*ceiling(length(n.nodes_vect)/2), units = "px", quality = 100, pointsize = pointsize,
         type = "windows")
    plot_comparison(comparison,
                    "R vs C++",
                    paste("K =", n.nodes_vect^2),
                    "log10(infinity norm)",
                    N_vect, results_names, colors_errors)
    dev.off()
  }
}


## Consistency with multivariate implementation ----
## |||||||||||||||||||||||||||||||||||||||||||||||||

# Pictures settings
width = 2000
height_per_row = 650
pointsize = 40

for(m in models[c(3,5)]){
  colors_errors <- rev(brewer.pal(3, "Reds"))
  for(i in X.index_vect_names){
    for(j in B.index_vect_names){
      comparison <- list()
      for(r in results_names)
        comparison[[r]] <- lapply(diff_m[[r]][[m]][[i]][[j]], log10)
      
      jpeg(paste(images.directory, "difference_m_",m,"_",i,"_",j, ".jpg", sep = ""),
           width = width, height = height_per_row*ceiling(length(n.nodes_vect)/2), units = "px", quality = 100, pointsize = pointsize,
           type = "windows")
      plot_comparison(comparison,
                      paste(strsplit(m, " ")[[1]][1]," vs Multivariate",
                            sep = ""),
                      paste("K =", n.nodes_vect^2),
                      "log10(infinity norm)",
                      N_vect, results_names, colors_errors)
      dev.off()
    }
  }
}


## Performances ----
## |||||||||||||||||

# Pictures settings
width = 2000
height_per_row = 1000
pointsize = 40

colors_errors <- rev(brewer.pal(3, "Reds"))
chosen <- c(1,2,6)
for(r in results_names){
  for(i in X.index_vect_names){
    for(j in B.index_vect_names){
      comparison <- list()
      for(m in models[chosen]){
        comparison[[m]] <- error[[r]][[m]][[i]][[j]]
      }
      
      jpeg(paste(images.directory, "error_",r,"_",i,"_",j, ".jpg", sep = ""),
           width = width, height = height_per_row*ceiling(length(n.nodes_vect)/2), units = "px", quality = 100, pointsize = pointsize,
           type = "windows")
      plot_comparison(comparison,
                      paste("RMSE", r,
                            sep = " "),
                      paste("K =", n.nodes_vect^2),
                      "",
                      N_vect, models[chosen], colors[chosen])
      dev.off()
    }
  }
}

colors_errors <- rev(brewer.pal(3, "Reds"))
chosen <- c(1,6)
for(r in c("B_hat")){
  for(i in X.index_vect_names){
    for(j in B.index_vect_names){
      comparison <- list()
      for(m in models[chosen]){
        comparison[[m]] <- error[[r]][[m]][[i]][[j]]
      }
      
      jpeg(paste(images.directory, "error_only_",r,"_",i,"_",j, ".jpg", sep = ""),
           width = width, height = height_per_row*ceiling(length(n.nodes_vect)/2), units = "px", quality = 100, pointsize = pointsize,
           type = "windows")
      plot_comparison(comparison,
                      paste("RMSE", r,
                            sep = " "),
                      paste("K =", n.nodes_vect^2),
                      "",
                      N_vect, models[chosen], colors[chosen])
      dev.off()
    }
  }
}




## Data Images ----
## ||||||||||||||||

# images.directory <- paste(images.directory, "data/", sep = "")
# if (!file.exists(images.directory)){
#   dir.create(images.directory)
# }
# 
# mesh <- data$mesh
# nodes <- mesh$nodes
# 
# # Pictures settings
# width = 800
# height = 900
# pointsize = 30
# 
# for(i in 1:2){
#   x <- generate_X(nodes, i)
#   jpeg(paste(images.directory, "X_us_", i, ".jpg", sep = ""),
#        width = width, height = height, units = "px", quality = 100, pointsize = pointsize,
#        type = "windows")
#   plot_field(nodes, nodes, x, range(x), "", TRUE)
#   dev.off()
# }
# 
# for(i in 1:5){
#   B <- generate_B(nodes, i)
#   jpeg(paste(images.directory, "B_us_", i, ".jpg", sep = ""),
#        width = width, height = height, units = "px", quality = 100, pointsize = pointsize,
#        type = "windows")
#   plot_field(nodes, nodes, B, range(B), "", TRUE)
#   dev.off()
# }
