# %%%%%%%%%%%%%%%%%%%% #
# %% PLS vs SIMPLS %% #
# %%%%%%%%%%%%%%%%%%% #

rm(list = ls())
graphics.off()


# ||||||||||||||
# Functions ----
# ||||||||||||||

load("PLSR.RData")
load("SIMPLS.RData")


plot_prediction <- function(Y_hat, Y, Y.labs) {
  
  N <- dim(Y)[1]
  L <- dim(Y)[2]
  par(mfrow = c(L,1))
  
  for(l in 1:L){
    Y_hat.range <- range(Y_hat[,l])
    Y.range <- range(Y[,l])
    plot(1:N, Y_hat[,l], pch = 8, main = Y.labs[l],
         xlab = "Realizations", ylab = Y.labs[l], xaxt = "n",
         ylim = c(min(Y_hat.range[1], Y.range[1]), max(Y_hat.range[2], Y.range[2])))
    points(1:N, Y[,l])
    axis(1, at = 1:N, labels = 1:N)
    grid()
  }
}

plot_beta <- function(Beta){
  S <- dim(Beta)[1]
  L <- dim(Beta)[2]
  par(mfrow = c(l,1))
  for(l in 1:L){
    image(matrix(as.numeric(Beta[,l]), sqrt(S), sqrt(S)), main = paste("Beta_",l, sep = '' ))
  }
}

compare_results <- function(plsr_1, plsr_2){
  
  elements_to_compare <- c("W", "V", "TT", "UU", "C", "D", "Beta", "Y_hat", "X_hat")
  
  for(elem in elements_to_compare){
    
    norm <- norm(plsr_1[[elem]] - plsr_2[[elem]], "I")
    if(norm < 1e-9)
      print(paste("Test on", elem, "passed!", sep = ' '))
    else
      print(paste("Test on", elem, "failed!", sep = ' '))
    print(paste("Norm:", norm, sep = ''))
    
  }
  
}


# |||||||||
# Data ----
# |||||||||

path1 <- "../../fdaPDE/test/data/models/FPLSR/2D_test0/test1/"
path2 <- "../../fdaPDE/test/data/models/FPLSR/2D_test0/test2/"
path3 <- "../../fdaPDE/test/data/models/FPLSR/2D_test0/test3/"
path4 <- "../../fdaPDE/test/data/models/FPLSR/2D_test0/test4/"
path5 <- "../../fdaPDE/test/data/models/FPLSR/2D_test0/test5/"
path6 <- "../../fdaPDE/test/data/models/FPLSR/2D_test0/test6/"

X = read.csv(paste(path1, "X.csv", sep = ''), header = TRUE)[,2:3601]

Y1 = matrix(read.csv(paste(path1, "Y.csv", sep = ''), header = TRUE)[,2], ncol = 1)
Y2 = matrix(read.csv(paste(path2, "Y.csv", sep = ''), header = TRUE)[,2], ncol = 1)
Y3 = matrix(read.csv(paste(path2, "Y.csv", sep = ''), header = TRUE)[,2], ncol = 1)
Y4 = matrix(read.csv(paste(path2, "Y.csv", sep = ''), header = TRUE)[,2], ncol = 1)
Y5 = matrix(read.csv(paste(path2, "Y.csv", sep = ''), header = TRUE)[,2], ncol = 1)
Y6 = matrix(read.csv(paste(path2, "Y.csv", sep = ''), header = TRUE)[,2], ncol = 1)

Y = cbind(Y1)
Y.labs = c("Y1")

## Dimensions ----
## |||||||||||||||

N <- dim(X)[1]
S <- dim(X)[2]
L <- dim(Y)[2]


## De-trending ----
## |||||||||||||||

Xc = scale(X, scale = FALSE)
X.mean = attr(Xc, "scaled:center")

Yc = scale(Y, scale = FALSE)
Y.mean = attr(Yc, "scaled:center")

# Mean values
Y.MEAN = matrix(Y.mean, ncol = L, nrow = N, byrow = TRUE)
X.MEAN = matrix(X.mean, ncol = S, nrow = N, byrow = TRUE)


## Plot ----
## |||||||||

# Y-spaces
par(mfrow = c(L,1))
for(l in 1:L){
  plot(1:N, Y[,l], main = Y.labs[l])
  abline(h = Y.mean[l], lty = 2, col = "black")
  grid()
}
# X-space
par(mfrow = c(1,1))
image(matrix(as.numeric(X[1,]), 60, 60), main = "X")


# |||||||||
# PLSR ----
# |||||||||

## Standard implementation ----
## ||||||||||||||||||||||||||||

plsr = PLSR(X, Y, 3)

Beta = plsr[["Beta"]]
Y_hat = plsr[["Y_hat"]]
X_hat = plsr[["X_hat"]]

plot_beta(Beta)
par(mfrow = c(1,1))
image(matrix(as.numeric(X_hat[1,]), 60, 60), main = "X_hat")
plot_prediction(Y_hat, Y, Y.labs)


## No Y deflation ----
## |||||||||||||||||||

plsr_noYDefl = PLSR(X, Y, 3, TRUE)

TT_noYDefl = plsr_noYDefl[["TT"]]
B_noYDefl = plsr_noYDefl[["B"]]

Beta_noYDefl = plsr_noYDefl[["Beta"]]
Y_hat_noYDefl = plsr_noYDefl[["Y_hat"]]
X_hat_noYDefl = plsr_noYDefl[["X_hat"]]

par(mfrow = c(1,1))
plot_beta(Beta_noYDefl)
image(matrix(as.numeric(X_hat[1,]), 60, 60), main = "X_hat")
plot_prediction(Y_hat, Y, Y.labs)

Y_hat_noYDefl_2 = TT_noYDefl %*% B_noYDefl %*% t(TT_noYDefl) %*% Yc + Y.MEAN
norm(Y_hat_noYDefl - Y_hat_noYDefl_2, "I")


# SIMPLS ----
# |||||||||||

simpls = SIMPLS(X, Y, 5)

TT_simpls = simpls[["TT"]]

Beta_simpls = simpls[["Beta"]]
Y_hat_simpls = simpls[["Y_hat"]]
X_hat_simpls = simpls[["X_hat"]]

h_simpls = simpls[["h"]]
varX_simpls = simpls[["varX"]]
varY_simpls = simpls[["varY"]]

par(mfrow = c(1,1))
plot_beta(Beta_simpls)
image(matrix(as.numeric(X_hat_simpls[1,]), 60, 60), main = "X_hat")
plot_prediction(Y_hat_simpls, Y, Y.labs)

par(mfrow = c(1, 1))
barplot(h_simpls, main = "Leverages", xlab = "Samples")

par(mfrow = c(1, 2))
plot(cumsum(varY_simpls), type = "b",
     main = "Explained Y variance",
     xlab = "# of PLS components",
     ylab = "Var(Y)")
abline(h = var(Y), col = "blue", lty = 2)
grid()
plot(cumsum(varX_simpls), type = "b",
     main = "Explained X variance",
     xlab = "# of PLS components",
     ylab = "Var(X)")
grid()

Y_hat_simpls_2 = TT_simpls %*% t(TT_simpls) %*% Yc + Y.MEAN
norm(Y_hat_simpls - Y_hat_simpls_2, "I")




# Comparison ----
# |||||||||||||||

## Univariate response ----
## ||||||||||||||||||||||||

Y = cbind(Y1)

plsr = PLSR(X, Y, 3)
plsr_noYDefl = PLSR(X, Y, 3, TRUE)
simpls_comparison = SIMPLS(X, Y, 3)

compare_results(plsr, plsr_noYDefl)
compare_results(plsr, simpls_comparison)


## Multivariate response ----
## ||||||||||||||||||||||||||

Y = cbind(Y1, Y2)

plsr = PLSR(X, Y, 3)
plsr_noYDefl = PLSR(X, Y, 3, TRUE)
simpls_comparison = SIMPLS(X, Y, 3)

compare_results(plsr, plsr_noYDefl)
compare_results(plsr, simpls_comparison)

# In the multivariate case PLSR and SIMPLS also give different final result.
# More tests are needed to enstablish which is the best
