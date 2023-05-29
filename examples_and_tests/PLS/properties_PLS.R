# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# %% Multivariate Principal Least Squares Regression %% #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

rm(list = ls())
graphics.off()


# ||||||||||||||
# Libraries ----
# ||||||||||||||



# ||||||||||||||
# Functions ----
# ||||||||||||||

load("PLSR.RData")
# load("PLSR_star.RData")


plot_prediction <- function(Y_hat, Y, N, Y.labs) {
  par(mfrow = c(2,1))
  Y_hat.x_range <- range(Y_hat[,1])
  Y.x_range <- range(Y[,1])
  plot(1:N, Y_hat[,1], pch = 8, main = Y.labs[1],
       xlab = "Realizations", ylab = Y.labs[1], xaxt = "n",
       ylim = c(min(Y_hat.x_range[1], Y.x_range[1]), max(Y_hat.x_range[2], Y.x_range[2])))
  points(1:N, Y[,1])
  axis(1, at = 1:N, labels = 1:N)
  grid()
  Y_hat.y_range <- range(Y_hat[,2])
  Y.y_range <- range(Y[,2])
  plot(1:N, Y_hat[,2], pch = 8, main = Y.labs[2],
       xlab = "Realizations", ylab = Y.labs[2], xaxt = "n",
       ylim = c(min(Y_hat.y_range[1], Y.y_range[1]), max(Y_hat.y_range[2], Y.y_range[2])))
  points(1:N, Y[,2])
  axis(1, at = 1:N, labels = 1:N)
  grid()
}


# |||||||||
# Data ----
# |||||||||

## Covariates ----
## |||||||||||||||

height <- c(195, 185, 156, 181, 164, 158, 185, 186, 187)
distance_to_work <- c(45, 30, 90, 45, 30, 0, 15, 105, 45)
X <- data.frame(height, distance_to_work)
X.labs <- c("Height", "Distance to work")

## Responses ----
## ||||||||||||||

weight <- c(95,73,63, 80, 70, 55, 89, 78, 91)
shoe_size <- c(44, 43, 36, 45, 40, 38, 44, 45, 43)
Y <- data.frame(weight, shoe_size)
Y.labs <- c("Weight", "Shoe size")


## Dimensions ----
## |||||||||||||||

N <- dim(X)[1]
S <- dim(X)[2]
L <- dim(Y)[2]


## Detrending ----
## |||||||||||||||

Xc = scale(X, scale = FALSE)
X.mean = attr(Xc, "scaled:center")
X.MEAN <- matrix(X.mean, ncol = S, nrow = N, byrow = TRUE)

Yc = scale(Y, scale = FALSE)
Y.mean = attr(Yc, "scaled:center")
Y.MEAN <- matrix(Y.mean, ncol = L, nrow = N, byrow = TRUE)


## Plot ----
## |||||||||

par(mfrow = c(1,2))

# Y-space
plot(Y, main = "Y", xlab = Y.labs[1], ylab = Y.labs[2], asp = 1)
abline(v = Y.mean[1], h = Y.mean[2], lty = 2, col = "black")
grid()

# X-space
plot(X, main = "X", xlab = X.labs[1], ylab = X.labs[2], asp = 1)
abline(v = X.mean[1], h = X.mean[2], lty = 2, col = "black",  lwd = 1)
grid()


## Plot centered ----
## ||||||||||||||||||

par(mfrow = c(1,2))

# Y-space
plot(Yc, main = "Y centered", xlab = Y.labs[1], ylab = Y.labs[2], asp = 1)
abline(v = 0, h = 0, lty = 2, col = "black")
grid()

# X-space
plot(Xc, main = "X centered", xlab = X.labs[1], ylab = X.labs[2], asp = 1)
abline(v =0, h = 0, lty = 2, col = "black",  lwd = 1)
grid()


# |||||||||||||||||||||||||||||||||||||||||
# Principal Component Regression - PCR ----
# |||||||||||||||||||||||||||||||||||||||||

## X-space PCA ----
## ||||||||||||||||

pca.X <- princomp(X)
X.loadings <- pca.X$loadings[,]
X.scores <- pca.X$scores[,]

## Y-space PCA ----
## ||||||||||||||||

pca.Y <- princomp(Y)
Y.loadings <- pca.Y$loadings[,]
Y.scores <- pca.Y$scores[,]

## Plot ----
## |||||||||

par(mfrow = c(2,2))

# Y-space
plot(Y, main = "Y", xlab = Y.labs[1], ylab = Y.labs[2], asp = 1)
abline(v = Y.mean[1], h = Y.mean[2], lty = 2, col = "grey")
abline(Y.mean[2] - Y.loadings[2,1]/Y.loadings[1,1]*Y.mean[1], Y.loadings[2,1]/Y.loadings[1,1], col = "blue", lwd = 2)
grid()

# X-scores vs Y-scores
plot(X.scores[,1], Y.scores[,1], main = "1st Principal Components", xlab = "X-scores", ylab = "Y-scores", pch = 8, asp = 1)
b1 = t(Y.scores[,1]) %*% X.scores[,1] / sum( X.scores[,2]^2)
abline(0, b1, col = "red", lwd = 2)
grid()

cov(Y.scores[,1], X.scores[,1])/sqrt(var(Y.scores[,1])*var(X.scores[,1]))

# Empty plot 
plot.new()

# X-space
plot(X, main = "X", xlab = X.labs[1], ylab = X.labs[2], asp = 1)
abline(v = X.mean[1], h = X.mean[2], lty = 2, col = "grey")
abline(X.mean[2] - X.loadings[2,1]/X.loadings[1,1]*X.mean[1], X.loadings[2,1]/X.loadings[1,1], col = "blue", lwd = 2)
grid()


# ||||||||||||||||||||||||||||||||||||||||||||
# Partial Least Squares Regression - PLSR ----
# ||||||||||||||||||||||||||||||||||||||||||||

plsr <- PLSR(X, Y, 2, TRUE)

W = plsr[["W"]]
V = plsr[["V"]]
TT = plsr[["TT"]]
UU = plsr[["UU"]]
B = plsr[["B"]]
C = plsr[["C"]]
D = plsr[["D"]]
EE = plsr[["EE"]]
FF = plsr[["FF"]]
XX = plsr[["XX"]]
YY = plsr[["YY"]]


## Properties ----
## |||||||||||||||

### Comparison: W, C ----
### |||||||||||||||||||||

# Comparison w1 vs c1
t(W[,1]) %*% C[,1]
sqrt(sum(W[,1]^2))
sqrt(sum(C[,1]^2))
# Comparison w2 vs c2
t(W[,2]) %*% C[,2]
sqrt(sum(W[,1]^2))
sqrt(sum(C[,2]^2))


### Comparison: D, V ----
### |||||||||||||||||||||

D
V %*% B


### Orthogonality of W ----
### |||||||||||||||||||||||

t(W[,1]) %*% W[,2] # YES


### Orthogonality of V ----
### |||||||||||||||||||||||

t(V[,1]) %*% V[,2] # NO


### Orthogonality of T ----
### |||||||||||||||||||||||

t(TT[,1]) %*% TT[,2]


### Orthogonality between W anc C ----
### ||||||||||||||||||||||||||||||||||

t(W[,1]) %*% C[,2]


### Orthogonality between W and E ----
### ||||||||||||||||||||||||||||||||||
 
t(W[,1]) %*% t(EE)
t(W[,2]) %*% t(EE)


### Orthogonality between T and E ----
### ||||||||||||||||||||||||||||||||||

t(TT[,1]) %*% EE
t(TT[,2]) %*% EE


## Plot, 1st latent component ---- 
## |||||||||||||||||||||||||||||||

par(mfrow = c(2,2))

# Y-space
plot(YY[[1]] + Y.MEAN, main = "Y", xlab = Y.labs[1], ylab = Y.labs[2], asp = 1)
abline(v = Y.mean[1], h = Y.mean[2], lty = 2, col = "grey")
abline(Y.mean[2] - Y.loadings[2,1]/Y.loadings[1,1]*Y.mean[1], Y.loadings[2,1]/Y.loadings[1,1], col = "blue", lwd = 1, lty = 2)
abline(Y.mean[2] - V[2,1]/V[1,1]*Y.mean[1], V[2,1]/V[1,1], col = "green", lwd = 2)
grid()

# X-scores vs Y-scores
plot(TT[,1], UU[,1], main = "1st Latent components", xlab = "t", ylab = "u", pch = 8, asp = 1)
abline(0, B[1,1], col = "red", lwd = 2)
grid()

cov(UU[,1], TT[,1])/sqrt(var(UU[,1])*var(TT[,1]))

# Void plot 
plot.new()

# X-space
plot(XX[[1]] + X.MEAN, main = "X", xlab = X.labs[1], ylab = X.labs[2], asp = 1)
abline(v = X.mean[1], h = X.mean[2], lty = 2, col = "grey")
abline(X.mean[2] - X.loadings[2,1]/X.loadings[1,1]*X.mean[1], X.loadings[2,1]/X.loadings[1,1], col = "blue", lwd = 1, lty = 2)
abline(X.mean[2] - W[2,1]/W[1,1]*X.mean[1], W[2,1]/W[1,1], col = "green", lwd = 2)
grid()


## Plot, 2nd latent component ---- 
## |||||||||||||||||||||||||||||||

par(mfrow = c(2,2))

# Y-space
plot(YY[[2]] + Y.MEAN, main = "Y - deflated", xlab = Y.labs[1], ylab = Y.labs[2], asp = 1)
abline(v = Y.mean[1], h = Y.mean[2], lty = 2, col = "grey")
abline(Y.mean[2] - V[2,2]/V[1,2]*Y.mean[1], V[2,2]/V[1,2], col = "green", lwd = 2)
grid()

# X-scores vs Y-scores
plot(TT[,2], UU[,2], main = "2nd Latent components", xlab = "t", ylab = "u", pch = 8, asp = 1)
abline(0, B[2,2], col = "red", lwd = 2)
grid()

cov(UU[,2], TT[,2])/sqrt(var(UU[,2])*var(TT[,2]))

# Void plot 
plot.new()

# X-space
plot(XX[[2]] + X.MEAN, main = "X - deflated", xlab = X.labs[1], ylab = X.labs[2], asp = 1)
abline(v = X.mean[1], h = X.mean[2], lty = 2, col = "grey")
abline(X.mean[2] - W[2,2]/W[1,2]*X.mean[1], W[2,2]/W[1,2], col = "green", lwd = 2)
grid()


## Deflation ----
## ||||||||||||||

par(mfrow = c(3,2))

Y.x_range = range(YY[[1]][,1])
Y.y_range = range(YY[[1]][,2])

X.x_range = range(XX[[1]][,1])
X.y_range = range(XX[[1]][,2])

# Y-space, initial
plot(YY[[1]], main = "Y, initial", xlab = Y.labs[1], ylab = Y.labs[2], pch = 8,
     xlim = Y.x_range, ylim = Y.y_range, asp=1)
abline(0, V[2,1]/V[1,1], col = "green", lwd = 2)
grid()

# X-space, initial
plot(XX[[1]], main = "X, initial", xlab = X.labs[1], ylab = X.labs[2], pch = 8,
     xlim = X.x_range, ylim = X.y_range, asp=1)
abline(0, W[2,1]/W[1,1], col = "green", lwd = 2)
grid()

# Y-space, 1st deflation
plot(YY[[2]], main = "Y, 1st deflation", xlab = Y.labs[1], ylab = Y.labs[2], pch = 8,
     xlim = Y.x_range, ylim = Y.y_range, asp=1)
abline(0, V[2,1]/V[1,1], col = "green", lwd = 1, lty = 2)
abline(0, V[2,2]/V[1,2], col = "green", lwd = 2)
grid()

# X-space, 1st deflation
plot(XX[[2]], main = "X, 1st deflation", xlab = X.labs[1], ylab = X.labs[2], pch = 8,
     xlim = X.x_range, ylim = X.y_range, asp=1)
abline(0, W[2,1]/W[1,1], col = "green", lwd = 1, lty = 2)
abline(0, W[2,2]/W[1,2], col = "green", lwd = 2)
grid()

# Y-space, 2nd deflation
plot(YY[[3]], main = "Y, 2nd deflation", xlab = Y.labs[1], ylab = Y.labs[2], pch = 8,
     xlim = Y.x_range, ylim = Y.y_range, asp=1)
abline(0, V[2,2]/V[1,2], col = "green", lwd = 1, lty = 2)
grid()

# X-space, 2nd deflation
plot(XX[[3]], main = "X, 2nd deflation", xlab = X.labs[1], ylab = X.labs[2], pch = 8,
     xlim = X.x_range, ylim = X.y_range, asp=1)
abline(0, W[2,2]/W[1,2], col = "green", lwd = 1, lty = 2)
grid()


## Prediction ----
## |||||||||||||||


### PLS decomposition ----
### |||||||||||||||||||||||

Y_h <- TT %*% t(D) + Y.MEAN
Y_dec <- TT %*% t(D) + Y.MEAN + YY[[3]] 

# Exact
par(mfrow = c(1,1))
plot(Y, main = "Y", col = rainbow(N), asp = 1,
     xlab = Y.labs[1], ylab = Y.labs[2])
points(Y_dec , pch = 8, col = rainbow(N))
grid()

# Approximated
plot_prediction(Y_h, Y, N, Y.labs)

sum((Y-Y_h)^2)/N


### PLS (Y-space) decomposition ----
### ||||||||||||||||||||||||||||||||

Y_hat_1 = UU %*% t(V) + Y.MEAN

plot_prediction(Y_hat_1, Y, N, Y.labs)

sum((Y-Y_hat_1)^2)/N

### PLS regression ----
### |||||||||||||||||||

Beta <- W %*% solve(t(C) %*% W, B %*% t(V))
Y_hat_2 <- Xc %*% W %*% solve(t(C) %*% W, B %*% t(V)) + Y.MEAN

plot_prediction(Y_hat_2, Y, N, Y.labs)

sum((Y-Y_hat_2)^2)/N

# # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# # Partial Least Squares Regression no deflation - PLSR* ----
# # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# 
# plsr_star <- PLSR_star(Xc, Yc, 2)
# 
# W_star = plsr_star[["W"]]
# V_star = plsr_star[["V"]]
# TT_star = plsr_star[["TT"]]
# UU_star = plsr_star[["UU"]]
# B_star = plsr_star[["B"]]
# C_star = plsr_star[["C"]]
# D_star = plsr_star[["D"]]
# EE_star = plsr_star[["EE"]]
# FF_star = plsr_star[["FF"]]
# XX_star = plsr_star[["XX"]]
# YY_star = plsr_star[["YY"]]
# 
# 
# ## Properties ----
# ## |||||||||||||||
# 
# ### Comparison: W_star, C_star ----
# ### |||||||||||||||||||||||||||||||
# 
# # Comparison w1 vs c1
# t(W_star[,1]) %*% C_star[,1]
# sqrt(sum(W_star[,1]^2))
# sqrt(sum(C_star[,1]^2))
# # Comparison w2 vs c2
# t(W_star[,2]) %*% C_star[,2]
# sqrt(sum(W_star[,2]^2))
# sqrt(sum(C_star[,2]^2))
# 
# 
# ### Comparison: D_star, V_star ----
# ### ||||||||||||||||||||||||||||||||
# 
# D_star
# V_star %*% B_star
# 
# 
# ### Orthogonality of W_star ----
# ### ||||||||||||||||||||||||||||
# 
# t(W_star[,1]) %*% W_star[,2]
# 
# 
# ### Orthogonality of W_star ----
# ### ||||||||||||||||||||||||||||
# 
# t(V_star[,1]) %*% V_star[,2]
# 
# 
# ### Orthogonality of T_star ----
# ### ||||||||||||||||||||||||||||
# 
# t(TT_star[,1]) %*% TT_star[,2]
# 
# 
# ### Orthogonality between W_star anc C_star ----
# ### ||||||||||||||||||||||||||||||||||||||||||||
# 
# t(W_star[,1]) %*% C_star[,2]
# 
# 
# ### Orthogonality between W_star and E_star ----
# ### ||||||||||||||||||||||||||||||||||||||||||||
# 
# t(W_star[,1]) %*% t(EE_star)
# t(W_star[,2]) %*% t(EE_star)
# 
# 
# ### Orthogonality between T_star and E_star ----
# ### ||||||||||||||||||||||||||||||||||||||||||||
# 
# t(TT_star[,1]) %*% EE_star
# t(TT_star[,2]) %*% EE_star
# 
# ## Plot, 1st latent component ----
# ## |||||||||||||||||||||||||||||||
# 
# par(mfrow = c(2,2))
# 
# # Y-space
# plot(YY_star[[1]] + Y.MEAN, main = "Y", xlab = Y.labs[1], ylab = Y.labs[2], asp = 1)
# abline(v = Y.mean[1], h = Y.mean[2], lty = 2, col = "grey")
# abline(Y.mean[2] - Y.loadings[2,1]/Y.loadings[1,1]*Y.mean[1], Y.loadings[2,1]/Y.loadings[1,1], col = "blue", lwd = 1, lty = 2)
# abline(Y.mean[2] - V_star[2,1]/V_star[1,1]*Y.mean[1], V_star[2,1]/V_star[1,1], col = "green", lwd = 2)
# grid()
# 
# # X-scores vs Y-scores
# plot(TT_star[,1], UU_star[,1], main = "1st Latent components", xlab = "t", ylab = "u", pch = 8, asp = 1)
# abline(0, B_star[1,1], col = "red", lwd = 2)
# grid()
# 
# cov(UU_star[,1], TT_star[,1])/sqrt(var(UU_star[,1])*var(TT_star[,1]))
# 
# # Void plot 
# plot.new()
# 
# # X-space
# plot(XX_star[[1]] + X.MEAN, main = "X", xlab = X.labs[1], ylab = X.labs[2], asp = 1)
# abline(v = X.mean[1], h = X.mean[2], lty = 2, col = "grey")
# abline(X.mean[2] - X.loadings[2,1]/X.loadings[1,1]*X.mean[1], X.loadings[2,1]/X.loadings[1,1], col = "blue", lwd = 1, lty = 2)
# abline(X.mean[2] - W_star[2,1]/W_star[1,1]*X.mean[1], W_star[2,1]/W_star[1,1], col = "green", lwd = 2)
# grid()
# 
# ## Plot, 2nd latent component ----
# ## |||||||||||||||||||||||||||||||
# 
# par(mfrow = c(2,2))
# 
# # Y-space
# plot(YY_star[[2]] + Y.MEAN, main = "Y", xlab = Y.labs[1], ylab = Y.labs[2], asp = 1)
# abline(v = Y.mean[1], h = Y.mean[2], lty = 2, col = "grey")
# abline(Y.mean[2] - V_star[2,2]/V_star[1,2]*Y.mean[1], V_star[2,2]/V_star[1,2], col = "green", lwd = 2)
# grid()
# 
# # X-scores vs Y-scores
# plot(TT_star[,2], UU_star[,2], main = "2nd Latent components", xlab = "t", ylab = "u", pch = 8, asp = 1)
# abline(0, B_star[2,2], col = "red", lwd = 2)
# grid()
# 
# cov(UU_star[,2], TT_star[,2])/sqrt(var(UU_star[,2])*var(TT_star[,2]))
# 
# # V_staroid plot 
# plot.new()
# 
# # X-space
# plot(XX_star[[2]] + X.MEAN, main = "X", xlab = X.labs[1], ylab = X.labs[2], asp = 1)
# abline(v = X.mean[1], h = X.mean[2], lty = 2, col = "grey")
# abline(X.mean[2] - W_star[2,2]/W_star[1,2]*X.mean[1], W_star[2,2]/W_star[1,2], col = "green", lwd = 2)
# grid()
# 
# 
# ## Deflation ----
# ## ||||||||||||||
# 
# par(mfrow = c(3,2))
# 
# Y.x_range = range(YY_star[[1]][,1])
# Y.y_range = range(YY_star[[1]][,2])
# 
# X.x_range = range(XX_star[[1]][,1])
# X.y_range = range(XX_star[[1]][,2])
# 
# # Y-space, initial
# plot(YY_star[[1]], main = "Y, initial", xlab = Y.labs[1], ylab = Y.labs[2], pch = 8,
#      xlim = Y.x_range, ylim = Y.y_range, asp=1)
# abline(0, V_star[2,1]/V_star[1,1], col = "green", lwd = 2)
# grid()
# 
# # X-space, initial
# plot(XX_star[[1]], main = "X, initial", xlab = X.labs[1], ylab = X.labs[2], pch = 8,
#      xlim = X.x_range, ylim = X.y_range, asp=1)
# abline(0, W_star[2,1]/W_star[1,1], col = "green", lwd = 2)
# grid()
# 
# # Y-space, 1st deflation
# plot(YY_star[[2]], main = "Y, 1st deflation", xlab = Y.labs[1], ylab = Y.labs[2], pch = 8,
#      xlim = Y.x_range, ylim = Y.y_range, asp=1)
# abline(0, V_star[2,1]/V_star[1,1], col = "green", lwd = 1, lty = 2)
# abline(0, V_star[2,2]/V_star[1,2], col = "green", lwd = 2)
# grid()
# 
# # X-space, 1st deflation
# plot(XX_star[[2]], main = "X, 1st deflation", xlab = X.labs[1], ylab = X.labs[2], pch = 8,
#      xlim = X.x_range, ylim = X.y_range, asp=1)
# abline(0, W_star[2,1]/W_star[1,1], col = "green", lwd = 1, lty = 2)
# abline(0, W_star[2,2]/W_star[1,2], col = "green", lwd = 2)
# grid()
# 
# # Y-space, 2nd deflation
# plot(YY_star[[3]], main = "Y, 2nd deflation", xlab = Y.labs[1], ylab = Y.labs[2], pch = 8,
#      xlim = Y.x_range, ylim = Y.y_range, asp=1)
# abline(0, V_star[2,2]/V_star[1,2], col = "green", lwd = 1, lty = 2)
# grid()
# 
# # X-space, 2nd deflation
# plot(XX_star[[3]], main = "X, 2nd deflation", xlab = X.labs[1], ylab = X.labs[2], pch = 8,
#      xlim = X.x_range, ylim = X.y_range, asp=1)
# abline(0, W_star[2,2]/W_star[1,2], col = "green", lwd = 1, lty = 2)
# grid()
# 
# 
# ## Prediction ----
# ## |||||||||||||||
# 
# 
# ### PLS decomposition ----
# ### |||||||||||||||||||||||
# 
# Y_h_star <- TT_star %*% t(D_star) + Y.MEAN
# Y_dec_star <- TT_star %*% t(D_star) + Y.MEAN + YY[[3]] 
# 
# # Exact
# par(mfrow = c(1,1))
# plot(Y, main = "Y", col = rainbow(N), asp = 1,
#      xlab = Y.labs[1], ylab = Y.labs[2])
# points(Y_dec_star , pch = 8, col = rainbow(N))
# grid()
# 
# # Approximated
# plot_prediction(Y_h_star, Y, N, Y.labs)
# 
# sum((Y-Y_h_star)^2)/N
# 
# 
# ### PLS (Y-space) decomposition ----
# ### ||||||||||||||||||||||||||||||||
# 
# Y_hat_1_star = UU_star %*% t(V_star) + Y.MEAN
# 
# plot_prediction(Y_hat_1_star, Y, N, Y.labs)
# 
# sum((Y-Y_hat_1_star)^2)
# 
# ### PLS regression ----
# ### |||||||||||||||||||
# 
# Beta_star <- W_star %*% solve(t(C_star) %*% W_star, B_star %*% t(V_star))
# Y_hat_2_star <- Xc %*% Beta_star + Y.MEAN
# 
# plot_prediction(Y_hat_2_star, Y, N, Y.labs)
# 
# sum((Y-Y_hat_2_star)^2)/N
# 
# 
# # Comparison PLSR, PLSR* ----
# # |||||||||||||||||||||||||||
# 
# 
# t(W) %*% W_star
# 
# t(V) %*% V_star
# 
# TT
# TT_star
# 
# UU
# UU_star
# 
# C
# C_star
# 
# D
# D_star
# 
# Beta
# Beta_star
# 
# Y_hat_2
# Y_hat_2_star
