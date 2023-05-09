# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# %% Multivariate Principal Least Squares Regression %% #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

rm(list = ls())
graphics.off()

# ||||||||||||||
# Functions ----
# ||||||||||||||

load("PLSR.RData")
load("PLSR_star.RData")


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

path1 <- "../../fdaPDE/test/data/models/FPLSR/2D_test0/test1/"
path2 <- "../../fdaPDE/test/data/models/FPLSR/2D_test0/test2/"
path3 <- "../../fdaPDE/test/data/models/FPLSR/2D_test0/test3/"

X = read.csv(paste(path1, "X.csv", sep = ''), header = TRUE)[,2:3601]

Y1 = matrix(read.csv(paste(path1, "Y.csv", sep = ''), header = TRUE)[,2], ncol = 1)
Y3 = matrix(read.csv(paste(path3, "Y.csv", sep = ''), header = TRUE)[,2], ncol = 1)
Y = cbind(Y1, Y3)

## De-trending ----
## |||||||||||||||

Xc = scale(X, scale = FALSE)
X.mean = attr(Xc, "scaled:center")

Yc = scale(Y, scale = FALSE)
Y.mean = attr(Yc, "scaled:center")

## Dimensions ----
## |||||||||||||||

N <- dim(X)[1]
S <- dim(X)[2]
L <- dim(Y)[2]


## Plot ----
## |||||||||


# Y-spaces
par(mfrow = c(2,1))
plot(1:N, Y[,1], main = "Y1")
abline(h = Y.mean[1], lty = 2, col = "black")
grid()
plot(1:N, Y[,2], main = "Y2")
abline(h = Y.mean[2], lty = 2, col = "black")
grid()

# X-space
par(mfrow = c(1,1))
image(matrix(as.numeric(X[1,]), 60, 60), main = "X")


# |||||||||
# PLSR ----
# |||||||||

Y.MEAN = matrix(Y.mean, ncol = L, nrow = N, byrow = TRUE)

plsr = PLSR(Xc, Yc, 2)

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

Beta <- W %*% solve(t(C) %*% W, B %*% t(V))
Y_hat_2 <- Xc %*% Beta + Y.MEAN

image(matrix(as.numeric(Beta[,1]), 60, 60), main = "X")
plot_prediction(Y_hat_2, Y, N, c("","",""))


# ||||||||||
# PLSR* ----
# ||||||||||

plsr_star = PLSR_star(Xc, Yc, 2)

W_star = plsr_star[["W"]]
V_star = plsr_star[["V"]]
TT_star = plsr_star[["TT"]]
UU_star = plsr_star[["UU"]]
B_star = plsr_star[["B"]]
C_star = plsr_star[["C"]]
D_star = plsr_star[["D"]]
EE_star = plsr_star[["EE"]]
FF_star = plsr_star[["FF"]]
XX_star = plsr_star[["XX"]]
YY_star = plsr_star[["YY"]]

Beta_star <- W_star %*% solve(t(C_star) %*% W_star, B_star %*% t(V_star))
Y_hat_2_star <- Xc %*% Beta_star + Y.MEAN

image(matrix(as.numeric(Beta_star[,1]), 60, 60), main = "X")
plot_prediction(Y_hat_2_star, Y, N, c("","",""))


