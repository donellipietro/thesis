# PLSR, NIPALS algorithm

PLSR_Donelli <- function(X, Y, A, deflation_Y = TRUE) {
  
  # Parameters
  toll <- 1e-6
  max_iter <- 100
  
  # Centering
  Xc = scale(X, scale = FALSE)
  X.mean = attr(Xc, "scaled:center")
  Yc = scale(Y, scale = FALSE)
  Y.mean = attr(Yc, "scaled:center")
  
  # Dimensions
  N <- dim(Xc)[1]
  S <- dim(Xc)[2]
  L <- dim(Yc)[2]
  
  # Mean values
  Y.MEAN = matrix(Y.mean, ncol = L, nrow = N, byrow = TRUE)
  X.MEAN = matrix(X.mean, ncol = S, nrow = N, byrow = TRUE)
  
  # Room for solutions
  TT <- matrix(0, nrow = N, ncol = A)
  
  # Extra
  C <- matrix(0, nrow = S, ncol = A)
  D <- matrix(0, nrow = L, ncol = A)
  
  # Initilization
  EE <- Xc
  FF <- Yc
  
  # Intermediate steps
  XX <- list()
  YY <- list()
  XX[[1]] <- EE
  YY[[1]] <- FF
  
  for(i in 1:A){
    
    # Normalization constants
    alpha <- norm(EE)
    beta <- norm(FF)
    
    # Initialization
    C.YX <- t(FF) %*% EE
    SVD <- svd(C.YX, nu = 1, nv = 1)
    c <- SVD$v # x-direction
    d <- SVD$u # y-direction
    
    c_change <- 1
    d_change <- 1
    delta_x <- 1
    delta_y <- 1
    n_iter <- 0
    
    while((c_change > toll || d_change > toll) && n_iter < max_iter){
      
      n_iter <- n_iter + 1
      
      # norm_coeff <- alpha * sum(c^2) + beta * sum(d^2)
      TT[,i] <- (alpha * EE %*% c + beta * FF %*% d)
      TT[,i] <- TT[,i]/norm(TT[,i],  "2")
      
      c <- t(EE) %*% TT[,i] # / sum(TT[,i]^2)
      d <- t(FF) %*% TT[,i] # / sum(TT[,i]^2)
      
      delta_x <- norm(c, "2")
      c <- c/delta_x
      delta_y <- norm(d, "2")
      d <- d/delta_y
    
      c_change <- norm(C[,i] - c, "2")
      d_change <- norm(D[,i] - d, "2")
        
      C[,i] <- c
      D[,i] <- d
      
    }
    
    C[,i] <- C[,i]*delta_x
    D[,i] <- D[,i]*delta_y
    
    print(delta_x)
    print(delta_y)
    
    print(n_iter)
    print(c_change)
    print(d_change)
    
    # Deflation
    EE <- EE - TT[,i] %*% t(C[,i])
    if(deflation_Y == TRUE)
      FF <- FF - TT[,i] %*% t(D[,i])

    
    # Intermediate steps
    XX[[i+1]] = EE
    YY[[i+1]] = FF
    
  }
  
  Beta <- C %*% solve(t(C) %*% C, t(D))
  
  Y_hat <- Xc %*% Beta + Y.MEAN
  X_hat <- TT %*% t(C) + X.MEAN 
  
  return(list(TT = TT,
              C = C,
              D = D,
              EE = EE,
              FF = FF,
              XX = XX,
              YY = YY,
              X.mean = X.mean,
              Y.mean = Y.mean,
              Beta = Beta,
              Y_hat = Y_hat,
              X_hat = X_hat))
  
}

save(PLSR_Donelli, file = "PLSR_Donelli.RData")