# PLSR, NIPALS algorithm

SIMPLS <- function(X, Y, A) {
  
  # X as a matrix
  X = as.matrix(X)
  
  # Dimensions
  N <- dim(X)[1]
  S <- dim(X)[2]
  L <- dim(Y)[2]
  
  # Centering
  Yc = scale(Y, scale = FALSE)
  Y.mean = attr(Yc, "scaled:center")
  Xc = scale(X, scale = FALSE)
  X.mean = attr(Xc, "scaled:center")
  
  # Mean values
  Y.MEAN = matrix(Y.mean, ncol = L, nrow = N, byrow = TRUE)
  X.MEAN = matrix(X.mean, ncol = S, nrow = N, byrow = TRUE)
  
  # Room for solutions
  R <- matrix(0, nrow = S, ncol = A)  # X block factor weights
  Q <- matrix(0, nrow = L, ncol = A)  # Y block factor weights
  TT <- matrix(0, nrow = N, ncol = 0) # X block factor scores
  UU <- matrix(0, nrow = N, ncol = A) # Y block factor scores
  P <- matrix(0, nrow = S, ncol = A)  # X block factor loadings
  V <- matrix(0, nrow = S, ncol = 0)  # Orthogonal loadings
  
  SS <- list()
  SS[[1]] <- t(X) %*% Yc
  
  for(i in 1:A){
    
    S <- SS[[i]]
    
    # Decomposition
    EIGEN <- eigen(t(S)%*%S)
    Q[,i] <- EIGEN$vectors[,1] # Y-block factor weights
    R[,i] <- S %*% Q[,i] # X-block factor weights
    
    TT <- cbind(TT, X %*% R[,i]) # X-block factor scores
    TT[,i] <- TT[,i] - mean(TT[,i])
    normt <- sqrt(sum(TT[,i]^2))
    TT[,i] <- TT[,i]/normt
    
    R[,i] <- R[,i]/normt
    
    P[,i] <- t(X) %*% TT[,i]
    Q[,i] <- t(Yc) %*% TT[,i]
    
    UU[,i] <- Yc %*% Q[,i]
    
    v <- P[,i]
    if(i > 1){
      v <- v - V %*% (t(V)%*%P[,i]) 
      UU[,i] <- UU[,i] - TT %*% (t(TT)%*%UU[,i]) 
    }
    V <- cbind(V, v/sqrt(sum(v^2)))
    
    SS[[i+1]] <- S - V[,i] %*% (t(V[,i]) %*% S)
    
  }
  
  Beta <- R %*% t(Q)
  
  Y_hat = Xc %*% Beta + Y.MEAN
  X_hat = TT %*% t(P) + X.MEAN
  
  h <- diag(TT %*% t(TT)) + 1/N
  varX <- diag(t(P) %*% P)/(N-1)
  varY <- diag(t(Q) %*% Q)/(N-1)
  
  return(list(W = R,
              V = Q,
              TT = TT,
              UU = UU,
              C = P,
              D = Q,
              V = V,
              SS = SS,
              Beta = Beta,
              Y_hat = Y_hat,
              X_hat = X_hat, 
              h = h,
              varX = varX,
              varY = varY))
  
}

save(SIMPLS, file = "SIMPLS.RData")