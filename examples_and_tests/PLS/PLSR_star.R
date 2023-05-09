# PLSR, no deflation

PLSR_star <- function(Xc, Yc, A) {
  
  # Dimensions
  N <- dim(Xc)[1]
  S <- dim(Xc)[2]
  L <- dim(Yc)[2]
  
  # Room for solutions
  W <- matrix(0, nrow = S, ncol = A)
  V <- matrix(0, nrow = L, ncol = A)
  TT <- matrix(0, nrow = N, ncol = A)
  UU <- matrix(0, nrow = N, ncol = A)
  B <- matrix(0, nrow = A, ncol = A)
  
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
  
  # Decomposition
  C.YX <- t(FF) %*% EE
  SVD <- svd(C.YX, nu = A, nv = A)
  
  for(i in 1:A){
    
    W[,i] <- SVD$v[,i] # x-direction
    if(L != 1){
      V[,i] <- SVD$u[,i] # y-direction
    }
    else{
      V[,i] = 1
    }
    
    TT[,i] <- EE %*% W[,i] # x-component
    UU[,i] <- FF %*% V[,i] # y-component
    
    tt = sum(TT[,i]^2)
    
    # Regression
    C[,i] <- t(EE) %*% TT[,i] / tt
    D[,i] <- t(FF) %*% TT[,i] / tt
    B[i,i] <- t(UU[,i]) %*% TT[,i] / tt
    
    # Deflation
    # EE <- EE - TT[,i] %*% t(C[,i])
    # FF <- FF - B[i,i] * TT[,i] %*% t(V[,i])
    
    # Intermediate steps
    XX[[i+1]] = EE
    YY[[i+1]] = FF
    
  }
  
  return(list(W = W,
              V = V,
              TT = TT,
              UU = UU,
              B = B,
              C = C,
              D = D,
              EE = EE,
              FF = FF,
              XX = XX,
              YY = YY))
  
}

save(PLSR_star, file = "PLSR_star.RData")