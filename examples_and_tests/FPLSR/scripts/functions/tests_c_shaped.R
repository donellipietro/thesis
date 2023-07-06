generate_X <- function(points, X.index){
  
  if(X.index == 1){
    a1 = stats::rnorm(1, mean = 1, sd = 0.2)
    a2 = stats::rnorm(1, mean = 1, sd = 0.2)
    z = a1*cos(2*pi*points[,1]) + a2*cos(2*pi*points[,2]) + 1
  }
  if(X.index == 2){
    a1 = stats::rnorm(1, mean = 1, sd = 0.1)
    a2 = stats::rnorm(1, mean = 1, sd = 0.1)
    a3 = stats::rnorm(1, mean = 1, sd = 0.1)
    f <- function(x, y, z = 1){
      coe <- function(x,y){
        1/2*sin(5*pi*x)*exp(-x^2)+1
      }
      a1*sin(2*pi*(coe(y,1)*x*cos(z-2)-a2*y*sin(z-2)))*cos(2*pi*(coe(y,1)*x*cos(z-2+pi/2)+a3*coe(x,1)*y*sin((z-2)*pi/2)))
    }
    # Exact solution (pointwise at nodes)
    z = f(points[,1], points[,2])
  }
  
  return(z)
  
}

generate_B <- function(nodes, B.index) {
  # Generate B(x, y):
  if (B.index == 1) {
    r <- 0.2
    z <- 5*exp(-((nodes[, 1] + 0.5)^2 + (nodes[, 2] + 0)^2)/( 2*r^2 )) +
         5*exp(-((nodes[, 1] - 3.0)^2 + (nodes[, 2] - 0.5)^2)/( 2*r^2 )) + 
         5*exp(-((nodes[, 1] - 3.0)^2 + (nodes[, 2] + 0.5)^2)/( 2*r^2 ))
  }else if (B.index == 2) {
    # top right corner
    r <- 0.2
    r.x <- 0.5
    r.y <- 0.2
    z <- 5*exp(-((nodes[, 1] + 0.5)^2 + (nodes[, 2] + 0)^2)/( 2*r^2 )) +
         5*exp(-((nodes[, 1] - 1.5)^2)/( 2*r.x^2 ) - ((nodes[, 2] - 0.5)^2)/( 2*r.y^2 )) + 
         5*exp(-((nodes[, 1] - 1.5)^2)/( 2*r.x^2 ) - ((nodes[, 2] + 0.5)^2)/( 2*r.y^2 ))
  }else if (B.index == 3) {
    r <- 0.2
    r.x <- 0.5
    r.y <- 0.2
    z <- 5*exp(-((nodes[, 1] - 1.5)^2)/( 2*r.x^2 ) - ((nodes[, 2] - 0.5)^2)/( 2*r.y^2 )) + 
         5*exp(-((nodes[, 1] - 1.5)^2)/( 2*r^2 ) - ((nodes[, 2] + 0.5)^2)/( 2*r^2 ))
  }else if (B.index == 4) {
    r.x <- 0.5
    r.y <- 0.2
    z <- 5*exp(-((nodes[, 1] - 1.5)^2)/( 2*r.x^2 ) - ((nodes[, 2] - 0.1)^2)/( 2*r.y^2 )) * (nodes[, 2]>0)
  }else if (B.index == 5) {
    r.x <- 2.0
    r.y <- 0.5
    z <- 5*exp(-((nodes[, 1] - 3.5)^2)/( 2*r.x^2 ) - ((nodes[, 2] - 0.5)^2)/( 2*r.y^2 )) * (nodes[,2] > 0 | nodes[,1] < 0) -
         5*exp(-((nodes[, 1] - 3.5)^2)/( 2*r.x^2 ) - ((nodes[, 2] + 0.5)^2)/( 2*r.y^2 )) * (nodes[,2] < 0 | nodes[,1] < 0)
  }
  B <- as.matrix(z)
  return(B)
}

save(generate_X, generate_B, file = "scripts/functions/tests_c_shaped.RData")


