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
    # centered at (0.5, 0.5):
    r  <-  0.2 # set the r parameter
    z  <-  5*exp(-((nodes[, 1] - 0.5 )^2 + (nodes[, 2] - 0.5 )^2)/( 2*r^2 ))
  }else if (B.index == 2) {
    # top right corner
    z <- 5*exp(-((nodes[, 1] - 0.75)^2 + (nodes[, 2] - 0.75)^2)/( 2*0.2^2 ))
  }else if (B.index == 3) {
    # bottom left + top right corner:
    z <- 5*exp(-((nodes[, 1] - 0.75)^2 + (nodes[, 2] - 0.75)^2)/( 2*0.25^2 )) +
      5*exp(-((nodes[, 1] - 0.1)^2 + (nodes[, 2] - 0.1)^2)/( 2*0.25^2 ))
  }else if (B.index == 4) {
    # semi-circumference:
    r  <-  0.2 # set the r parameter
    z  <-  5*exp(-((nodes[, 1] - 0.5 )^2 + (nodes[, 2] )^2)/( 2*r^2 ))
  }else if (B.index == 5) {
    # monkey saddle
    z = ((nodes[, 1]*4 - 2)^3 - 3*(nodes[, 1]*4 - 2)*((nodes[, 2]*4 - 2)^2))
  }
  B <- as.matrix(z)
  return(B)
}

save(generate_X, generate_B, file = "scripts/functions/tests_unit_square.RData")


