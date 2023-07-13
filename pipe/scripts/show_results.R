

load("data/locations/mesh_1.RData")
f <- read.csv("results/f.csv", header = FALSE)
#f <- f[mesh.unrolled$extra == 0]

FEMbasis <- create.FEM.basis(mesh.unrolled)


# plot(FEM(as.numeric(f), FEMbasis))
# points3d(sensors, col = "purple", pch = ".", cex = 4)

sensors.real <- sensors.unrolled[sensors.unrolled[,"extra"]==0, 1:2]

FEMFunction.new <- FEM(coeff = as.numeric(f), FEMbasis =  FEMbasis)
plot(FEMFunction.new)

sensor.estimate <- eval.FEM(FEMFunction.new, locations = sensors.real)


X <- read.csv("data/X_alpha10_40psi.csv", header = TRUE)[2:950]





par(mfrow = c(1,1))
plot(sensor.estimate, type = "l")
for(i in 1:1){
  points(as.numeric(X[i,]), type = "l") #, col = rainbow(10)[1])
  # for(j in 2:10){
  #   points(as.numeric(distance), as.numeric(strains[[experiments.names[i]]][j,]), type = "l", col = rainbow(10)[i])
  # }
  # abline(v = c(1, 495, 581, 581+949)*0.01, col = "red", lwd = 2, lty = 2)
}
points(sensor.estimate, type = "l", col = "blue")
