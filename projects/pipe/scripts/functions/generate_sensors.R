generate_sensors <- function(R, alpha, windings, step) {
  w = 2*pi
  pitch <- 2*pi*R*tan(alpha)
  conv = step/sqrt(R^2*w^2+pitch^2)
  
  t = seq(0, windings, length = 1000)
  K = 0:(windings/conv)
  
  helic <- data.frame(x = R*cos(w*t),
                      y = R*sin(w*t),
                      z = pitch*t)
  sensors <- data.frame(x = R*cos(w*conv*K),
                        y = R*sin(w*conv*K),
                        z = pitch*conv*K)
  
  return(sensors)
}

save(generate_sensors, file = "scripts/functions/generate_sensors.RData")