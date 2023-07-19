rm(list = ls())
graphics.off()


# ||||||||||||||
# Libraries ----
# ||||||||||||||

# library(softImpute)


# ||||||||||||||||
# Import data ----
# ||||||||||||||||

UNROLLED = TRUE

# Parameters
step = 0.01
start <- 1000
end <- 2629

# Data
DATA <- read.csv("../data/Test2_strain.csv")
pressure <- DATA[, dim(DATA)[2]]
strains_all <- DATA[, start:end]
distance <- step*(1:dim(strains_all)[2])

# Esperiments
experiments.pressure <- unique(pressure)
experiments.number <- length(experiments.pressure)
experiments.names <- paste(experiments.pressure, "psi", sep = '')
experiments.indexes <- list()
for(i in 1:experiments.number){
  experiments.indexes[[experiments.names[i]]] <- which(pressure == experiments.pressure[i])
}

experiment.type.names <- c("alpha10","alpha20")
experiment.type.indexes <- list("alpha10" = c(565, 565+949-1), "alpha20" = c(10, 10+495-1))

if(UNROLLED == TRUE){
  # Sensors location
  sensors <- list()
  for(i in 1:length(experiment.type.names))
    sensors[[experiment.type.names[i]]] <- read.csv(paste("../data/fdaPDE_data/locations/locations_unrolled_", i ,".csv", sep = ''))
}

# Strains
strains <- list()
for(i in 1:experiments.number){
  strains[[experiments.names[i]]] <- strains_all[experiments.indexes[[experiments.names[i]]],]
}

par(mfrow = c(2,1))
for(i in 1:experiments.number){
  plot(as.numeric(distance), as.numeric(strains[[experiments.names[i]]][1,]),
       type = "l") #, col = rainbow(10)[1])
  # for(j in 2:10){
  #   points(as.numeric(distance), as.numeric(strains[[experiments.names[i]]][j,]), type = "l", col = rainbow(10)[i])
  # }
  abline(v = c(experiment.type.indexes[[1]], experiment.type.indexes[[2]])*0.01, col = "red", lwd = 2, lty = 2)
}


# |||||||||||||||||||
# NaN imputation ----
# |||||||||||||||||||

for(i in 1:experiments.number){
  
  strains[[experiments.names[i]]][is.na(strains[[experiments.names[i]]])] = 0
  
  print(sum(is.na(strains[[experiments.names[i]]])))
  
  # fit <- softImpute(strains[[experiments.names[i]]],
  #                   rank=min(dim(strains[[experiments.names[i]]]))-1,
  #                   lambda=30)
  # strains[[experiments.names[i]]] <- complete(as.numeric(strains[[experiments.names[i]]]), fit)
  
}

# ||||||||||||||||
# Export data ----
# ||||||||||||||||

dir <- "../data/fdaPDE_data/data"
if (!file.exists(dir)){
  dir.create(dir)
}


for(type in experiment.type.names){
  for(pressure in experiments.names){
    range <- experiment.type.indexes[[type]]
    temp <- strains[[pressure]][,range[1]:range[2]]
    unrolled <- ""
    if(UNROLLED == TRUE){
      temp <- temp[,sensors[[type]]$id]
      unrolled <- "unrolled"
    }
    write.csv(temp, file = paste(paste("../data/fdaPDE_data/data/X",unrolled,type,pressure, sep = '_'), ".csv", sep = ""))
  }
}




