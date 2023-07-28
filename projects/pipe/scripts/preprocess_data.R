# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Crack Detection Project: Preprocess data %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()
def.par = par()

cat("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
cat("\n%% Crack Detection Project: Preprocess data %%")
cat("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")


# ||||||||||||||
# Libraries ----
# ||||||||||||||

# |||||||||||||||
# Parameters ----
# |||||||||||||||

# Mesh
step = 0.01

# Data
test.name = "Test1"

# Data trimming
start <- 1000
end <- 2629

# Plots
PLOT = TRUE


# ||||||||||||||||
# Import data ----
# ||||||||||||||||

cat("\n")
cat("\n# |||||||||||||||||||||||")
cat("\n# Importing raw data ----")
cat("\n# |||||||||||||||||||||||\n\n")

# Data
DATA <- read.csv(paste("data/", test.name, "_strain.csv", sep = ""))
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

# Strains
strains <- list()
for(i in 1:experiments.number){
  strains[[experiments.names[i]]] <- strains_all[experiments.indexes[[experiments.names[i]]],]
}


# |||||||||||||||||||
# NaN imputation ----
# |||||||||||||||||||

cat("\n")
cat("\n# |||||||||||||||||||")
cat("\n# NaN imputation ----")
cat("\n# |||||||||||||||||||\n")


for(i in 1:experiments.number){
  
  strains[[experiments.names[i]]][is.na(strains[[experiments.names[i]]])] = 0
  
  # print(sum(is.na(strains[[experiments.names[i]]])))
  
  # fit <- softImpute(strains[[experiments.names[i]]],
  #                   rank=min(dim(strains[[experiments.names[i]]]))-1,
  #                   lambda=30)
  # strains[[experiments.names[i]]] <- complete(as.numeric(strains[[experiments.names[i]]]), fit)
  
}


# |||||||||||||||||||
# Visualize data ----
# |||||||||||||||||||

cat("\n")
cat("\n# ||||||||||||||||||||")
cat("\n# Visualize data ----")
cat("\n# ||||||||||||||||||||\n\n")


if(PLOT){
  par(mfrow = c(2,1))
  for(i in 1:experiments.number){
    ylim = range(strains[[experiments.names[i]]])
    N <- dim(strains[[experiments.names[i]]])[1]
    plot(NA, xlim = range(distance), ylim = ylim, main = experiments.names[i])
    grid()
    for(j in 1:N){
      points(as.numeric(distance), as.numeric(strains[[experiments.names[i]]][j,]), type = "l", col = rainbow(N)[j])
    }
    abline(v = c(experiment.type.indexes[[1]], experiment.type.indexes[[2]])*0.01, col = "red", lwd = 2, lty = 2)
  }
}


# ||||||||||||||||
# Export data ----
# ||||||||||||||||

cat("\n")
cat("\n# |||||||||||||||||||")
cat("\n# Exporting data ----")
cat("\n# |||||||||||||||||||\n")

dir <- "data/strains/"
if (!file.exists(dir)){
  dir.create(dir)
}

dir <- paste(dir, test.name, "/", sep = "")
if (!file.exists(dir)){
  dir.create(dir)
}

for(type in experiment.type.names){
  for(pressure in experiments.names){
    range <- experiment.type.indexes[[type]]
    temp <- strains[[pressure]][,range[1]:range[2]]
    write.csv(temp, file = paste(dir ,paste("X", type, pressure, sep = '_'), ".csv", sep = ""))
  }
}




