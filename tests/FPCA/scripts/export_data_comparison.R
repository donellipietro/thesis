rm(list = ls())
graphics.off()
def.par = par()

dir <- paste("data/fPCA_vs_fPCACS/comparison_test/", sep = '')
if (!file.exists(dir)){
  dir.create(dir)
}

################################################################################

N <- 400
load("data/fPCA_vs_fPCACS/data_nnodes10.RData")

X <- data_centered[1:N, ]
locations <- locations

dir <- paste("data/fPCA_vs_fPCACS/comparison_test/", "mesh10/", sep = '')
if (!file.exists(dir)){
  dir.create(dir)
}

write.csv(mesh_data$boundary, paste(dir, "boundary.csv", sep = ''))
write.csv(mesh_data$edges, paste(dir, "edges.csv", sep = ''))
write.csv(mesh_data$elements, paste(dir, "elements.csv", sep = ''))
write.csv(mesh_data$neigh, paste(dir, "neigh.csv", sep = ''))
write.csv(mesh_data$nodes, paste(dir, "points.csv", sep = ''))

dir <- paste("data/fPCA_vs_fPCACS/comparison_test/", "2D_test_comparison/", sep = '')
if (!file.exists(dir)){
  dir.create(dir)
}

write.csv(X, paste(dir, "X10.csv", sep = ''))
write.csv(locations, paste(dir, "locations10.csv", sep = ''))


################################################################################

N <- 50
load("data/fPCA_vs_fPCACS/data_nnodes60.RData")

X <- data_centered[1:N, ]
locations <- locations

dir <- paste("data/fPCA_vs_fPCACS/comparison_test/", "mesh60/", sep = '')
if (!file.exists(dir)){
  dir.create(dir)
}

write.csv(mesh_data$boundary, paste(dir, "boundary.csv", sep = ''))
write.csv(mesh_data$edges, paste(dir, "edges.csv", sep = ''))
write.csv(mesh_data$elements, paste(dir, "elements.csv", sep = ''))
write.csv(mesh_data$neigh, paste(dir, "neigh.csv", sep = ''))
write.csv(mesh_data$nodes, paste(dir, "points.csv", sep = ''))

dir <- paste("data/fPCA_vs_fPCACS/comparison_test/", "2D_test_comparison/", sep = '')
if (!file.exists(dir)){
  dir.create(dir)
}

write.csv(X, paste(dir, "X60.csv", sep = ''))
write.csv(locations, paste(dir, "locations60.csv", sep = ''))

