library(plot3D)


load("../data/data_fdaPDE.RData")
load("../data/load_score_pietro.RData")
load("functions/plot_field.RData")

nodes <- mesh$nodes

plot_field(nodes, X_new, load1, range(load1), "Load 1", TRUE)
plot_field(nodes, X_new, load2, range(load2), "Load 2", TRUE)
plot_field(nodes, X_new, load3, range(load3), "Load 3", TRUE)
