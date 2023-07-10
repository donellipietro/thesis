# install dependencies (execute this only once, if you not have Rcpp installed yet)
install.packages("Rcpp")
install.packages("RcppEigen")

# load Rcpp library
library(Rcpp)

setwd("~pietrodonelli/shared-folder/thesis/fdaPDE/wrappers/R/")
# update RCppExports.cpp
compileAttributes(".")
# install fdaPDE
install.packages(".", type="source", repos=NULL)


