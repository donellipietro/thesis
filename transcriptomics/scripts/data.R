rm(list=ls())
graphics.off()


DATA <- read.csv("../data/hippo_all.csv", header = TRUE)

summary(DATA)
# colnames(DATA)
colnames.cell_id <- colnames(DATA)[1]
colnames.gene <- colnames(DATA)[2:160]
colnames.coord <- colnames(DATA)[161:162]
colnames.animal_id <- colnames(DATA)[170]
colnames.subject <- colnames(DATA)[166]


subjects <- unique(DATA[, colnames.subject])





mask <- read.csv("../data/masks/BFN-Corn-37-bottom-brain-mask.csv", header = FALSE)
      

plot(DATA[DATA[, colnames.subject] == subjects[7], colnames.coord],
     asp = 1, pch = '.')
points(mask, asp = 1, col = "red") 

