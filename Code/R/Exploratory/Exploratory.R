library(R.matlab)
library(tidyverse)

project_path <- "/Users/mizuhirosuzuki/Dropbox/MFdiffusion_replication/"
file.path(project_path, "datav4.0/Matlab Replication/India Networks/adjacencymatrix.mat")

adjacency_mat <- readMat(file.path(project_path, "datav4.0/Matlab Replication/India Networks/adjacencymatrix.mat"))

leader_mat <- read_tsv(file.path(project_path, "datav4.0/Matlab Replication/India Networks/HHhasALeader1.csv"), col_names = FALSE)

takeup_mat <- read_tsv(file.path(project_path, "datav4.0/Matlab Replication/India Networks/MF1.csv"), col_names = FALSE)

giant_mat <- read_tsv(file.path(project_path, "datav4.0/Matlab Replication/India Networks/inGiant1.csv"), col_names = FALSE)

covariate_mat <- read_tsv(file.path(project_path, "datav4.0/Matlab Replication/India Networks/hhcovariates1.csv"), col_names = FALSE)

A <- adjacency_mat$X[[1]][[1]]
Sec <- (A %*% A > 0)

for (i in 1:nrow(A)){
  Sec[i,i] <- 0
}
