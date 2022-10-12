# This script handles the
# clusterisation of an integrated samples

# Date: 06/09/22 By Alix Silvert Based
# on
# https://satijalab.org/seurat/articles/hashing_vignette.html
# Compiled in 2022-01-11

###INITIATING ###

library(Matrix)
library(ggplot2)
library(Seurat)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)

rds_file_path = args[1]
output = args[2]
outputDirectory = args[3]

###Loading data
seurat_file = loadRDS(rds_file_path)


###Saving seurat file###
savreRDS(seurat_file, file=outputDirectory)
