# This script handles the normalisation
# of data

# Date: 03/10/22 By Alix Silvert

################## INITIATING ###
library(Seurat)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)

rds_file_path = args[1]
samples_info_path = args[2]
output = args[3]
outputDirectory = args[4]

seurat_data = readRDS(file = rds_file_path)

################### NORMALISING ###
