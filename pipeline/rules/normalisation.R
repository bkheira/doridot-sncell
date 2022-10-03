# This script handles the normalisation of data

#Date: 03/10/22
### By Alix Silvert

##################
### INITIATING ###
##################
library(Seurat)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

rds_file_path = args[1]
samples_info_path = args[2]
output = args[3]
outputDirectory = args[4]
cell_cycle_markers_path = args[5]

seurat_data = readRDS(file = rds_file_path)
cell_cycle_markers = readRDS(cell_cycle_markers_path)

###################
### NORMALISING ###
###################

# We are gonne use a normalization based on SCTransform while correcting for
# phase cycle (as this is not of interest to us)

split_seurat = SplitObject(seurat_data, split.by = "sample")
#This crate a list of sub object, each will be normalised on their own.
print(cell_cycle_markers)
