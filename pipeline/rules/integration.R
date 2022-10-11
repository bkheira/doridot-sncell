#this script aims at integrating several files that have been demultiplexed
#This script is intended to be used inside a snakemake preparing_pipeline

#by ALix Silvert 10/10/2022


####INITIATING####
## Load libraries
library(Seurat)

#Get files
args = commandArgs(trailingOnly = TRUE)

samples_info_path = args[1]
output = args[2]
rds_files_paths = args[3:length(args)] #Note : this is a vector of lenght "Number of files to analyse"

for (i in rds_files_paths){
  print(i)
}

all_seurat_objects = unlist(sapply(rds_files_paths,readRDS), recursive = FALSE)

integ_features <- SelectIntegrationFeatures(object.list = all_seurat_objects)
integ_anchors <- FindIntegrationAnchors(object.list = all_seurat_objects)
integrated_data <- IntegrateData(anchorset = integ_anchors)

saveRDS(integrated_data, file = output)
