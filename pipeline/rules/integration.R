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
outputDir = args[3]
rds_files_paths = args[4:length(args)] #Note : this is a vector of lenght "Number of files to analyse"


all_seurat_objects = unlist(sapply(rds_files_paths,readRDS), recursive = FALSE)

integ_features <- SelectIntegrationFeatures(object.list = all_seurat_objects)
integ_anchors <- FindIntegrationAnchors(object.list = all_seurat_objects)
integrated_data <- IntegrateData(anchorset = integ_anchors)

#Print UMAP colored with sample names
integrated_data<- ScaleData(integrated_data, verbose = FALSE)
integrated_data<- RunPCA(integrated_data, verbose = FALSE)
integrated_data<- RunUMAP(integrated_data, reduction = "pca", dims = 1:30)

pdf(file.path(outputDir, "UMAP_postIntegration.pdf"),
    width = 5, height = 5)
print(DimPlot(integrated_data))
dev.off()



saveRDS(integrated_data, file = output)
