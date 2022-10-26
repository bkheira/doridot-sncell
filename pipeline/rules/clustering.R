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

library(future)

plan("multicore", workers = 4)

args = commandArgs(trailingOnly = TRUE)

rds_file_path = args[1]
output = args[2]
outputDirectory = args[3]

###Loading data
seurat_file = readRDS(rds_file_path)

###Normalisation of data###
seurat_file = SCTransform(seurat_file, vars.to.regress = c("percent.mt", "G2M.Score", "S.Score"))


###Clustering###
#We will follow the Seurat default
seurat_file <- RunPCA(seurat_file)
seurat_file <- FindNeighbors(seurat_file)
seurat_file <- RunUMAP(seurat_file, dims = 1:10, future.seed=TRUE)

seurat_file <- FindClusters(seurat_file, future.seed=TRUE)



###Print UMAP of Clusters###
pdf(file.path(outputDirectory, "UMAP_clusterisation.pdf"),
    width = 5, height = 5)
print(DimPlot(seurat_file))
dev.off()

###Get Clusters identifiers###
markers <- FindAllMarkers(seurat_file, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers, file = file.path(outputDirectory, "clusters_markers.tsv"), sep="\t")

###Get Heatmap of markers###
top10Markers <- markers %>%
  group_by(cluster) %>%
  top_n(n=10, wt = avg_log2FC)

pdf(file.path(outputDirectory, "Clusters_Heatmap.pdf"),
      width = 5, height = 5)
DoHeatmap(seurat_file, features = top10Markers$gene) + NoLegend()
dev.off()


###Saving seurat file###
saveRDS(seurat_file, file=output)
