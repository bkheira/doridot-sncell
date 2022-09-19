# This script handles the demultiplexing of one snRNAseq run

#Date: 06/09/22
### By Alix Silvert
### Based on https://satijalab.org/seurat/articles/hashing_vignette.html Compiled in 2022-01-11

##################
### INITIATING ###
##################

library(Matrix)
library(ggplot2)
library(Seurat)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

rds_file_path = args[1]
samples_info_path = args[2]
output = args[3]
outputDirectory = args[4]

samples_informations = read.delim(samples_info_path)

### Load filtered data ###
print(rds_file_path)
seurat_object = readRDS(rds_file_path)

seurat_object <- NormalizeData(seurat_object, assay = "HTO", normalization.method = "CLR")
seurat_object <- HTODemux(seurat_object, assay = "HTO", positive.quantile = 0.99)

### Change Idents to the correct names ###
new_names = samples_informations$sample_code
names(new_names) = samples_informations$sample_name
seurat_object <- RenameIdents(seurat_object, new_names)



### Output DEMUX statistics ###
table_statistics <-
  seurat_object@meta.data %>%
  select(hash.ID) %>%
  table

write.table(table_statistics, file.path(outputDirectory, "number_cells_per_HTO.txt"))

### Output QC per sample and negative and doublets ###
plotQCs <- function(sampleName, srtObject, outputDir){
  data_to_print = srtObject@meta.data
  data_to_print <- data_to_print %>%
      filter(hash.ID == sampleName)
  p <- ggplot(data_to_print, aes(y = nCount_RNA, x = "nCount_RNA"))
  p <- p + theme_bw()
  p <- p + xlab(sampleName)
  p <- p + geom_violin(fill = "#F07DEA", color = NA)
  p <- p + geom_jitter( width = 0.5, alpha = 0.1, size = 0.5)
  p <- p + theme(
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

  pdf(file.path(outputDir, paste("nCount_postFilters_", sampleName, ".pdf", sep="")), width = 3, height = 5)
  print(p)
  dev.off()

  p <- ggplot(data_to_print, aes(y = nFeature_RNA, x = "nFeature_RNA"))
  p <- p + theme_bw()
  p <- p + xlab(sampleName)
  p <- p + geom_violin(fill = "#AF0171", color = NA)
  p <- p + geom_jitter( width = 0.5, alpha = 0.1, size = 0.5)
  p <- p + theme(
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

  pdf(file.path(outputDir, paste("nFeature_postFilters_", sampleName, ".pdf", sep="")), width = 3, height = 5)
  print(p)
  dev.off()


  p <- ggplot(data_to_print, aes(y = percent.mt, x = "percent.mt"))
  p <- p + theme_bw()
  p <- p + xlab(sampleName)
  p <- p + geom_violin(fill = "#7FB77E", color = NA)
  p <- p + geom_jitter( width = 0.5, alpha = 0.1, size = 0.5)
  p <- p + theme(
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

  pdf(file.path(outputDir, paste("percent_mithochondrial_postFilters_", sampleName, ".pdf", sep="")), width = 3, height = 5)
  print(p)
  dev.off()

  p <- ggplot(data_to_print, aes(x = nCount_RNA, y = nFeature_RNA))
  p <- p + theme_bw()
  p <- p + geom_point(color = "#A460ED")


  pdf(file.path(outputDir, paste("nCount_nFeature_comparison_postFilters_", sampleName, ".pdf", sep="")), width = 5, height = 5)
  print(p)
  dev.off()
}

lapply(samples_informations$sample_name, plotQCs, srtObject = seurat_object, outputDir = outputDirectory)


### Calculate a tSNE embedding of the HTO data



seurat_object_copy <- subset(seurat_object, idents = "Negative", invert = TRUE)
DefaultAssay(seurat_object_copy) <- "HTO"
seurat_object_copy <- ScaleData(seurat_object_copy, features = rownames(seurat_object_copy),
    verbose = FALSE)
seurat_object_copy<- RunPCA(seurat_object_copy, features = rownames(seurat_object_copy), approx = FALSE)
seurat_object_copy<- RunTSNE(seurat_object_copy, dims = 1:nrow(samples_informations), perplexity = 100, check_duplicates = FALSE)
p <- DimPlot(seurat_object_copy)

pdf(file.path(outputDirectory, "HTO_tSNE_plot.pdf"), width = 7, height = 7)
print(p)
dev.off()

DefaultAssay(seurat_object_copy) <- "RNA"


###Clustering of HTO
p <- HTOHeatmap(seurat_object_copy, assay = "HTO", ncells = 500)

pdf(file.path(outputDirectory, "HTO_heatmap_plot.pdf"), width = 7, height = 7)
print(p)
dev.off()


### Remove cells from doublets and negative  and plot a tSNE (on a copy) ###
seurat_object <- subset(seurat_object, idents = "Negative", invert = TRUE)
seurat_object <- subset(seurat_object, idents = "Doublet", invert = TRUE)

seurat_object_tSNE <- FindVariableFeatures(seurat_object)
seurat_object_tSNE <- ScaleData(seurat_object_tSNE, features = VariableFeatures(seurat_object_tSNE))
seurat_object_tSNE <- RunPCA(seurat_object_tSNE)
seurat_object_tSNE <- RunTSNE(seurat_object_tSNE, dims = 1:5, perplexity = 100)
p <- DimPlot(seurat_object_tSNE)

pdf(file.path(outputDirectory, "RNA_tSNE_plot.pdf"), width = 7, height = 7)
print(p)
dev.off()

### Save results ### 
saveRDS(seurat_object, file = output)
