# This script handles the demultiplexing of one snRNAseq run

#Date: 06/09/22
### By Alix Silvert
### Based on https://satijalab.org/seurat/articles/hashing_vignette.html Compiled in 2022-01-11


library(Matrix)
library(data.table)
library(ggplot2)
library(Seurat)

args = commandArgs(trailingOnly=TRUE)

print(args)
print(args[4])


# samples_info_path = args[1]
# output = args[2]
# filterthreshold= c(args[3], args[4])


samples_info_path = "/Users/asilvert/Programmation/Doridot/doridot-sncell/pipeline/outputs/pouloup.tsv"
output = "blaaaah"
filterthreshold= c(20, 200)


samples_informations = fread(samples_info_path, header = TRUE)
print(samples_informations)

### Load filtered data ###

barcode.path <- paste(samples_informations$path_to_mtx_filtered_directory[1], "barcodes.tsv.gz", sep="")
feature.path <- paste(samples_informations$path_to_mtx_filtered_directory[1], "features.tsv.gz", sep="")
matrix.path <- paste(samples_informations$path_to_mtx_filtered_directory[1], "matrix.mtx.gz", sep="")

print(barcode.path)
print(matrix.path)

#Loading data
expression_matrix <- ReadMtx(
  mtx = matrix.path, features = feature.path,
  cells = barcode.path
)

expression_matrix_HTO = expression_matrix[samples_informations$sample_name,]
expression_matrix_RNA = expression_matrix[!(rownames(expression_matrix) %in% samples_informations$sample_name),]

seurat_object = CreateSeuratObject(counts = expression_matrix_RNA)

seurat_object[["HTO"]] <- CreateAssayObject(counts = expression_matrix_HTO)



#We'll try No and With filters to check
NF_seurat = seurat_object
NF_seurat <- NormalizeData(NF_seurat, assay = "HTO", normalization.method = "CLR")
NF_seurat <- HTODemux(NF_seurat, assay = "HTO", positive.quantile = 0.99)

table(NF_seurat$HTO_classification.global)
FeatureScatter(NF_seurat, feature1 = "hto_2-2", feature2 = "hto_5-2")
HTOHeatmap(NF_seurat, assay = "HTO", ncells = 5000)


write.table(NF_seurat@meta.data[, c("HTO_classification", "HTO_classification.global")], file = "no_filters_classification.tsv", sep="\t", quote = F)


#With Filters
YF_seurat = seurat_object
YF_seurat[["percent.mt"]] <- PercentageFeatureSet(YF_seurat, pattern = "^mt-")

YF_seurat <- subset(YF_seurat, subset = nFeature_RNA > 800 & nFeature_RNA < 7500 & percent.mt < 5)
VlnPlot(YF_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
YF_seurat <- NormalizeData(YF_seurat, assay = "HTO", normalization.method = "CLR")
YF_seurat <- HTODemux(YF_seurat, assay = "HTO", positive.quantile = 0.99)

table(YF_seurat$HTO_classification.global)
FeatureScatter(YF_seurat, feature1 = "hto_2-2", feature2 = "hto_5-2")
HTOHeatmap(YF_seurat, assay = "HTO", ncells = 5000)

write.table(YF_seurat@meta.data[, c("HTO_classification", "HTO_classification.global")], file = "filters_classification.tsv", sep="\t", quote = F)


#GetAssayData(object = YF_seurat, slot = "data", assay = "HTO")

# mat <- readMM(file = matrix.path)
# feature.names = read.delim(feature.path,
#                            header = FALSE,
#                            stringsAsFactors = FALSE)
# barcode.names = read.delim(barcode.path,
#                            header = FALSE,
#                            stringsAsFactors = FALSE)
# colnames(mat) = barcode.names$V1
# rownames(mat) = feature.names$V1
