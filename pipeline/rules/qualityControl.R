# 12 September 2022 QC script By Alix
# Silvert This script is meant to be
# used inside a specific snakemake
# pipeline, if you are modifying it,
# please make sure it is a copy or your
# modifications are appropriate for the
# pipeline

################## INITIATING ###

## Load libraries
library(Seurat)
library(ggplot2)

## Treating inputs
args = commandArgs(trailingOnly = TRUE)

samples_info_path = args[1]
output = args[2]
MIN_RNA_PER_DROPLET = as.numeric(args[3])
MAX_RNA_PER_DROPLET = as.numeric(args[4])
MIN_FEATURES_PER_DROPLET = as.numeric(args[5])
MAX_FEATURES_PER_DROPLET = as.numeric(args[6])
MAX_PERCENT_MT_READS = as.numeric(args[7])
SPECIE = args[8]
OUTPUT_DIRECTORY = args[9]



## Loading data
samples_informations = read.delim(samples_info_path,
    header = TRUE)

sample_group_name = samples_informations$multiplexing_group_name[1]

barcode.path <- paste(samples_informations$path_to_mtx_filtered_directory[1],
    "barcodes.tsv.gz", sep = "")
feature.path <- paste(samples_informations$path_to_mtx_filtered_directory[1],
    "features.tsv.gz", sep = "")
matrix.path <- paste(samples_informations$path_to_mtx_filtered_directory[1],
    "matrix.mtx.gz", sep = "")

expression_matrix <- ReadMtx(mtx = matrix.path,
    features = feature.path, cells = barcode.path)

expression_matrix_HTO = expression_matrix[samples_informations$sample_name,
    ]
expression_matrix_RNA = expression_matrix[!(rownames(expression_matrix) %in%
    samples_informations$sample_name), ]

seurat_object = CreateSeuratObject(counts = expression_matrix_RNA)
seurat_object[["HTO"]] <- CreateAssayObject(counts = expression_matrix_HTO)

################## RUN THE QC ###
################## Computing
################## mithochondrial DNA
################## Percent
if (SPECIE == "Human") {
    seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object,
        pattern = "^MT-", assay = "RNA")
} else if (SPECIE == "Mouse") {
    seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object,
        pattern = "^mt-", assay = "RNA")
} else {
    stop("'SPECIE' Should have a value of either 'Human' or 'Mouse'")
}

## Printing pre-filters plots with
## visualisation of filters
data_to_print = seurat_object@meta.data

p <- ggplot(data_to_print, aes(y = nCount_RNA,
    x = "nCount_RNA"))
p <- p + theme_bw()
p <- p + xlab(sample_group_name)
p <- p + geom_violin(fill = "#F07DEA", color = NA)
p <- p + geom_jitter(width = 0.5, alpha = 0.1,
    size = 0.5)
p <- p + geom_hline(yintercept = MIN_RNA_PER_DROPLET,
    linetype = "dashed", color = "red")
p <- p + geom_hline(yintercept = MAX_RNA_PER_DROPLET,
    linetype = "dashed", color = "red")
p <- p + theme(axis.text.x = element_blank(),
    axis.ticks.x = element_blank())

pdf(file.path(OUTPUT_DIRECTORY, "nCount_preFilters.pdf"),
    width = 3, height = 5)
print(p)
dev.off()

p <- ggplot(data_to_print, aes(y = nFeature_RNA,
    x = "nFeature_RNA"))
p <- p + theme_bw()
p <- p + xlab(sample_group_name)
p <- p + geom_violin(fill = "#AF0171", color = NA)
p <- p + geom_jitter(width = 0.5, alpha = 0.1,
    size = 0.5)
p <- p + geom_hline(yintercept = MIN_FEATURES_PER_DROPLET,
    linetype = "dashed", color = "red")
p <- p + geom_hline(yintercept = MAX_FEATURES_PER_DROPLET,
    linetype = "dashed", color = "red")
p <- p + theme(axis.text.x = element_blank(),
    axis.ticks.x = element_blank())

pdf(file.path(OUTPUT_DIRECTORY, "nFeature_preFilters.pdf"),
    width = 3, height = 5)
print(p)
dev.off()


p <- ggplot(data_to_print, aes(y = percent.mt,
    x = "percent.mt"))
p <- p + theme_bw()
p <- p + xlab(sample_group_name)
p <- p + geom_violin(fill = "#7FB77E", color = NA)
p <- p + geom_jitter(width = 0.5, alpha = 0.1,
    size = 0.5)
p <- p + geom_hline(yintercept = MAX_PERCENT_MT_READS,
    linetype = "dashed", color = "red")
p <- p + theme(axis.text.x = element_blank(),
    axis.ticks.x = element_blank())

pdf(file.path(OUTPUT_DIRECTORY, "percent_mitochondrialDNA_preFilters.pdf"),
    width = 3, height = 5)
print(p)
dev.off()

p <- ggplot(data_to_print, aes(x = nCount_RNA,
    y = nFeature_RNA))
p <- p + theme_bw()
p <- p + geom_point(color = "#A460ED")
p <- p + geom_hline(yintercept = MIN_FEATURES_PER_DROPLET,
    linetype = "dashed", color = "red")
p <- p + geom_hline(yintercept = MAX_FEATURES_PER_DROPLET,
    linetype = "dashed", color = "red")
p <- p + geom_vline(xintercept = MIN_RNA_PER_DROPLET,
    linetype = "dashed", color = "red")
p <- p + geom_vline(xintercept = MAX_RNA_PER_DROPLET,
    linetype = "dashed", color = "red")

pdf(file.path(OUTPUT_DIRECTORY, "nCout_nFeatures_comparison_preFilters.pdf"),
    width = 5, height = 5)
print(p)
dev.off()

## Actually running the filters
seurat_object <- subset(seurat_object, subset = nCount_RNA >
    MIN_RNA_PER_DROPLET & nCount_RNA < MAX_RNA_PER_DROPLET)
seurat_object <- subset(seurat_object, subset = nFeature_RNA >
    MIN_FEATURES_PER_DROPLET & nFeature_RNA <
    MAX_FEATURES_PER_DROPLET)
seurat_object <- subset(seurat_object, subset = percent.mt <
    MAX_PERCENT_MT_READS)

# Printing the rest
data_to_print = seurat_object@meta.data

p <- ggplot(data_to_print, aes(y = nCount_RNA,
    x = "nCount_RNA"))
p <- p + theme_bw()
p <- p + xlab(sample_group_name)
p <- p + geom_violin(fill = "#F07DEA", color = NA)
p <- p + geom_jitter(width = 0.5, alpha = 0.1,
    size = 0.5)
p <- p + theme(axis.text.x = element_blank(),
    axis.ticks.x = element_blank())

pdf(file.path(OUTPUT_DIRECTORY, "nCount_postFilters.pdf"),
    width = 3, height = 5)
print(p)
dev.off()

p <- ggplot(data_to_print, aes(y = nFeature_RNA,
    x = "nFeature_RNA"))
p <- p + theme_bw()
p <- p + xlab(sample_group_name)
p <- p + geom_violin(fill = "#AF0171", color = NA)
p <- p + geom_jitter(width = 0.5, alpha = 0.1,
    size = 0.5)
p <- p + theme(axis.text.x = element_blank(),
    axis.ticks.x = element_blank())

pdf(file.path(OUTPUT_DIRECTORY, "nFeature_postFilters.pdf"),
    width = 3, height = 5)
print(p)
dev.off()


p <- ggplot(data_to_print, aes(y = percent.mt,
    x = "percent.mt"))
p <- p + theme_bw()
p <- p + xlab(sample_group_name)
p <- p + geom_violin(fill = "#7FB77E", color = NA)
p <- p + geom_jitter(width = 0.5, alpha = 0.1,
    size = 0.5)
p <- p + theme(axis.text.x = element_blank(),
    axis.ticks.x = element_blank())

pdf(file.path(OUTPUT_DIRECTORY, "percent_mitochondrialDNA_postFilters.pdf"),
    width = 3, height = 5)
print(p)
dev.off()

p <- ggplot(data_to_print, aes(x = nCount_RNA,
    y = nFeature_RNA))
p <- p + theme_bw()
p <- p + geom_point(color = "#A460ED")
p <- p + geom_hline(yintercept = MIN_FEATURES_PER_DROPLET,
    linetype = "dashed", color = "red")
p <- p + geom_hline(yintercept = MAX_FEATURES_PER_DROPLET,
    linetype = "dashed", color = "red")
p <- p + geom_vline(xintercept = MIN_RNA_PER_DROPLET,
    linetype = "dashed", color = "red")
p <- p + geom_vline(xintercept = MAX_RNA_PER_DROPLET,
    linetype = "dashed", color = "red")

pdf(file.path(OUTPUT_DIRECTORY, "nCout_nFeatures_comparison_postFilters.pdf"),
    width = 5, height = 5)
print(p)
dev.off()

## Save the Rds file
saveRDS(seurat_object, file = output)
