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
cell_cycle_markers_path = args[4]

seurat_data = readRDS(file = rds_file_path)
cell_cycle_markers = readRDS(cell_cycle_markers_path)

######################################
### GETTING THE DATA TO REGRESS ON ###
######################################


GPhase <- cell_cycle_markers %>%
  filter(phase == "G2/M")
SPhase <- cell_cycle_markers %>%
  filter(phase == "S")

###################
### NORMALISING ###
###################

# We are gonne use a normalization based on SCTransform while correcting for
# phase cycle (as this is not of interest to us)


split_seurat = SplitObject(seurat_data, split.by = "HTO_classification")


#We split the seurat object before normalisation because of possible
#pre sequencing factors
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=GPhase, s.features=SPhase)
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("percent.mt", "G2M.Score", "S.Score"))
}


saveRDS(split_seurat, file = output)
