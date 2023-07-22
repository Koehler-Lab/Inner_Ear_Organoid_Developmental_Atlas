snapshot <- "/lab-share/RC-Data-Science-e2/Public/DST_software/pipeline_R_snapshot/r4.1-20220803-bch-dst-v1/"
.libPaths(snapshot)
source('/lab-share/RC-Data-Science-e2/Public/Daniel/rFunctions.R')
library(Seurat)
library(SeuratDisk)

# Set WD
setwd('/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler')

# Loading data
load('P147D100_Annotated.RData')
P147D100@meta.data <- as.data.frame.array(P147D100@meta.data)
P147D100 <- ScaleData(P147D100, rownames(P147D100))

SaveH5Seurat(P147D100, filename = 'P147D100.h5Seurat')
Convert('P147D100.h5Seurat', dest = "h5ad")
