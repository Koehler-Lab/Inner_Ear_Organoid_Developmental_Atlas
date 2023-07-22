require(Seurat)
snapshot <- "/lab-share/RC-Data-Science-e2/Public/DST_software/pipeline_R_snapshot/r4.1-20220803-bch-dst-v1/"
.libPaths(snapshot)
source('/lab-share/RC-Data-Science-e2/Public/Daniel/rFunctions.R')
library(rBCS)

# Setting WD
setwd('/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler')

# Loading Data
load('/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler/P147D100_Annotated.RData')

# Saving Object
ExportSeurat(P147D100, 'P147D100_UMAP2D.bcs', overwrite = TRUE)
