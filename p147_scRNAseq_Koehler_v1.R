# Libraries
snapshot <- "/lab-share/RC-Data-Science-e2/Public/DST_software/pipeline_R_snapshot/r4.1-20220803-bch-dst-v1/"
.libPaths(snapshot)
library(Seurat)
library(ggplot2)
library(Nebulosa)
library(patchwork)
library(pbapply)

tenColors <- c('#a9d6e5', '#89c2d9', '#61a5c2', '#468faf', '#2c7da0', '#2a6f97', '#014f86', '#01497c', '#013a63', '#012a4a')
#c('#ffff3f', '#eeef20', '#dddf00', '#d4d700', '#bfd200','#aacc00', '#80b918', '#55a630', '#2b9348', '#007f5f')
source('/lab-share/RC-Data-Science-e2/Public/Daniel/rFunctions.R')

# Setting WD
setwd('/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler')

# Loading the data
samplesList <- list.files('data/')
samplesData <- lapply(samplesList, function(sample){
    X <- Read10X(paste0('data/', sample))
    png(paste0('results/QC_', sample, '.png'), width = 2000, height = 1000, res = 300)
    P <- plotMatrixQC(X) + plot_annotation(title = sample)
    print(P)
    dev.off()
    X <- CreateSeuratObject(X, project = sample)
    return(X)
})

# Merging data
for(sample in seq_along(samplesData)[-1]){
    samplesData[[1]] <- merge(samplesData[[1]], samplesData[[sample]])
    samplesData[[sample]] <- new('Seurat')
}
P147 <- samplesData[[1]]
Idents(P147) <- factor(P147$orig.ident, c('D0', 'D3', 'D6', 'D8', 'D10', 
                          'D13', 'D18', 'D24', 'D30', 'D36'))

# Removing objects
rm(samplesList, samplesData)

# QC Integrated Dataset
png('results/QC1_P147.png', width = 2000, height = 1000, res = 300)
plotMatrixQC(P147@assays$RNA@counts) + plot_annotation(title = 'P147 - Before QC')
dev.off()

# Identifying and removinf doublets
P147$isDoublet <- isDoublet(P147)
# 2887 Doublets will be removed
P147 <- P147[,!P147$isDoublet]

# Quality controls
P147 <- scQC(P147)

# QC Integrated Dataset
png('results/QC2_P147.png', width = 2000, height = 1000, res = 300)
plotMatrixQC(P147@assays$RNA@counts) + plot_annotation(title = 'P147 - After QC')
dev.off()

# Data Processing
P147 <- NormalizeData(P147)
P147 <- FindVariableFeatures(P147, selection.method = 'vst', nfeatures = 2000)
P147 <- ScaleData(P147)
set.seed(1)
P147 <- RunPCA(P147, npcs = 100)

# Defining the number of PCS
png('results/elbowPlot_P147.png', width = 2000, height = 2000, res = 300)
ElbowPlot(P147, 100) +
    theme_light() + 
    theme(panel.grid = element_blank(), 
        panel.border = element_blank())
dev.off()

# UMAP Before Integation
set.seed(1)
P147 <- RunUMAP(P147, reduction = 'pca', dims = 1:10)

png('results/nointegratedUMAP.png', width = 2000, height = 2000, res = 300)
UMAPPlot(P147) + 
    scale_color_manual(values = tenColors) + 
    theme_light() + 
    xlab('UMAP 1') + 
    ylab('UMAP 2') + 
    theme(panel.grid = element_blank(), 
        panel.border = element_blank())
dev.off()

# Data Integration
set.seed(1)
P147 <- harmony::RunHarmony(P147, 'orig.ident', max.iter.harmony = 1e3, dims = 1:100)

# UMAP After Integration
set.seed(1)
P147 <- RunUMAP(P147, reduction = 'harmony', dims = 1:10)

png('results/integratedUMAP.png', width = 2000, height = 2000, res = 300)
UMAPPlot(P147) + 
    scale_color_manual(values = tenColors) + 
    theme_light() + 
    xlab('UMAP 1') + 
    ylab('UMAP 2') + 
    theme(panel.grid = element_blank(), 
        panel.border = element_blank())
dev.off()

# Spectral clustering
set.seed(1)
P147 <- FindNeighbors(P147, reduction = 'harmony', dims = 1:10)
set.seed(1)
P147 <- FindClusters(P147, resolution = 0.5)

# Clustered UMAP
png('results/clusteredUMAP.png', width = 2000, height = 2000, res = 300)
UMAPPlot(P147, label = TRUE, repel = TRUE) + 
    theme_light() + 
    xlab('UMAP 1') + 
    ylab('UMAP 2') + 
    theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        legend.position = 'None')
dev.off()

# Cluster Markers
clusterMarkers <- doMarkersPlot(P147)

# Plotting the markers
png('results/clusterMarkers.png', width = 4000, height = 2000, res = 300)
print(clusterMarkers) +
        theme(panel.grid = element_blank(), 
        panel.border = element_blank())
dev.off()

# Plotting provided markers
providedMarkers <- sort(c('POU5F1', 'PAX6', 'NEFM', 'OC90', 'PAX2',  
                    'SOX10', 'PAX3', 'TWIST1', 'PRRX1',
                    'GSC', 'TGFBI', 'WNT6', 'HAND1', 'TFAP2A',
                    'GPM6A', 'TWIST2', 
                    'NKX2-5', 'FGF8', 'SIX1', 'PAX8', 'VGLL2',
                    'PRPH', 'MPZ', 'COL2A1', 
                    'OTOF'))
P147 <- ScaleData(P147, features = providedMarkers)
png('results/providedMarkers.png', width = 2000, height = 2000, res = 300)
DotPlot(P147, features = providedMarkers, dot.min = 0.01, col.min = 0) +
        theme_light() + 
        theme(panel.grid = element_blank(), 
                panel.border = element_blank()) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 3))

dev.off()

# Nebulosa Plots
png('results/nebulosaMarkers.png', width = 4000, height = 3200, res = 300)
# , dot.min = 0.01, col.min = 0
plot_density(P147, providedMarkers, reduction = 'umap') & 
            theme_light() & 
            theme(panel.grid = element_blank(), 
                            panel.border = element_blank(),
                            legend.key.width= unit(0.5, 'mm'),
                            plot.title = element_text(face = 3)) &
            scale_color_gradient2(mid = '#f4f5f6', high = '#012a4a', midpoint = 0)
dev.off()
# png('results/sox2UMAP.png', width = 2000, height = 2000, res = 300)
# FeaturePlot(P147, features = 'SOX2', order = TRUE) + 
#     theme_light() + 
#     xlab('UMAP 1') + 
#     ylab('UMAP 2') + 
#     theme(panel.grid = element_blank(), 
#     panel.border = element_blank())
# dev.off()

# png('results/sox2UMAP.png', width = 1200, height = 1000, res = 300)
# plot_density(P147, 'SOX2', reduction = 'umap') & 
#             theme_light() & 
#             theme(panel.grid = element_blank(), 
#                             panel.border = element_blank(),
#                             legend.key.width= unit(0.5, 'mm'),
#                             plot.title = element_text(face = 3)) &
#             scale_color_gradient2(mid = '#f4f5f6', high = '#012a4a', midpoint = 0)
# dev.off()
