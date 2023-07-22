# Libraries
snapshot <- "/lab-share/RC-Data-Science-e2/Public/DST_software/pipeline_R_snapshot/r4.1-20220803-bch-dst-v1/"
.libPaths(snapshot)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(Nebulosa)
library(patchwork)
library(pbapply)

tenColors <- c('#a9d6e5', '#89c2d9', '#61a5c2', '#468faf', '#2c7da0', '#2a6f97', '#014f86', '#01497c', '#013a63', '#012a4a')
#c('#ffff3f', '#eeef20', '#dddf00', '#d4d700', '#bfd200','#aacc00', '#80b918', '#55a630', '#2b9348', '#007f5f')
source('/lab-share/RC-Data-Science-e2/Public/Daniel/rFunctions.R')

# Setting WD
setwd('/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler')

# # Loading the data
# samplesList <- list.files('data/')
# samplesData <- lapply(samplesList, function(sample){
#     X <- Read10X(paste0('data/', sample))
#     png(paste0('results/QC_', sample, '.png'), width = 2000, height = 1000, res = 300)
#     P <- plotMatrixQC(X) + plot_annotation(title = sample)
#     print(P)
#     dev.off()
#     X <- CreateSeuratObject(X, project = sample)
#     return(X)
# })

# # Merging data
# for(sample in seq_along(samplesData)[-1]){
#     samplesData[[1]] <- merge(samplesData[[1]], samplesData[[sample]])
#     samplesData[[sample]] <- new('Seurat')
# }
# P147 <- samplesData[[1]]
# Idents(P147) <- factor(P147$orig.ident, c('D0', 'D3', 'D6', 'D8', 'D10', 
#                           'D13', 'D18', 'D24', 'D30', 'D36'))


# # Preprocessing
# # P147D10 <- doPreprocessingHuman(P147, ID = 'P147', outFolder = 'results/D10', nDims = 10, resValue = 0.5)
# P147D30 <- doPreprocessingHuman(P147, ID = 'P147', outFolder = 'results/D30', nDims = 30, resValue = 0.5)
# P147D50 <- doPreprocessingHuman(P147, ID = 'P147', outFolder = 'results/D50', nDims = 50, resValue = 0.5)
# P147D100 <- doPreprocessingHuman(P147, ID = 'P147', outFolder = 'results/D100', nDims = 100, resValue = 0.5)

doPostProcessing <- function(X, ID, outFolder = 'results/'){
        # Replacing coloring of timepoints
        png(paste0(outFolder, '/UMAP2_', ID, '.png'), , width = 2000, height = 2000, res = 300)
        attr(X, 'U2') <- attr(X, 'U2') + 
            scale_color_manual(values = tenColors) + 
            theme_light() + 
            xlab('UMAP 1') + 
            ylab('UMAP 2') + 
            theme(panel.grid = element_blank(), 
                panel.border = element_blank())
        print(attr(X, 'U2'))
        dev.off()

        # Plotting provided markers
        providedMarkers <- sort(c('POU5F1', 'PAX6', 'NEFM', 'OC90', 'PAX2',  
                            'SOX10', 'PAX3', 'TWIST1', 'PRRX1',
                            'GSC', 'TGFBI', 'WNT6', 'HAND1', 'TFAP2A',
                            'GPM6A', 'TWIST2', 
                            'NKX2-5', 'FGF8', 'SIX1', 'PAX8', 'VGLL2',
                            'PRPH', 'MPZ', 'COL2A1', 
                            'OTOF'))
        
        X <- ScaleData(X, features = providedMarkers)
        png(paste0(outFolder, '/DP_', ID, '.png'), width = 2000, height = 2000, res = 300)
        attr(X, 'DP') <- DotPlot(X, features = providedMarkers, dot.min = 0.01, col.min = 0) +
                theme_light() + 
                theme(panel.grid = element_blank(), 
                        panel.border = element_blank()) +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 3))
        print(attr(X, 'DP'))
        dev.off()

        # Nebulosa Plots
        png(paste0(outFolder, '/N_', ID, '.png'), width = 4000, height = 3200, res = 300)
        attr(X, 'N') <- plot_density(X, providedMarkers, reduction = 'umap') & 
                    theme_light() & 
                    theme(panel.grid = element_blank(), 
                                    panel.border = element_blank(),
                                    legend.key.width= unit(0.5, 'mm'),
                                    plot.title = element_text(face = 3)) &
                    scale_color_gradient2(mid = '#f4f5f6', high = '#012a4a', midpoint = 0)
        print(attr(X, 'N'))
        dev.off()

        return(X)
}


# Loading objects
load('results/D10/P147.RData')
P147D10 <- doPostProcessing(X, ID = 'P147', outFolder = 'results/D10')

load('results/D30/P147.RData')
P147D30 <- doPostProcessing(X, ID = 'P147', outFolder = 'results/D30')

load('results/D50/P147.RData')
P147D50 <- doPostProcessing(X, ID = 'P147', outFolder = 'results/D50')
P147D50 <- RunUMAP(X, reduction = 'harmony', dims = seq_len(50), n.components= 3)

# Generating H5 files
SaveH5Seurat(P147D50, filename = 'P147D50', overwrite = TRUE)
Convert("P147D50.h5seurat", dest = "h5ad", overwrite = TRUE)

# Generating rBCS file
library(rBCS)
ExportSeurat(P147D50, "P147D50.bcs",overwrite=TRUE)


load('results/D100/P147.RData')
P147D100 <- doPostProcessing(X, ID = 'P147', outFolder = 'results/D100')
P147D100 <- RunUMAP(X, reduction = 'harmony', dims = seq_len(100), n.components= 3)

# Generating H5 files
SaveH5Seurat(P147D100, filename = 'P147D100', overwrite = TRUE)
Convert("P147D100.h5seurat", dest = "h5ad", overwrite = TRUE)

# Generating rBCS file
library(rBCS)
ExportSeurat(P147D100, "P147D100.bcs",overwrite=TRUE)

# Adding Annotations
## Loading Object
load('results/D100/P147.RData')
P147D100 <- X

## Loading Annotations
KKAnnotations <- readRDS('KKannot2.rds')

## Transfering Annotations
P147D100$CT <- KKAnnotations@meta.data[,5]
rm(KKAnnotations)

## Short ID
shortID <- strsplit(as.vector(P147D100$CT), '\\(')
shortID <- unlist(lapply(shortID, function(X){X[2]}))
shortID <- gsub('\\)', '', shortID)
P147D100$ID <- shortID

# Setting Short ID as identity
Idents(P147D100) <- factor(P147D100$ID, sort(unique(P147D100$ID)))

# Writting UMAP
png('t.png', width = 2000, height = 2000, res = 300)
UMAPPlot(P147D100, label = TRUE) + 
theme_light() +
theme(legend.position = 'None', 
        panel.grid = element_blank(),
        panel.border = element_blank()) +
xlab('UMAP 1') +
ylab('UMAP 2')
dev.off()

CM <- doMarkersPlot(P147D100)

png('t1.png', width =5000, height = 2000, res = 300)
CM + xlab('Cell Type')
dev.off()

png('t2.png', width =2000, height = 2500, res = 300)
doNumericUMAP(P147D100)
dev.off()