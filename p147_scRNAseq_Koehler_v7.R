# Libraries
library(patchwork)
library(Seurat)
snapshot <- "/lab-share/RC-Data-Science-e2/Public/DST_software/pipeline_R_snapshot/r4.1-20220803-bch-dst-v1/"
.libPaths(snapshot)
library(SeuratDisk)
library(ggplot2)
library(Nebulosa)
library(pbapply)
library(ggrepel)
library(harmony)
library(dplyr)


tpColors <- tenColors <- c('#a9d6e5', '#89c2d9', '#61a5c2', '#468faf', '#2c7da0', '#2a6f97', '#014f86', '#01497c', '#013a63', '#012a4a')
ctColors <- c('#AF0508', '#F1070B', '#F94144', 
'#9A3C09', '#D4520C', '#F3722C', 
'#9B5805', '#D67907', '#F8961E', 
'#B48007', '#F6AF0A', '#F9C74F', 
'#527735', '#70A349', '#90BE6D', 
'#26604F', '#35846C', '#43AA8B', 
'#314352', '#445C71', '#577590')
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

# doPostProcessing <- function(X, ID, outFolder = 'results/'){
#         # Replacing coloring of timepoints
#         png(paste0(outFolder, '/UMAP2_', ID, '.png'), , width = 2000, height = 2000, res = 300)
#         attr(X, 'U2') <- attr(X, 'U2') + 
#             scale_color_manual(values = tenColors) + 
#             theme_light() + 
#             xlab('UMAP 1') + 
#             ylab('UMAP 2') + 
#             theme(panel.grid = element_blank(), 
#                 panel.border = element_blank())
#         print(attr(X, 'U2'))
#         dev.off()

#         # Plotting provided markers
#         providedMarkers <- sort(c('POU5F1', 'PAX6', 'NEFM', 'OC90', 'PAX2',  
#                             'SOX10', 'PAX3', 'TWIST1', 'PRRX1',
#                             'GSC', 'TGFBI', 'WNT6', 'HAND1', 'TFAP2A',
#                             'GPM6A', 'TWIST2', 
#                             'NKX2-5', 'FGF8', 'SIX1', 'PAX8', 'VGLL2',
#                             'PRPH', 'MPZ', 'COL2A1', 
#                             'OTOF'))
        
#         X <- ScaleData(X, features = providedMarkers)
#         png(paste0(outFolder, '/DP_', ID, '.png'), width = 2000, height = 2000, res = 300)
#         attr(X, 'DP') <- DotPlot(X, features = providedMarkers, dot.min = 0.01, col.min = 0) +
#                 theme_light() + 
#                 theme(panel.grid = element_blank(), 
#                         panel.border = element_blank()) +
#                 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 3))
#         print(attr(X, 'DP'))
#         dev.off()

#         # Nebulosa Plots
#         png(paste0(outFolder, '/N_', ID, '.png'), width = 4000, height = 3200, res = 300)
#         attr(X, 'N') <- plot_density(X, providedMarkers, reduction = 'umap') & 
#                     theme_light() & 
#                     theme(panel.grid = element_blank(), 
#                                     panel.border = element_blank(),
#                                     legend.key.width= unit(0.5, 'mm'),
#                                     plot.title = element_text(face = 3)) &
#                     scale_color_gradient2(mid = '#f4f5f6', high = '#012a4a', midpoint = 0)
#         print(attr(X, 'N'))
#         dev.off()

#         return(X)
# }

# Adding Annotations
## Loading Object
load('results/D100/P147.RData')
P147D100 <- X

## Loading Annotations
KKAnnotations <- readRDS('KKannot2.rds')

## Transfering Annotations
P147D100 <- P147D100[,colnames(KKAnnotations)]
P147D100$CT <- KKAnnotations@meta.data[,5]

## Adding Hair Cells
set.seed(1)
P147D100 <- FindClusters(P147D100, resolution = 1)
P147D100$CT[P147D100$RNA_snn_res.1 == 33] <- 'Hair Cells (HC)'

## Short ID
shortID <- strsplit(as.vector(P147D100$CT), '\\(')
shortID <- unlist(lapply(shortID, function(X){X[2]}))
shortID <- gsub('\\)', '', shortID)
P147D100$ID <- shortID

# Setting Short ID as identity
Idents(P147D100) <- factor(P147D100$ID, sort(unique(P147D100$ID)))

# Labels
LabelPosition <- lapply(split(as.data.frame(P147D100@reductions$umap@cell.embeddings), Idents(P147D100)), function(X){apply(X,2,median)})
Labels <- names(LabelPosition)
LabelPosition <- do.call(rbind.data.frame, LabelPosition)
colnames(LabelPosition) <- c('X', 'Y')
LabelPosition$Label <- Labels

Label1 <- unique(P147D100$CT)
Label2 <- unique(P147D100$ID)
matchColors2 <- ctColors
names(matchColors2) <- Label2

Idents(P147D100) <- factor(Idents(P147D100), Label2)
# save(P147D100, file = 'P147D100_Annotated.RData')
# load('P147D100_Annotated.RData')

for(tp in c('D13', 'D18', 'D30')){
        for(ct in c('OEP', 'OEPD')){
                DE <- FindMarkers(P147D100[,P147D100$orig.ident %in% tp], ct, logfc.threshold = 0.5, only.pos = TRUE)
                outFile <- paste0('OEP-OEPD/ctMarkers_', ct, '-', tp, '.csv')
                write.csv(DE, outFile)
        }
}

## Velocity without added info
# Plotting
load('scVelocity_Results.RData')
umap_arrows <- arrows$garrows %>%
    as.data.frame() %>%
    mutate(x2 = x0 + (x1 - x0) * 1,
           y2 = y0 + (y1 - y0) * 1)


png('t.png', width = 1500, height = 1500, res = 300)
UMAPPlot(P147D100) +
theme_void() +
theme(legend.position = 'None', 
        panel.grid = element_blank(),
        panel.border = element_blank()) +
        geom_curve(data = umap_arrows, curvature = 0.2, angle = 45,
                 aes(x = x0, xend = x2, y = y0, yend = y2),
                 size = .3,
                 arrow = arrow(length = unit(1, "points"), type = "closed"),
                 colour = "grey20", alpha = 0.8) +
                 xlab('UMAP 1') +
xlab('UMAP 1') +
ylab('UMAP 2') +
scale_color_manual(values = matchColors2[levels(Idents(P147D100))])
dev.off()

png('t0.png', width = 1500, height = 1500, res = 300)
UMAPPlot(P147D100) +
theme_void() +
theme(legend.position = 'None', 
        panel.grid = element_blank(),
        panel.border = element_blank()) +
        # geom_curve(data = umap_arrows, curvature = 0.2, angle = 45,
        #          aes(x = x0, xend = x2, y = y0, yend = y2),
        #          size = .3,
        #          arrow = arrow(length = unit(1, "points"), type = "closed"),
        #          colour = "grey20", alpha = 0.8) +
                 xlab('UMAP 1') +
xlab('UMAP 1') +
ylab('UMAP 2') +
scale_color_manual(values = matchColors2[levels(Idents(P147D100))])
dev.off()

##### Figure 5
D <- P147D100@reductions$umap@cell.embeddings
D <- as.data.frame.array(D)
D$ID <- P147D100$ID
D$Color <- matchColors2[D$ID]
D$TP <- P147D100$orig.ident
D$Color[!D$ID %in% c('NE', 'NEP', 'NCC', 'aPPE', 'pPPE', 'OEPD', 'OEP', 'NEU')] <- 'gray95'
D <- D[order(D$Color, decreasing = TRUE),]
P5 <- ggplot(D, aes(UMAP_1, UMAP_2)) +
geom_point(color = D$Color, size = 0.01) +
theme_light() +
theme(legend.position = 'None', 
        panel.grid = element_blank(),
        panel.border = element_blank()) +
xlab('UMAP 1') +
ylab('UMAP 2') +
scale_color_manual(values = matchColors2[levels(Idents(P147D100[,P147D100$orig.ident == 'D3']))]) +
geom_text_repel(data = LabelPosition[LabelPosition$Label %in% c('NE', 'NEP', 'NCC', 'aPPE', 'pPPE', 'OEPD', 'OEP', 'NEU'),], 
mapping = aes(X,Y, label = Label), 
bg.color = 'white', 
box.padding = 0,   
min.segment.length = 0) 

png('F5-UMAP.png', width = 1500, height = 1500, res = 300)
print(P5)
dev.off()

F5Data <- P147D100[,P147D100$ID %in% c('NE', 'NEP', 'NCC', 'aPPE', 'pPPE', 'OEPD', 'OEP', 'NEU')]

Idents(F5Data) <- factor(gsub('D', '',F5Data$orig.ident), c(0, 3, 6, 8, 10, 13, 18, 24, 30, 36)) 
png('F5-DP.png', width = 1500, height = 950, res = 300)
F5Data <- ScaleData(F5Data, c('SOX2','NEUROG1', 'NEUROD1', 'POU4F1', 'ISL1', 'STMN2', 'PRPH'))
DotPlot(F5Data, 
features= rev(c('SOX2','NEUROG1', 'NEUROD1', 'POU4F1', 'ISL1', 'STMN2', 'PRPH')), 
dot.min = 0.01, col.min = 0, scale.by = 'size') +
theme_light() + 
coord_flip() + 
theme(panel.grid = element_blank(), panel.border = element_blank()) + 
theme(axis.text.y = element_text(face = 3)) +
theme(legend.position = 'right') +
scale_color_gradient2(mid = 'gray90', high = '#F94144', midpoint = 0) +
ylab('Days') +
xlab(NULL)
dev.off()

Idents(F5Data) <- F5Data$ID
set.seed(1)
F5Data <- RunUMAP(F5Data, reduction = 'harmony', dims = 1:5)
save(F5Data, file = 'SS3.RData')

load('SS3.RData')
load('scVelocity_SS3.RData')
library(dplyr)

F5LabelPosition <- lapply(split(as.data.frame(F5Data@reductions$umap@cell.embeddings), Idents(F5Data)), function(X){apply(X,2,median)})
F5Labels <- names(F5LabelPosition)
F5LabelPosition <- do.call(rbind.data.frame, F5LabelPosition)
colnames(F5LabelPosition) <- c('X', 'Y')
F5LabelPosition$Label <- F5Labels

umap_arrows <- arrows$garrows %>%
    as.data.frame() %>%
    mutate(x2 = x0 + (x1 - x0) * 2,
           y2 = y0 + (y1 - y0) * 2)

png('F5-R.png', width = 1100, height = 1100, res = 300)
UMAPPlot(F5Data) +
theme_light() +
theme(legend.position = 'None', 
        panel.grid = element_blank(),
        panel.border = element_blank()) +
xlab('UMAP 1') +
ylab('UMAP 2') +
scale_color_manual(values = matchColors2[levels(Idents(F5Data))]) +
geom_curve(data = umap_arrows, curvature = 0.2, angle = 25,
                 aes(x = x0, xend = x2, y = y0, yend = y2),
                 size = .2,
                 arrow = arrow(length = unit(1, "points"), type = "closed"),
                 colour = "grey20", alpha = 0.8) +
geom_text_repel(data = F5LabelPosition, 
mapping = aes(X,Y, label = Label), 
bg.color = 'white', 
box.padding = 0,   
min.segment.length = 0) 
dev.off()

load('SS3.RData')
load('scVelo_SS3.RData')

png('umapSS3-2.png', width = 1500, height = 1500, res = 300)
UMAPPlot(F5Data) +
theme_light() +
theme(legend.position = 'bottom', 
        panel.grid = element_blank(),
        panel.border = element_blank()) +
geom_segment(data=VF, 
        mapping=aes(
            x=start.1, 
            y=start.2, 
            xend=end.1, 
            yend=end.2), 
        size = 0.2,
        lineend = 'round', 
        linejoin = 'round',
        arrow=arrow(length=unit(0.5, "mm"))) +
xlab('UMAP 1') + 
ylab('UMAP 2') + 
scale_color_manual(values = matchColors2[levels(Idents(F5Data))])
dev.off()


png('F5-T.png', width = 1500, height = 1600, res = 300)
Idents(F5Data) <- factor(gsub('D', '',F5Data$orig.ident), c(0, 3, 6, 8, 10, 13, 18, 24, 30, 36)) 
levels(Idents(F5Data)) <- paste0('Day ', c(0, 3, 6, 8, 10, 13, 18, 24, 30, 36))
UMAPPlot(F5Data) +
theme_light() +
theme(legend.position = 'bottom', 
        panel.grid = element_blank(),
        panel.border = element_blank()) +
xlab('UMAP 1') +
ylab('UMAP 2') +
scale_color_manual(values = tenColors) +
geom_text_repel(data = F5LabelPosition, 
mapping = aes(X,Y, label = Label), 
bg.color = 'white', 
box.padding = 0,   
min.segment.length = 0) 
dev.off()

png('IODA-MKI67.png', width = 1200, height = 1000, res = 300)
plot_density(P147D100, 'MKI67', reduction = 'umap') + 
theme_void() +
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.width= unit(1, 'mm')) +
scale_color_gradient2(mid = 'gray90', high = '#012a4a', midpoint = 0)
dev.off()

png('IODA-POU5F1-1.png', width = 1200, height = 1000, res = 300)
plot_density(P147D100, 'POU5F1', reduction = 'umap') + 
theme_void() +
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.width= unit(1, 'mm')) +
scale_color_gradient2(mid = 'gray90', high = '#012a4a', midpoint = 0)
dev.off()

png('IODA-FOXH1-1.png', width = 1200, height = 1000, res = 300)
plot_density(P147D100, 'FOXH1', reduction = 'umap') + 
theme_void() +
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.width= unit(1, 'mm')) +
scale_color_gradient2(mid = 'gray90', high = '#012a4a', midpoint = 0)
dev.off()

png('IODA-PAX8-1.png', width = 1200, height = 1000, res = 300)
plot_density(P147D100, 'PAX8', reduction = 'umap') + 
theme_void() +
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.width= unit(1, 'mm')) +
scale_color_gradient2(mid = 'gray90', high = '#012a4a', midpoint = 0)
dev.off()

png('IODA-ISL1-1.png', width = 1200, height = 1000, res = 300)
plot_density(P147D100, 'ISL1', reduction = 'umap') + 
theme_void() +
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.width= unit(1, 'mm')) +
scale_color_gradient2(mid = 'gray90', high = '#012a4a', midpoint = 0)
dev.off()

png('IODA-TBX2-1.png', width = 1200, height = 1000, res = 300)
plot_density(P147D100, 'TBX2', reduction = 'umap') + 
theme_void() +
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.width= unit(1, 'mm')) +
scale_color_gradient2(mid = 'gray90', high = '#012a4a', midpoint = 0)
dev.off()

png('IODA-OTOR-1.png', width = 1200, height = 1000, res = 300)
plot_density(P147D100, 'OTOR', reduction = 'umap') + 
theme_void() +
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.width= unit(1, 'mm')) +
scale_color_gradient2(mid = 'gray90', high = '#012a4a', midpoint = 0)
dev.off()

png('IODA-POU3F4-1.png', width = 1200, height = 1000, res = 300)
plot_density(P147D100, 'POU3F4', reduction = 'umap') + 
theme_void() +
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.width= unit(1, 'mm')) +
scale_color_gradient2(mid = 'gray90', high = '#012a4a', midpoint = 0)
dev.off()

png('IODA-TUBB3-1.png', width = 1200, height = 1000, res = 300)
plot_density(P147D100, 'TUBB3', reduction = 'umap') + 
theme_void() +
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.width= unit(1, 'mm')) +
scale_color_gradient2(mid = 'gray90', high = '#012a4a', midpoint = 0)
dev.off()

png('IODA-POU5F1-2.png', width = 2000, height = 500, res = 300)
Idents(P147D100) <- factor(gsub('D', '',P147D100$orig.ident), c(0, 3, 6, 8, 10, 13, 18, 24, 30, 36)) 
levels(Idents(P147D100)) <- paste0('Day ', c(0, 3, 6, 8, 10, 13, 18, 24, 30, 36))

VlnPlot(P147D100, features = c('POU5F1'), pt.size = 0) +
theme_minimal() +
theme(legend.position = 'None',
plot.title = element_text(face = 3),
panel.grid = element_blank()) +
scale_fill_manual(values = tenColors) +
xlab(NULL)
dev.off()

png('IODA-FOXH1-2.png', width = 2000, height = 500, res = 300)
Idents(P147D100) <- factor(gsub('D', '',P147D100$orig.ident), c(0, 3, 6, 8, 10, 13, 18, 24, 30, 36)) 
levels(Idents(P147D100)) <- paste0('Day ', c(0, 3, 6, 8, 10, 13, 18, 24, 30, 36))

VlnPlot(P147D100, features = c('FOXH1'), pt.size = 0) +
theme_minimal() +
theme(legend.position = 'None',
plot.title = element_text(face = 3),
panel.grid = element_blank()) +
scale_fill_manual(values = tenColors) +
xlab(NULL)
dev.off()

Idents(P147D100) <- P147D100$ID
png('IODA-BMP.png', width = 2000, height = 600, res = 300)
DotPlot(P147D100, 
features = rev(c('BMP2','BMP3','BMP4','BMP5', 'BMP7')),
dot.min = 0.01) +
coord_flip() + 
theme_minimal() +
theme(legend.position = 'bottom',
legend.key.height = unit(1, 'mm'),
plot.title = element_text(face = 3),
panel.grid = element_blank()) +
xlab(NULL) +
ylab(NULL) + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
axis.text.y = element_text(face = 3)) + 
scale_color_gradient2(low = 'gray90', mid = rgb(1,1,1,0), high = '#012a4a', midpoint = 0)
dev.off()

Idents(P147D100) <- P147D100$ID
png('IODA-HOX.png', width = 2000, height = 2800, res = 300)
DotPlot(P147D100, 
features = sort(rownames(P147D100)[grepl('^HOX', rownames(P147D100))], decreasing = TRUE),
dot.min = 0.01, scale.by = 'size') +
coord_flip() + 
theme_minimal() +
theme(legend.position = 'bottom',
legend.key.height = unit(1, 'mm'),
plot.title = element_text(face = 3),
panel.grid = element_blank()) +
xlab(NULL) +
ylab(NULL) + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
axis.text.y = element_text(face = 3)) + 
scale_color_gradient2(low = 'gray90', mid = rgb(1,1,1,0), high = '#012a4a', midpoint = 0)
dev.off()

Idents(P147D100) <- P147D100$ID
png('IODA-HOX-1.png', width = 1750, height = 2800, res = 300)
DotPlot(P147D100, 
features = rev(c('HOXA1', 'HOXB1', 'HOXD1', 'HOXA2', 'HOXB2', 
'HOXA3', 'HOXB3', 'HOXD3', 'HOXB4', 'HOXC4', 'HOXA5',
'HOXC5', 'HOXB6', 'HOXC6', 'HOXA7', 'HOXB7', 'HOXB8',
'HOXC8', 'HOXD8', 'HOXA9', 'HOXB9', 'HOXC9', 'HOXD9',
'HOXA10', 'HOXC10', 'HOXA11', 'HOXC11', 'HOXA13', 'HOXB13',
'HOXC13', 'HOXD13')),
dot.min = 0.01, col.min = 0, scale.by = 'size') +
coord_flip() + 
theme_minimal() +
theme(legend.position = 'right',
# legend.key.height = unit(1, 'mm'),
plot.title = element_text(face = 3),
panel.grid = element_blank()) +
xlab(NULL) +
ylab(NULL) + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
axis.text.y = element_text(face = 3)) + 
scale_color_gradient2(low = 'gray90', mid = rgb(1,1,1,0), high = '#012a4a', midpoint = 0)
dev.off()

png('IODA-HOX-2.png', width = 1500, height = 2800, res = 300)
Idents(P147D100) <- factor(gsub('D', '',P147D100$orig.ident), c(0, 3, 6, 8, 10, 13, 18, 24, 30, 36)) 
levels(Idents(P147D100)) <- paste0('Day ', c(0, 3, 6, 8, 10, 13, 18, 24, 30, 36))

DotPlot(P147D100, 
features = rev(c('HOXA1', 'HOXB1', 'HOXD1', 'HOXA2', 'HOXB2', 
'HOXA3', 'HOXB3', 'HOXD3', 'HOXB4', 'HOXC4', 'HOXA5',
'HOXC5', 'HOXB6', 'HOXC6', 'HOXA7', 'HOXB7', 'HOXB8',
'HOXC8', 'HOXD8', 'HOXA9', 'HOXB9', 'HOXC9', 'HOXD9',
'HOXA10', 'HOXC10', 'HOXA11', 'HOXC11', 'HOXA13', 'HOXB13',
'HOXC13', 'HOXD13')),
dot.min = 0.01, col.min = 0, scale.by = 'size') +
coord_flip() + 
theme_minimal() +
theme(legend.position = 'right',
# legend.key.height = unit(1, 'mm'),
plot.title = element_text(face = 3),
panel.grid = element_blank()) +
xlab(NULL) +
ylab(NULL) + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
axis.text.y = element_text(face = 3)) + 
scale_color_gradient2(low = 'gray90', mid = rgb(1,1,1,0), high = '#012a4a', midpoint = 0)
dev.off()

png('IODA-LR-1.png', width = 1500, height = 4000, res = 300)
Idents(P147D100) <- factor(gsub('D', '',P147D100$orig.ident), c(0, 3, 6, 8, 10, 13, 18, 24, 30, 36)) 
levels(Idents(P147D100)) <- paste0('Day ', c(0, 3, 6, 8, 10, 13, 18, 24, 30, 36))

DotPlot(P147D100, 
features = rev(c('BMP2', 'BMP3', 'BMP4', 'BMP5', 'BMP7','BMP8B',
'TGFB2', 'TGFB3','MSX1', 'MSX2', 'BAMBI', 'BMPER', 'ID1', 'ID2', 
'ID3', 'ID4', 'ZEB1', 'ZEB2', 'TWIST1', 'TWIST2', 'TGFBR2', 'FGF1',
'FGF3', 'FGF4', 'FGF7', 'FGF8', 'FGF9', 'FGF11', 'FGF13', 'FGF17',
'FGF18', 'FGF19', 'FGF20', 'SPRY1', 'SPRY2', 'DUSP6', 'ETV4', 'FOS', 
'MYC', 'JUNB', 'FGFR3', 'DHH', 'PTCH1', 'HHIP', 'DLK1', 'DLK2',
'DLL1', 'DLL3', 'JAG2', 'HEY1', 'HEY2', 'HES1', 'HES5', 'NRARP',
'NOTCH1', 'ALDH1A1', 'ALDH1A2', 'RDH10', 'CYP26A1', 'ZNF703', 
'RARG', 'WNT1', 'WNT2', 'WNT3', 'WNT3A', 'WNT4', 'WNT5A', 'WNT6', 
'WNT7B', 'WNT11', 'WNT16','LEF1', 'TCF4', 'SP5', 'CCND1', 'DKK1', 
'FZD7', 'FZD9')),
dot.min = 0.01, col.min = 0, scale.by = 'size') +
coord_flip() + 
theme_minimal() +
theme(legend.position = 'right',
# legend.key.height = unit(1, 'mm'),
plot.title = element_text(face = 3),
panel.grid = element_blank()) +
xlab(NULL) +
ylab(NULL) + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
axis.text.y = element_text(face = 3)) + 
scale_color_gradient2(low = 'gray90', mid = rgb(1,1,1,0), high = '#012a4a', midpoint = 0)
dev.off()

png('IODA-LRP-1.png', width = 1500, height = 1500, res = 300)
Idents(P147D100) <- factor(gsub('D', '',P147D100$orig.ident), c(0, 3, 6, 8, 10, 13, 18, 24, 30, 36)) 
levels(Idents(P147D100)) <- paste0('Day ', c(0, 3, 6, 8, 10, 13, 18, 24, 30, 36))

DotPlot(P147D100, 
features = rev(c('BMP2', 'BMP3', 'BMP4', 'BMP5', 'BMP7','BMP8B',
'TGFB2', 'TGFB3','MSX1', 'MSX2', 'BAMBI', 'BMPER', 'ID1', 'ID2', 
'ID3', 'ID4', 'ZEB1', 'ZEB2', 'TWIST1', 'TWIST2', 'TGFBR2')),
dot.min = 0.01, col.min = 0, scale.by = 'size') +
coord_flip() + 
theme_minimal() +
theme(legend.position = 'right',
# legend.key.height = unit(1, 'mm'),
plot.title = element_text(face = 3),
panel.grid = element_blank()) +
xlab(NULL) +
ylab(NULL) + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
axis.text.y = element_text(face = 3)) + 
scale_color_gradient2(low = 'gray90', mid = rgb(1,1,1,0), high = '#012a4a', midpoint = 0)
dev.off()

png('IODA-LRP-2.png', width = 1500, height = 1500, res = 300)
Idents(P147D100) <- factor(gsub('D', '',P147D100$orig.ident), c(0, 3, 6, 8, 10, 13, 18, 24, 30, 36)) 
levels(Idents(P147D100)) <- paste0('Day ', c(0, 3, 6, 8, 10, 13, 18, 24, 30, 36))

DotPlot(P147D100, 
features = rev(c('FGF1',
'FGF3', 'FGF4', 'FGF7', 'FGF8', 'FGF9', 'FGF11', 'FGF13', 'FGF17',
'FGF18', 'FGF19', 'FGF20', 'SPRY1', 'SPRY2', 'DUSP6', 'ETV4', 'FOS', 
'MYC', 'JUNB', 'FGFR3')),
dot.min = 0.01, col.min = 0, scale.by = 'size') +
coord_flip() + 
theme_minimal() +
theme(legend.position = 'right',
# legend.key.height = unit(1, 'mm'),
plot.title = element_text(face = 3),
panel.grid = element_blank()) +
xlab(NULL) +
ylab(NULL) + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
axis.text.y = element_text(face = 3)) + 
scale_color_gradient2(low = 'gray90', mid = rgb(1,1,1,0), high = '#012a4a', midpoint = 0)
dev.off()

png('IODA-LRP-3.png', width = 1500, height = 1000, res = 300)
Idents(P147D100) <- factor(gsub('D', '',P147D100$orig.ident), c(0, 3, 6, 8, 10, 13, 18, 24, 30, 36)) 
levels(Idents(P147D100)) <- paste0('Day ', c(0, 3, 6, 8, 10, 13, 18, 24, 30, 36))

DotPlot(P147D100, 
features = rev(c('DHH', 'PTCH1', 'HHIP')),
dot.min = 0.01, col.min = 0, scale.by = 'size') +
coord_flip() + 
theme_minimal() +
theme(legend.position = 'right',
# legend.key.height = unit(1, 'mm'),
plot.title = element_text(face = 3),
panel.grid = element_blank()) +
xlab(NULL) +
ylab(NULL) + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
axis.text.y = element_text(face = 3)) + 
scale_color_gradient2(low = 'gray90', mid = rgb(1,1,1,0), high = '#012a4a', midpoint = 0)
dev.off()

png('IODA-LRP-4.png', width = 1500, height = 1000, res = 300)
Idents(P147D100) <- factor(gsub('D', '',P147D100$orig.ident), c(0, 3, 6, 8, 10, 13, 18, 24, 30, 36)) 
levels(Idents(P147D100)) <- paste0('Day ', c(0, 3, 6, 8, 10, 13, 18, 24, 30, 36))

DotPlot(P147D100, 
features = rev(c('DLK1', 'DLK2',
'DLL1', 'DLL3', 'JAG2', 'HEY1', 'HEY2', 'HES1', 'HES5', 'NRARP',
'NOTCH1')),
dot.min = 0.01, col.min = 0, scale.by = 'size') +
coord_flip() + 
theme_minimal() +
theme(legend.position = 'right',
# legend.key.height = unit(1, 'mm'),
plot.title = element_text(face = 3),
panel.grid = element_blank()) +
xlab(NULL) +
ylab(NULL) + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
axis.text.y = element_text(face = 3)) + 
scale_color_gradient2(low = 'gray90', mid = rgb(1,1,1,0), high = '#012a4a', midpoint = 0)
dev.off()

png('IODA-LRP-5.png', width = 1500, height = 1000, res = 300)
Idents(P147D100) <- factor(gsub('D', '',P147D100$orig.ident), c(0, 3, 6, 8, 10, 13, 18, 24, 30, 36)) 
levels(Idents(P147D100)) <- paste0('Day ', c(0, 3, 6, 8, 10, 13, 18, 24, 30, 36))

DotPlot(P147D100, 
features = rev(c('ALDH1A1', 'ALDH1A2', 'RDH10', 'CYP26A1', 'ZNF703', 
'RARG')),
dot.min = 0.01, col.min = 0, scale.by = 'size') +
coord_flip() + 
theme_minimal() +
theme(legend.position = 'right',
# legend.key.height = unit(1, 'mm'),
plot.title = element_text(face = 3),
panel.grid = element_blank()) +
xlab(NULL) +
ylab(NULL) + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
axis.text.y = element_text(face = 3)) + 
scale_color_gradient2(low = 'gray90', mid = rgb(1,1,1,0), high = '#012a4a', midpoint = 0)
dev.off()

png('IODA-LRP-6.png', width = 1500, height = 1500, res = 300)
Idents(P147D100) <- factor(gsub('D', '',P147D100$orig.ident), c(0, 3, 6, 8, 10, 13, 18, 24, 30, 36)) 
levels(Idents(P147D100)) <- paste0('Day ', c(0, 3, 6, 8, 10, 13, 18, 24, 30, 36))

DotPlot(P147D100, 
features = rev(c('WNT1', 'WNT2', 'WNT3', 'WNT3A', 'WNT4', 'WNT5A', 'WNT6', 
'WNT7B', 'WNT11', 'WNT16','LEF1', 'TCF4', 'SP5', 'CCND1', 'DKK1', 
'FZD7', 'FZD9')),
dot.min = 0.01, col.min = 0, scale.by = 'size') +
coord_flip() + 
theme_minimal() +
theme(legend.position = 'right',
# legend.key.height = unit(1, 'mm'),
plot.title = element_text(face = 3),
panel.grid = element_blank()) +
xlab(NULL) +
ylab(NULL) + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
axis.text.y = element_text(face = 3)) + 
scale_color_gradient2(low = 'gray90', mid = rgb(1,1,1,0), high = '#012a4a', midpoint = 0)
dev.off()

png('IODA-LR-2.png', height = 1000, width = 4500, res = 300)
Idents(P147D100) <- factor(gsub('D', '',P147D100$orig.ident), rev(c(0, 3, 6, 8, 10, 13, 18, 24, 30, 36)))
levels(Idents(P147D100)) <- paste0('Day ', rev(c(0, 3, 6, 8, 10, 13, 18, 24, 30, 36)))

DotPlot(P147D100, 
features = rev(c('BMP2', 'BMP3', 'BMP4', 'BMP5', 'BMP7','BMP8B',
'TGFB2', 'TGFB3','MSX1', 'MSX2', 'BAMBI', 'BMPER', 'ID1', 'ID2', 
'ID3', 'ID4', 'ZEB1', 'ZEB2', 'TWIST1', 'TWIST2', 'TGFBR2', 'FGF1',
'FGF3', 'FGF4', 'FGF7', 'FGF8', 'FGF9', 'FGF11', 'FGF13', 'FGF17',
'FGF18', 'FGF19', 'FGF20', 'SPRY1', 'SPRY2', 'DUSP6', 'ETV4', 'FOS', 
'MYC', 'JUNB', 'FGFR3', 'DHH', 'PTCH1', 'HHIP', 'DLK1', 'DLK2',
'DLL1', 'DLL3', 'JAG2', 'HEY1', 'HEY2', 'HES1', 'HES5', 'NRARP',
'NOTCH1', 'ALDH1A1', 'ALDH1A2', 'RDH10', 'CYP26A1', 'ZNF703', 
'RARG', 'WNT1', 'WNT2', 'WNT3', 'WNT3A', 'WNT4', 'WNT5A', 'WNT6', 
'WNT7B', 'WNT11', 'WNT16','LEF1', 'TCF4', 'SP5', 'CCND1', 'DKK1', 
'FZD7', 'FZD9')),
dot.min = 0.01, col.min = 0, scale.by = 'size') +
theme_minimal() +
theme(legend.position = 'right',
# legend.key.height = unit(1, 'mm'),
plot.title = element_text(face = 3),
panel.grid = element_blank()) +
xlab(NULL) +
ylab(NULL) + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,face = 3)) + 
scale_color_gradient2(low = 'gray90', mid = rgb(1,1,1,0), high = '#012a4a', midpoint = 0)
dev.off()

png('IODA-LR-3.png', width = 2000, height = 4000, res = 300)
Idents(P147D100) <- P147D100$ID

DotPlot(P147D100, 
features = rev(c('BMP2', 'BMP3', 'BMP4', 'BMP5', 'BMP7','BMP8B',
'TGFB2', 'TGFB3','MSX1', 'MSX2', 'BAMBI', 'BMPER', 'ID1', 'ID2', 
'ID3', 'ID4', 'ZEB1', 'ZEB2', 'TWIST1', 'TWIST2', 'TGFBR2', 'FGF1',
'FGF3', 'FGF4', 'FGF7', 'FGF8', 'FGF9', 'FGF11', 'FGF13', 'FGF17',
'FGF18', 'FGF19', 'FGF20', 'SPRY1', 'SPRY2', 'DUSP6', 'ETV4', 'FOS', 
'MYC', 'JUNB', 'FGFR3', 'DHH', 'PTCH1', 'HHIP', 'DLK1', 'DLK2',
'DLL1', 'DLL3', 'JAG2', 'HEY1', 'HEY2', 'HES1', 'HES5', 'NRARP',
'NOTCH1', 'ALDH1A1', 'ALDH1A2', 'RDH10', 'CYP26A1', 'ZNF703', 
'RARG', 'WNT1', 'WNT2', 'WNT3', 'WNT3A', 'WNT4', 'WNT5A', 'WNT6', 
'WNT7B', 'WNT11', 'WNT16','LEF1', 'TCF4', 'SP5', 'CCND1', 'DKK1', 
'FZD7', 'FZD9')),
dot.min = 0.01, col.min = 0, scale.by = 'size') +
coord_flip() + 
theme_minimal() +
theme(legend.position = 'right',
# legend.key.height = unit(1, 'mm'),
plot.title = element_text(face = 3),
panel.grid = element_blank()) +
xlab(NULL) +
ylab(NULL) + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
axis.text.y = element_text(face = 3)) + 
scale_color_gradient2(low = 'gray90', mid = rgb(1,1,1,0), high = '#012a4a', midpoint = 0)
dev.off()

png('IODA-LRP-7.png', width = 2000, height = 1500, res = 300)
Idents(P147D100) <- P147D100$ID

DotPlot(P147D100, 
features = rev(c('BMP2', 'BMP3', 'BMP4', 'BMP5', 'BMP7','BMP8B',
'TGFB2', 'TGFB3','MSX1', 'MSX2', 'BAMBI', 'BMPER', 'ID1', 'ID2', 
'ID3', 'ID4', 'ZEB1', 'ZEB2', 'TWIST1', 'TWIST2', 'TGFBR2')),
dot.min = 0.01, col.min = 0, scale.by = 'size') +
coord_flip() + 
theme_minimal() +
theme(legend.position = 'right',
# legend.key.height = unit(1, 'mm'),
plot.title = element_text(face = 3),
panel.grid = element_blank()) +
xlab(NULL) +
ylab(NULL) + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
axis.text.y = element_text(face = 3)) + 
scale_color_gradient2(low = 'gray90', mid = rgb(1,1,1,0), high = '#012a4a', midpoint = 0)
dev.off()

png('IODA-LRP-8.png', width = 2000, height = 1500, res = 300)
Idents(P147D100) <- P147D100$ID

DotPlot(P147D100, 
features = rev(c('FGF1',
'FGF3', 'FGF4', 'FGF7', 'FGF8', 'FGF9', 'FGF11', 'FGF13', 'FGF17',
'FGF18', 'FGF19', 'FGF20', 'SPRY1', 'SPRY2', 'DUSP6', 'ETV4', 'FOS', 
'MYC', 'JUNB', 'FGFR3')),
dot.min = 0.01, col.min = 0, scale.by = 'size') +
coord_flip() + 
theme_minimal() +
theme(legend.position = 'right',
# legend.key.height = unit(1, 'mm'),
plot.title = element_text(face = 3),
panel.grid = element_blank()) +
xlab(NULL) +
ylab(NULL) + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
axis.text.y = element_text(face = 3)) + 
scale_color_gradient2(low = 'gray90', mid = rgb(1,1,1,0), high = '#012a4a', midpoint = 0)
dev.off()

png('IODA-LRP-9.png', width = 2000, height = 1000, res = 300)
Idents(P147D100) <- P147D100$ID

DotPlot(P147D100, 
features = rev(c('DHH', 'PTCH1', 'HHIP')),
dot.min = 0.01, col.min = 0, scale.by = 'size') +
coord_flip() + 
theme_minimal() +
theme(legend.position = 'right',
# legend.key.height = unit(1, 'mm'),
plot.title = element_text(face = 3),
panel.grid = element_blank()) +
xlab(NULL) +
ylab(NULL) + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
axis.text.y = element_text(face = 3)) + 
scale_color_gradient2(low = 'gray90', mid = rgb(1,1,1,0), high = '#012a4a', midpoint = 0)
dev.off()

png('IODA-LRP-10.png', width = 2000, height = 1000, res = 300)
Idents(P147D100) <- P147D100$ID

DotPlot(P147D100, 
features = rev(c('DLK1', 'DLK2',
'DLL1', 'DLL3', 'JAG2', 'HEY1', 'HEY2', 'HES1', 'HES5', 'NRARP',
'NOTCH1')),
dot.min = 0.01, col.min = 0, scale.by = 'size') +
coord_flip() + 
theme_minimal() +
theme(legend.position = 'right',
# legend.key.height = unit(1, 'mm'),
plot.title = element_text(face = 3),
panel.grid = element_blank()) +
xlab(NULL) +
ylab(NULL) + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
axis.text.y = element_text(face = 3)) + 
scale_color_gradient2(low = 'gray90', mid = rgb(1,1,1,0), high = '#012a4a', midpoint = 0)
dev.off()

png('IODA-LRP-11.png', width = 2000, height = 1000, res = 300)
Idents(P147D100) <- P147D100$ID

DotPlot(P147D100, 
features = rev(c('ALDH1A1', 'ALDH1A2', 'RDH10', 'CYP26A1', 'ZNF703', 
'RARG')),
dot.min = 0.01, col.min = 0, scale.by = 'size') +
coord_flip() + 
theme_minimal() +
theme(legend.position = 'right',
# legend.key.height = unit(1, 'mm'),
plot.title = element_text(face = 3),
panel.grid = element_blank()) +
xlab(NULL) +
ylab(NULL) + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
axis.text.y = element_text(face = 3)) + 
scale_color_gradient2(low = 'gray90', mid = rgb(1,1,1,0), high = '#012a4a', midpoint = 0)
dev.off()

png('IODA-LRP-12.png', width = 2000, height = 1500, res = 300)
Idents(P147D100) <- P147D100$ID

DotPlot(P147D100, 
features = rev(c('WNT1', 'WNT2', 'WNT3', 'WNT3A', 'WNT4', 'WNT5A', 'WNT6', 
'WNT7B', 'WNT11', 'WNT16','LEF1', 'TCF4', 'SP5', 'CCND1', 'DKK1', 
'FZD7', 'FZD9')),
dot.min = 0.01, col.min = 0, scale.by = 'size') +
coord_flip() + 
theme_minimal() +
theme(legend.position = 'right',
# legend.key.height = unit(1, 'mm'),
plot.title = element_text(face = 3),
panel.grid = element_blank()) +
xlab(NULL) +
ylab(NULL) + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
axis.text.y = element_text(face = 3)) + 
scale_color_gradient2(low = 'gray90', mid = rgb(1,1,1,0), high = '#012a4a', midpoint = 0)
dev.off()


png('IODA-LR-4.png', height = 1500, width = 4500, res = 300)
Idents(P147D100) <- P147D100$ID

DotPlot(P147D100, 
features = rev(c('BMP2', 'BMP3', 'BMP4', 'BMP5', 'BMP7','BMP8B',
'TGFB2', 'TGFB3','MSX1', 'MSX2', 'BAMBI', 'BMPER', 'ID1', 'ID2', 
'ID3', 'ID4', 'ZEB1', 'ZEB2', 'TWIST1', 'TWIST2', 'TGFBR2', 'FGF1',
'FGF3', 'FGF4', 'FGF7', 'FGF8', 'FGF9', 'FGF11', 'FGF13', 'FGF17',
'FGF18', 'FGF19', 'FGF20', 'SPRY1', 'SPRY2', 'DUSP6', 'ETV4', 'FOS', 
'MYC', 'JUNB', 'FGFR3', 'DHH', 'PTCH1', 'HHIP', 'DLK1', 'DLK2',
'DLL1', 'DLL3', 'JAG2', 'HEY1', 'HEY2', 'HES1', 'HES5', 'NRARP',
'NOTCH1', 'ALDH1A1', 'ALDH1A2', 'RDH10', 'CYP26A1', 'ZNF703', 
'RARG', 'WNT1', 'WNT2', 'WNT3', 'WNT3A', 'WNT4', 'WNT5A', 'WNT6', 
'WNT7B', 'WNT11', 'WNT16','LEF1', 'TCF4', 'SP5', 'CCND1', 'DKK1', 
'FZD7', 'FZD9')),
dot.min = 0.01, col.min = 0, scale.by = 'size') +
theme_minimal() +
theme(legend.position = 'right',
# legend.key.height = unit(1, 'mm'),
plot.title = element_text(face = 3),
panel.grid = element_blank()) +
xlab(NULL) +
ylab(NULL) + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,face = 3)) + 
scale_color_gradient2(low = 'gray90', mid = rgb(1,1,1,0), high = '#012a4a', midpoint = 0)
dev.off()


png('IODA-CICI-2.png', width = 2500, height = 750, res = 300)
Idents(P147D100) <- factor(gsub('D', '',P147D100$orig.ident), c(0, 3, 6, 8, 10, 13, 18, 24, 30, 36)) 
levels(Idents(P147D100)) <- paste0('Day ', c(0, 3, 6, 8, 10, 13, 18, 24, 30, 36))

DotPlot(P147D100, 
features = rev(c('GATA3', 'OTX2')),
dot.min = 0.01, col.min = 0, scale.by = 'size') +
coord_flip() + 
theme_minimal() +
theme(legend.position = 'bottom',
legend.key.height = unit(1, 'mm'),
plot.title = element_text(face = 3),
panel.grid = element_blank()) +
xlab(NULL) +
ylab(NULL) + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
axis.text.y = element_text(face = 3)) + 
scale_color_gradient2(low = 'gray90', mid = rgb(1,1,1,0), high = '#012a4a', midpoint = 0)
dev.off()

sapply(c('HOXA1', 'HOXB1', 'HOXD1', 'HOXA2', 'HOXB2', 
'HOXA3', 'HOXB3', 'HOXD3', 'HOXB4', 'HOXC4', 'HOXA5',
'HOXC5', 'HOXB6', 'HOXC6', 'HOXA7', 'HOXB7', 'HOXB8',
'HOXC8', 'HOXD8', 'HOXA9', 'HOXB9', 'HOXC9', 'HOXD9',
'HOXA10', 'HOXC10', 'HOXA11', 'HOXC11', 'HOXA13', 'HOXB13',
'HOXC13', 'HOXD13', 'GATA3', 'OTX2', 'NKX2-5', 'VGLL2'), function(id){
      png(paste0('IODA-',id,'-1.png'), width = 1200, height = 1000, res = 300)
        P <- plot_density(P147D100, id, reduction = 'umap') + 
        theme_void() +
        theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(face = 3),
        legend.key.width= unit(1, 'mm')) +
        scale_color_gradient2(mid = 'gray90', high = '#012a4a', midpoint = 0)
        print(P)
        dev.off()  
})


Idents(P147D100) <- factor(gsub('D', '',P147D100$orig.ident), c(0, 3, 6, 8, 10, 13, 18, 24, 30, 36)) 
levels(Idents(P147D100)) <- paste0('Day ', c(0, 3, 6, 8, 10, 13, 18, 24, 30, 36))

plotPathwayAvg <- function(geneList, plotTitle){
        E <- colMeans(P147D100@assays$RNA@data[c(geneList),])
        df <- data.frame(E = E, G = Idents(P147D100))
        P <- ggplot(df, aes(G,E, fill = G)) +
        geom_violin() +
        geom_boxplot(width = 0.05, fill = 'gray90',   outlier.color = NA) +
        theme_minimal() +
        theme(legend.position = 'None',
        plot.subtitle = element_text(face = 3),
        panel.grid = element_blank()) +
        scale_fill_manual(values = tenColors) +
        xlab(NULL) +
        ylab('Average Expression') +
        labs(title = plotTitle, subtitle = paste0(geneList, collapse = ', '))
        return(P)
}

png('IODA-BMPPathway.png', width = 2000, height = 600, res = 300)
plotPathwayAvg(
        geneList = c('ID1', 'ID2', 'ID3', 'ID4', 'MSX1', 'MSX2', 'BAMBI', 'BMPER'),
        plotTitle = 'BMP Target Genes')
dev.off()

png('IODA-FGFPathway.png', width = 2000, height = 600, res = 300)
plotPathwayAvg(
        geneList = c('SPRY1', 'SPRY2', 'DUSP6', 'ETV4'),
        plotTitle = 'FGF Target Genes')
dev.off()

png('IODA-WNTPathway.png', width = 2000, height = 600, res = 300)
plotPathwayAvg(
        geneList = c('LEF1', 'TCF4', 'SP5', 'CCND1', 'DKK1'),
        plotTitle = 'WNT Target Genes')
dev.off()

Idents(P147D100) <- P147D100$ID
load('/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler/P147D100_Annotated.RData')
F5Data <- P147D100

png('F5-M.png', height = 500*2, width = 500*5, res = 300)
F5PA1 <- plot_density(F5Data, c('PHOX2B'), reduction = 'umap', size = 0.01)
F5PA2 <- plot_density(F5Data, c('NEUROG2'), reduction = 'umap', size = 0.01)
F5PB1 <- plot_density(F5Data, c('PRDM12'), reduction = 'umap', size = 0.01)
F5PB2 <- plot_density(F5Data, c('ADAMTS6'), reduction = 'umap', size = 0.01)
F5PC1 <- plot_density(F5Data, c('GPM6A'), reduction = 'umap', size = 0.01)
F5PC2 <- plot_density(F5Data, c('POU3F2'), reduction = 'umap', size = 0.01)
F5PD1 <- plot_density(F5Data, c('FGF8'), reduction = 'umap', size = 0.01)
F5PD2 <- plot_density(F5Data, c('PAX2'), reduction = 'umap', size = 0.01)
F5PE1 <- plot_density(F5Data, c('S100A1'), reduction = 'umap', size = 0.01)
F5PE2 <- plot_density(F5Data, c('S100A13'), reduction = 'umap', size = 0.01)

((F5PA1/F5PA2)|(F5PB1/F5PB2)|(F5PC1/F5PC2)|(F5PD1/F5PD2)|(F5PE1/F5PE2))&
theme_void() &
xlim(c(3,11)) & 
ylim(c(-11,-4)) & 
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.height= unit(1, 'mm'), 
legend.position = 'None') &
scale_color_gradient2(mid = 'gray90', high = '#F94144', midpoint = 0)
dev.off()

F5PA <- plot_density(F5Data, c('PHOX2B', 'NEUROG2'), reduction = 'umap', size = 0.01, joint = TRUE)[[3]] + labs(subtitle = 'Epibranchial Neurons')
F5PB <- plot_density(F5Data, c('PRDM12', 'ADAMTS6'), reduction = 'umap', size = 0.01, joint = TRUE)[[3]] + labs(subtitle = 'NC Neurons')
F5PC <- plot_density(F5Data, c('GPM6A', 'POU3F2'), reduction = 'umap', size = 0.01, joint = TRUE)[[3]] + labs(subtitle = 'CNS Neurons')
F5PD <- plot_density(F5Data, c('FGF8', 'PAX2'), reduction = 'umap', size = 0.01, joint = TRUE)[[3]] + labs(subtitle = 'Immature Otic Neurons')
F5PE <- plot_density(F5Data, c('S100A1', 'S100A13'), reduction = 'umap', size = 0.01, joint = TRUE)[[3]] + labs(subtitle = 'Otic Neurons')

png('F5-M2.png', width = 750*5, height = 850, res = 300)
(F5PA | F5PB | F5PC | F5PD | F5PE)&
theme_void() &
xlim(c(3,11)) & 
ylim(c(-11,-4)) & 
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.height= unit(1, 'mm'), 
legend.position = 'None') &
scale_color_gradient2(mid = 'gray90', high = '#F94144', midpoint = 0)
dev.off()

png('F5-M3.png', width = 750*2, height = 750*2, res = 300)
plot_density(F5Data, c('NEUROG1', 'NEUROD1', 'ISL1', 'PRPH'), reduction = 'umap', size = 0.01) &
theme_void() &
xlim(c(3,11)) & 
ylim(c(-11,-4)) & 
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.height= unit(1, 'mm'), 
legend.position = 'None') &
scale_color_gradient2(mid = 'gray90', high = '#F94144', midpoint = 0)
dev.off()


png('F5-M4.png', width = 750*1, height = 750*1, res = 300)
plot_density(F5Data, c('PAX3'), reduction = 'umap', size = 0.01) &
theme_void() &
# xlim(c(3,11)) & 
# ylim(c(-11,-4)) & 
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.height= unit(1, 'mm'), 
legend.position = 'None') &
scale_color_gradient2(mid = 'gray90', high = '#F94144', midpoint = 0)
dev.off()


# Writting UMAP
png('results/annotatedUMAP.png', width = 2000, height = 2000, res = 300)
PA1 <- UMAPPlot(P147D100) + 
theme_light() +
theme(legend.position = 'None', 
        panel.grid = element_blank(),
        panel.border = element_blank()) +
xlab('UMAP 1') +
ylab('UMAP 2') +
scale_color_manual(values = ctColors) +
geom_text_repel(data = LabelPosition, 
mapping = aes(X,Y, label = Label), 
bg.color = 'white', 
box.padding = 0,   
min.segment.length = 0)
print(PA1)
dev.off()

CP <- table(P147D100$orig.ident, P147D100$CT)
# CP <- CP[,levels(Idents(P147D100))]
CP <- (CP/rowSums(CP)) * 100
CP <- reshape2::melt(CP)
colnames(CP) <- c('TimePoint', 'CellType', 'Proportion')
CP$TimePoint <-  gsub('D', 'Day ', CP$TimePoint)
CP$TimePoint <- factor(CP$TimePoint, rev(c('Day 0', 'Day 3', 
'Day 6', 'Day 8', 'Day 10', 'Day 13', 'Day 18', 'Day 24', 
'Day 30', 'Day 36')))
CP$CellType <- factor(CP$CellType, Label1)

png('results/sampleProportions.png', width = 4000, height = 1800, res = 300)
PA2 <- ggplot(CP, aes(Proportion, TimePoint, fill = CellType)) +
geom_bar(stat = 'identity') +
scale_fill_manual(values = ctColors) +
theme_light() +
theme(panel.border = element_blank(), panel.grid = element_blank()) +
guides(fill=guide_legend(ncol=1)) +
xlab('Proportion (%)') +
ylab('Time Point') +
labs(fill = 'Cell Type')
print(PA2)
dev.off()

png('results/sampleProportionsVOID.png', width = 4000, height = 1800, res = 300)
PA2 <- ggplot(CP, aes(Proportion, TimePoint, fill = CellType)) +
geom_bar(stat = 'identity') +
scale_fill_manual(values = ctColors) +
theme_void() +
theme(panel.border = element_blank(), 
panel.grid = element_blank(),
legend.position = 'None') +
guides(fill=guide_legend(ncol=1)) +
xlab('Proportion (%)') +
ylab('Time Point') +
labs(fill = 'Cell Type')
print(PA2)
dev.off()


png('PA.png', width = 5000, height = 1800, res = 300)
PA1 + PA2
dev.off()

PC <- plot_density(P147D100, 
c('POU5F1', 'OC90', 'NEFM', 'TFAP2C'), 
reduction = 'umap', size = 0.05) &
theme_light() +
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.width= unit(1, 'mm')) &
scale_color_gradient2(mid = 'gray90', high = '#012a4a', midpoint = 0)

png('PC.png', width = 2500*0.8, height = 2000*0.7, res = 300)
print(PC)
dev.off()

PC <- plot_density(P147D100, 
c('POU5F1', 'OC90', 'NEFM', 'TFAP2C',
'PRRX1', 'CDH1', 'EPCAM', 'TFAP2A', 
'GATA2', 'GATA3', 'DLX3', 'DLX5', 'MYO7A',
'OTX2', 'OTX1', 'GBX2', 'PAX8', 'SIX4', 
'EYA1', 'SIX1'), 
reduction = 'umap', size = 0.05) &
theme_void() +
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.width= unit(1, 'mm')) &
scale_color_gradient2(mid = 'gray90', high = '#012a4a', midpoint = 0)

png('PC2.png', width = 750*5, height = 600*4, res = 300)
print(PC)
dev.off()

Idents(P147D100) <- factor(Idents(P147D100), (c('PSC1', 'PSC2', 'NE', 'cSE', 'ctSE', 'cySE', 'NCC', 'aPPE', 'pPPE', 
'bEPI', 'NEU', 'OEPD', 'OEP', 'cME', 'fME', 'FIB', 'MEL', 'spEPI', 'MYO', 'NEP', 'HC')))
PD <- doMarkersPlot(P147D100, n = 5)
PD$data$id <- factor(PD$data$id, rev(c('PSC1', 'PSC2', 'NE', 'cSE', 'ctSE', 'cySE', 'NCC', 'aPPE', 'pPPE', 
'bEPI', 'NEU', 'OEPD', 'OEP', 'cME', 'fME', 'FIB', 'MEL', 'spEPI', 'MYO', 'NEP', 'HC')))
PD <- PD + 
scale_color_gradient2(mid = rgb(1,1,1,0), high = '#012a4a') +
labs(fill = 'A', color = 'B') +
theme(panel.border = element_blank(), panel.grid = element_blank()) +
xlab('Genes') + ylab('Cell Types')

png('PD2.png', width = 6250, height = 1500, res = 300)
print(PD)
dev.off()

Idents(P147D100) <- factor(Idents(P147D100), sort(levels(Idents(P147D100)), decreasing = FALSE))
PD <- doMarkersPlot(P147D100, n = 5)
PD <- PD + 
scale_color_gradient2(mid = rgb(1,1,1,0), high = '#012a4a') +
labs(fill = 'A', color = 'B') +
theme(panel.border = element_blank(), panel.grid = element_blank()) +
xlab('Genes') + ylab('Cell Types')

png('PD.png', width = 4000, height = 1300, res = 300)
print(PD)
dev.off()

Idents(P147D100) <- gsub('D', 'Day ', P147D100$orig.ident)
Idents(P147D100) <- factor(Idents(P147D100), (c('Day 0', 'Day 3', 
'Day 6', 'Day 8', 'Day 10', 'Day 13', 'Day 18', 'Day 24', 
'Day 30', 'Day 36')))
png('PE.png', width = 2400*0.7, height = 2000*0.7, res = 300)
PE <- UMAPPlot(P147D100) + theme_void() +
theme(panel.border = element_blank(), panel.grid = element_blank()) +
scale_color_manual(values = tenColors) +
xlab('UMAP 1') +
ylab('UMAP 2')
print(PE)
dev.off()

png('LP.png', width = 5000, height = 1500, res = 300)
PD + PE + plot_layout(widths = c(2.5,1))
dev.off()


Idents(P147D100) <- factor(P147D100$ID, Label2)


D <- P147D100@reductions$umap@cell.embeddings
D <- as.data.frame.array(D)
D$ID <- P147D100$ID
D$Color <- matchColors2[D$ID]
D$TP <- P147D100$orig.ident
D$Color[D$TP != 'D0'] <- 'gray95'
D <- D[order(D$Color, decreasing = TRUE),]
P2BA <- ggplot(D, aes(UMAP_1, UMAP_2)) +
geom_point(color = D$Color, size = 0.01) +
theme_light() +
theme(legend.position = 'None', 
        panel.grid = element_blank(),
        panel.border = element_blank()) +
xlab('UMAP 1') +
ylab('UMAP 2') +
scale_color_manual(values = matchColors2[levels(Idents(P147D100[,P147D100$orig.ident == 'D0']))]) +
geom_text_repel(data = LabelPosition[LabelPosition$Label %in% c('PSC1', 'PSC2'),], 
mapping = aes(X,Y, label = Label), 
bg.color = 'white', 
box.padding = 0,   
min.segment.length = 0) +
labs(title = 'Day 0')
D <- P147D100@reductions$umap@cell.embeddings
D <- as.data.frame.array(D)
D$ID <- P147D100$ID
D$Color <- matchColors2[D$ID]
D$TP <- P147D100$orig.ident
D$Color[D$TP != 'D3'] <- 'gray95'
D <- D[order(D$Color, decreasing = TRUE),]
P2BB <- ggplot(D, aes(UMAP_1, UMAP_2)) +
geom_point(color = D$Color, size = 0.01) +
theme_light() +
theme(legend.position = 'None', 
        panel.grid = element_blank(),
        panel.border = element_blank()) +
xlab('UMAP 1') +
ylab('UMAP 2') +
scale_color_manual(values = matchColors2[levels(Idents(P147D100[,P147D100$orig.ident == 'D3']))]) +
geom_text_repel(data = LabelPosition[LabelPosition$Label %in% c('ctSE', 'cSE', 'NE'),], 
mapping = aes(X,Y, label = Label), 
bg.color = 'white', 
box.padding = 0,   
min.segment.length = 0) +
labs(title = 'Day 3')
P2B <- P2BA + P2BB

png('P2B.png', width = 2000, height = 900, res = 300)
print(P2B)
dev.off()

P2DData <- P147D100[,P147D100$orig.ident == 'D3' & P147D100$ID %in% c('NE', 'cSE')]

png('P2D.png', width = 2500, height = 500, res = 300)
VlnPlot(P2DData, 
features = c('TFAP2A', 'KRT18', 'KRT8', 'CDX2', 'HAND1', 'CDH1'), 
pt.size = 0, ncol = 6) &
theme_minimal() &
theme(panel.grid = element_blank(), legend.position = 'None') &
xlab(NULL) & 
scale_fill_manual(values = c('#9A3C09','#9B5805')) &
theme(plot.title = element_text(face = 3))
dev.off()

png('P2E.png', width = 2500, height = 500, res = 300)
VlnPlot(P2DData, 
features = c('PAX6', 'WNT1', 'NES', 'ZIC2', 'ZIC5', 'CDH2'), 
pt.size = 0, ncol = 6) &
theme_minimal() &
theme(panel.grid = element_blank(), legend.position = 'None') &
xlab(NULL) & 
scale_fill_manual(values = c('#9A3C09','#9B5805')) &
theme(plot.title = element_text(face = 3))
dev.off()


png('P2F.png', width = 2500, height = 400, res = 300)
DotPlot(P2DData, 
features = c('GABRP', 'VTCN1', 'S100A10', 'GADD45G', 'ACTC1', 'IGFBP3')) +
theme_light() +
theme(panel.grid = element_blank(),
panel.border = element_blank(), 
legend.position = 'bottom', 
legend.key.height = unit(1,'mm'),
axis.text.x = element_text(face = 3)) +
xlab(NULL) +
ylab(NULL) +
scale_color_gradient2(low = 'gray90', mid = rgb(1,1,1,0), high = '#9B5805', midpoint = 0)
dev.off()

# Requested on 04/03
P2DData <- P147D100[,P147D100$orig.ident == 'D3' & P147D100$ID %in% c('NE', 'cSE', 'ctSE')]

png('P2F-1.png', width = 2500, height = 500, res = 300)
DotPlot(P2DData,
features = c('GABRP', 'VTCN1', 'S100A10', 'GADD45G', 'ACTC1', 'IGFBP3')) +
theme_light() +
theme(panel.grid = element_blank(),
panel.border = element_blank(), 
legend.position = 'bottom', 
legend.key.height = unit(1,'mm'),
axis.text.x = element_text(face = 3)) +
xlab(NULL) +
ylab(NULL) +
scale_color_gradient2(low = rgb(1,1,1,0), mid = rgb(155/255,88/255,5/255,.25), high = '#9B5805', midpoint = 0)
dev.off()


# 
# CM <- doMarkersPlot(P147D100)

# png('t1.png', width =5000, height = 2000, res = 300)
# CM + xlab('Cell Type')
# dev.off()


png('P3C.png', width = 3000, height = 1500, res = 300)
D <- P147D100@reductions$umap@cell.embeddings
D <- as.data.frame.array(D)
D$ID <- P147D100$ID
D$Color <- matchColors2[D$ID]
D$TP <- P147D100$orig.ident
D$Color[D$TP != 'D6'] <- 'gray95'
D <- D[order(D$Color, decreasing = TRUE),]
sLabels <- table(D$ID[D$Color != 'gray95'])
sLabels <- names(sLabels[sLabels >= 50])
UMAPD6 <- ggplot(D, aes(UMAP_1, UMAP_2)) +
geom_point(color = D$Color, size = 0.01) +
theme_light() +
theme(legend.position = 'None', 
        panel.grid = element_blank(),
        panel.border = element_blank()) +
xlab('UMAP 1') +
ylab('UMAP 2') +
scale_color_manual(values = matchColors2[levels(Idents(P147D100[,P147D100$orig.ident == 'D6']))]) +
geom_text_repel(data = LabelPosition[LabelPosition$Label %in% sLabels,], 
mapping = aes(X,Y, label = Label), 
bg.color = 'white', 
box.padding = 0,   
min.segment.length = 0) +
labs(title = 'Day 6')

D <- P147D100@reductions$umap@cell.embeddings
D <- as.data.frame.array(D)
D$ID <- P147D100$ID
D$Color <- matchColors2[D$ID]
D$TP <- P147D100$orig.ident
D$Color[D$TP != 'D8'] <- 'gray95'
D <- D[order(D$Color, decreasing = TRUE),]
sLabels <- table(D$ID[D$Color != 'gray95'])
sLabels <- names(sLabels[sLabels >= 50])
UMAPD8 <- ggplot(D, aes(UMAP_1, UMAP_2)) +
geom_point(color = D$Color, size = 0.01) +
theme_light() +
theme(legend.position = 'None', 
        panel.grid = element_blank(),
        panel.border = element_blank()) +
xlab('UMAP 1') +
ylab('UMAP 2') +
scale_color_manual(values = matchColors2[levels(Idents(P147D100[,P147D100$orig.ident == 'D8']))]) +
geom_text_repel(data = LabelPosition[LabelPosition$Label %in% sLabels,], 
mapping = aes(X,Y, label = Label), 
bg.color = 'white', 
box.padding = 0,   
min.segment.length = 0) +
labs(title = 'Day 8')

UMAPD6 + UMAPD8
dev.off()

CP <- table(P147D100$orig.ident, P147D100$ID)
CP <- round((CP/rowSums(CP)) * 100,1)
CP <- CP[c('D6', 'D8'),]
CP <- reshape2::melt(CP)
CP <- CP[CP$value > 0,]
colnames(CP) <- c('Day', 'CT', 'Proportion')
CP$CT <- factor(as.vector(CP$CT), unique(as.vector(CP$CT)))

png('P3D.png', width = 3000, height = 600, res = 300)
ggplot(CP, aes(Proportion, Day, fill = CT)) +
geom_bar(stat = 'identity') +
theme_light() +
theme(legend.position = 'bottom', 
        panel.grid = element_blank(),
        panel.border = element_blank()) +
        scale_fill_manual(values = matchColors2[levels(CP$CT)]) +
        guides(fill=guide_legend(nrow=2)) +
        labs(fill = 'Cell Type')
dev.off()


# load('P147_scVeloResults.RData')

# library(scater)
# reducedDim(P147ALL, "UMAP") <- E

# png('AAAA.png', width = 1500, height = 1500, res = 300)
# M <- velociraptor::plotVelocityStream(
#     P147ALL, 
#     P147EM, 
#     color_by = matchColors2[levels(Idents(P147D100))],
#     use.dimred = 2, 
#     grid.resolution = 90,
#     scale = TRUE,
#     stream.L = 1,
#     stream.res = 100,
#     stream.width = 1,
#     arrow.length = 0.2,
#     arrow.angle = 25
#     )
# print(M)
# dev.off()


load('P147_scVeloResults.RData')
library(dplyr)
arrows <- P147G2
colnames(arrows) <- c('x0', 'y0', 'x1', 'y1')
umap_arrows <- arrows %>%
    as.data.frame() %>%
    mutate(x2 = x0 + (x1 - x0) * .5,
           y2 = y0 + (y1 - y0) * .5)

png('V.png', width = 1800, height = 1800, res = 300)
UMAPPlot(P147D100, label = FALSE) +
scale_color_manual(values = matchColors2[levels(Idents(P147D100))]) +
theme_light() +
theme(legend.position = 'None', 
        panel.grid = element_blank(),
        panel.border = element_blank()) +
        geom_curve(data = umap_arrows, curvature = 0.2, angle = 45,
                 aes(x = x0, xend = x2, y = y0, yend = y2),
                 size = .3,
                 arrow = arrow(length = unit(1, "points"), type = "closed"),
                 colour = "grey20", alpha = 0.8) +
                 xlab('UMAP 1') +
ylab('UMAP 2') +
geom_text_repel(data = LabelPosition, 
mapping = aes(X,Y, label = Label), 
bg.color = 'white', 
box.padding = 0,   
min.segment.length = 0)
dev.off()

D6D8 <- P147D100[,P147D100$orig.ident %in% c('D6', 'D8')]
SS1 <- D6D8[,D6D8$ID %in% c('cSE', 'ctSE', 'bEPI', 'spEPI', 'aPPE', 'pPPE', 'OEPD')]
SS1 <- FindVariableFeatures(SS1)
SS1 <- RunPCA(SS1)
SS1 <- RunHarmony(SS1, 'orig.ident')
SS1 <- RunUMAP(SS1, reduction = 'harmony', dims = 1:30)
save(SS1, file = 'SS1.RData')


SS2 <- D6D8[,D6D8$ID %in% c('NE', 'NEP', 'NEU', 'NCC')]
SS2 <- FindVariableFeatures(SS2)
SS2 <- RunPCA(SS2)
SS2 <- RunHarmony(SS2, 'orig.ident')
SS2 <- RunUMAP(SS2, reduction = 'harmony', dims = 1:3)
save(SS2, file = 'SS2.RData')

# png('P3F.png', width = 3000, height = 1250, res = 300)
# plot_density(P147D100[,P147D100$orig.ident %in% c('D6', 'D8')], c('OTX1', 'GBX2'), reduction = 'umap', size = 0.05) &
# theme_light() &
# theme(panel.grid = element_blank(), 
# panel.border = element_blank(),
# plot.title = element_text(face = 3),
# legend.key.width= unit(1, 'mm')) &
# scale_color_gradient2(mid = 'gray90', high = '#012a4a', midpoint = 0)
# dev.off()
load('SS1.RData')

png('P3F.png', width = 1000*2*0.8, height = 1000*0.8, res = 300)
plot_density(SS1, c('OTX1', 'GBX2'), reduction = 'umap', size = 0.05) &
theme_light() &
theme(legend.position = 'None', 
panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.width= unit(1, 'mm')) &
scale_color_gradient2(mid = 'gray90', high = '#012a4a', midpoint = 0)
dev.off()

png('P3F-1.png', width = 1000*1, height = 1000*0.8, res = 300)
plot_density(SS1, c('SOX2'), reduction = 'umap', size = 0.05) &
theme_light() &
theme(
# legend.position = 'None', 
panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.width= unit(1, 'mm')) &
scale_color_gradient2(mid = 'gray90', high = '#012a4a', midpoint = 0)
dev.off()
 
# 'KRT4', 'HAND1', 
png('P3F-2.png', width = 1000*3, height = 1000*3*0.8, res = 300)
plot_density(SS1, c('KRT19', 'KRT8','KRT5', 'WNT6', 'TP63', 'KRT1', 'GABRP', 'HAND1', 'IGFBP3'), reduction = 'umap', size = 0.05) &
theme_light() &
theme(
# legend.position = 'None', 
panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.width= unit(1, 'mm')) &
scale_color_gradient2(mid = 'gray90', high = '#012a4a', midpoint = 0)
dev.off()


# P3IA <- plot_density(P147D100[,P147D100$orig.ident %in% c('D6', 'D8')], c('IRX1', 'SIX1', 'SOX11'), reduction = 'umap', size = 0.05) &
# theme_light() &
# theme(panel.grid = element_blank(), 
# panel.border = element_blank(),
# plot.title = element_text(face = 3),
# legend.key.height= unit(1, 'mm'), 
# legend.position = 'None') &
# scale_color_gradient2(mid = 'gray90', high = '#012a4a', midpoint = 0)

# P3IB <- plot_density(P147D100[,P147D100$orig.ident %in% c('D6', 'D8')], c('PAX8', 'FGF8', 'PAX2'), reduction = 'umap', size = 0.05) &
# theme_light() &
# theme(panel.grid = element_blank(), 
# panel.border = element_blank(),
# plot.title = element_text(face = 3),
# legend.key.height= unit(1, 'mm'), 
# legend.position = 'None') &
# scale_color_gradient2(mid = 'gray90', high = '#90BE6D', midpoint = 0)

# P3IC <- plot_density(P147D100[,P147D100$orig.ident %in% c('D6', 'D8')], c('OTX2', 'TLX1', 'PAX3'), reduction = 'umap', size = 0.05) &
# theme_light() &
# theme(panel.grid = element_blank(), 
# panel.border = element_blank(),
# plot.title = element_text(face = 3),
# legend.key.height= unit(1, 'mm'), 
# legend.position = 'None') &
# scale_color_gradient2(mid = 'gray90', high = '#B48007', midpoint = 0)


P3IA <- plot_density(SS1, c('IRX1', 'SIX1', 'SOX11'), reduction = 'umap', size = 0.05) &
theme_light() &
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.height= unit(1, 'mm'), 
legend.position = 'None') &
scale_color_gradient2(mid = 'gray90', high = '#012a4a', midpoint = 0)

P3IB <- plot_density(SS1, c('PAX8', 'FGF8', 'PAX2'), reduction = 'umap', size = 0.05) &
theme_light() &
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.height= unit(1, 'mm'), 
legend.position = 'None') &
scale_color_gradient2(mid = 'gray90', high = '#90BE6D', midpoint = 0)

P3IC <- plot_density(SS1, c('OTX2', 'TLX1', 'PAX3'), reduction = 'umap', size = 0.05) &
theme_light() &
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.height= unit(1, 'mm'), 
legend.position = 'None') &
scale_color_gradient2(mid = 'gray90', high = '#B48007', midpoint = 0)

png('P3I.png', width = 1000*3*0.8, height = 1000*3*0.8, res = 300)
P3IA/P3IB/P3IC
dev.off()

P3IA1 <- plot_density(SS1, c('IRX1', 'SIX1', 'SOX11'), reduction = 'umap', size = 0.05, joint = TRUE) &
theme_light() &
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.height= unit(1, 'mm'), 
legend.position = 'None') &
scale_color_gradient2(mid = 'gray90', high = '#012a4a', midpoint = 0)

P3IB1 <- plot_density(SS1, c('PAX8', 'FGF8', 'PAX2'), reduction = 'umap', size = 0.05, joint = TRUE) &
theme_light() &
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.height= unit(1, 'mm'), 
legend.position = 'None') &
scale_color_gradient2(mid = 'gray90', high = '#90BE6D', midpoint = 0)

P3IC1 <- plot_density(SS1, c('OTX2', 'TLX1', 'PAX3'), reduction = 'umap', size = 0.05, joint = TRUE) &
theme_light() &
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.height= unit(1, 'mm'), 
legend.position = 'None') &
scale_color_gradient2(mid = 'gray90', high = '#B48007', midpoint = 0)

# P3IA1 <- plot_density(SS1, c('IRX1', 'SIX1', 'SOX11'), reduction = 'umap', size = 0.05, joint = TRUE) &
# theme_light() &
# theme(panel.grid = element_blank(), 
# panel.border = element_blank(),
# plot.title = element_text(face = 3),
# legend.key.height= unit(1, 'mm'), 
# legend.position = 'None') &
# scale_color_gradient2(mid = 'gray90', high = '#012a4a', midpoint = 0)

# P3IB1 <- plot_density(SS1, c('PAX8', 'FGF8', 'PAX2'), reduction = 'umap', size = 0.05, joint = TRUE) &
# theme_light() &
# theme(panel.grid = element_blank(), 
# panel.border = element_blank(),
# plot.title = element_text(face = 3),
# legend.key.height= unit(1, 'mm'), 
# legend.position = 'None') &
# scale_color_gradient2(mid = 'gray90', high = '#90BE6D', midpoint = 0)

# P3IC1 <- plot_density(SS1, c('OTX2', 'TLX1', 'PAX3'), reduction = 'umap', size = 0.05, joint = TRUE) &
# theme_light() &
# theme(panel.grid = element_blank(), 
# panel.border = element_blank(),
# plot.title = element_text(face = 3),
# legend.key.height= unit(1, 'mm'), 
# legend.position = 'None') &
# scale_color_gradient2(mid = 'gray90', high = '#B48007', midpoint = 0)

png('P3IALL.png', height = 1000*3*0.8, width = 1000*0.8, res = 300)
P3IA1[[4]] / P3IB1[[4]] / P3IC1[[4]]
dev.off()

load('SS1.RData')
load('scVelocity_SS1.RData')
library(dplyr)
umap_arrows <- arrows$garrows %>%
    as.data.frame() %>%
    mutate(x2 = x0 + (x1 - x0) * 2,
           y2 = y0 + (y1 - y0) * 2)
png('umapSS1.png', width = 1500, height = 1500, res = 300)
UMAPPlot(SS1) +
theme_light() +
theme(legend.position = 'bottom', 
        panel.grid = element_blank(),
        panel.border = element_blank()) +
        geom_curve(data = umap_arrows, curvature = 0.2, angle = 45,
                 aes(x = x0, xend = x2, y = y0, yend = y2),
                 size = .2,
                 arrow = arrow(length = unit(1, "points"), type = "closed"),
                 colour = "grey20", alpha = 0.8) +
                 xlab('UMAP 1') +
xlab('UMAP 1') +
ylab('UMAP 2') +
scale_color_manual(values = matchColors2[levels(Idents(SS1))])
dev.off()

load('SS1.RData')
load('scVelo_SS1.RData')

png('umapSS1-2.png', width = 1500, height = 1500, res = 300)
UMAPPlot(SS1) +
theme_light() +
theme(legend.position = 'bottom', 
        panel.grid = element_blank(),
        panel.border = element_blank()) +
        geom_segment(data=VF, 
        mapping=aes(
            x=start.1, 
            y=start.2, 
            xend=end.1, 
            yend=end.2), 
        size = 0.2,
        lineend = 'round', 
        linejoin = 'round',
        arrow=arrow(length=unit(0.5, "mm"))) +
        xlab('UMAP 1') + 
        ylab('UMAP 2') + 
scale_color_manual(values = matchColors2[levels(Idents(SS1))])
dev.off()

load('SS2.RData')
load('scVelocity_SS2.RData')
library(dplyr)
umap_arrows <- arrows$garrows %>%
    as.data.frame() %>%
    mutate(x2 = x0 + (x1 - x0) * 2,
           y2 = y0 + (y1 - y0) * 2)

png('umapSS2.png', width = 1500, height = 1500, res = 300)
UMAPPlot(SS2) +
theme_light() +
theme(legend.position = 'bottom', 
        panel.grid = element_blank(),
        panel.border = element_blank()) +
        geom_curve(data = umap_arrows, curvature = 0.2, angle = 45,
                 aes(x = x0, xend = x2, y = y0, yend = y2),
                 size = .2,
                 arrow = arrow(length = unit(1, "points"), type = "closed"),
                 colour = "grey20", alpha = 0.8) +
                 xlab('UMAP 1') +
xlab('UMAP 1') +
ylab('UMAP 2') +
scale_color_manual(values = matchColors2[levels(Idents(SS2))])
dev.off()

load('SS2.RData')
load('scVelo_SS2.RData')
png('umapSS2-2.png', width = 1500, height = 1500, res = 300)
UMAPPlot(SS2) +
theme_light() +
theme(legend.position = 'bottom', 
        panel.grid = element_blank(),
        panel.border = element_blank()) +
        geom_segment(data=VF, 
        mapping=aes(
            x=start.1, 
            y=start.2, 
            xend=end.1, 
            yend=end.2), 
        size = 0.2,
        lineend = 'round', 
        linejoin = 'round',
        arrow=arrow(length=unit(0.5, "mm"))) +
        xlab('UMAP 1') + 
        ylab('UMAP 2') + 
scale_color_manual(values = matchColors2[levels(Idents(SS2))])
dev.off()

png('P3O.png', width = 800*3, height = 800*2, res = 300)
PA <- plot_density(SS2, c('POU3F2', 'GDF7', 'PAX6'), reduction = 'umap') &
theme_light() &
theme(legend.position = 'None', 
panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.width= unit(1, 'mm')) &
scale_color_gradient2(mid = 'gray90', high = '#9A3C09', midpoint = 0)
PB <- plot_density(SS2, c('SOX10', 'FOXD3', 'SNAI2'), reduction = 'umap') &
theme_light() &
theme(legend.position = 'None', 
panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.width= unit(1, 'mm')) &
scale_color_gradient2(mid = 'gray90', high = '#F6AF0A', midpoint = 0)
PA/PB 
dev.off()


png('P3SOX2.png', width = 800, height = 800, res = 300)
PA <- plot_density(SS2, c('SOX2'), reduction = 'umap') &
theme_light() &
theme(legend.position = 'None', 
panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.width= unit(1, 'mm')) &
scale_color_gradient2(mid = 'gray90', high = '#012a4a', midpoint = 0)
PA
dev.off()

## Cell Communication Analyses
library('CellChat')
set.seed(1)
D6D8 <- subset(D6D8, idents = levels(CP$CT))

CCD6D8 <- doCellComHuman(D6D8)
resultsCCD6D8 <- netVisual_bubble(CCD6D8, return.data = TRUE)[[1]]
write.csv(resultsCCD6D8, 'CCD6D8.csv', row.names = FALSE)

labelsKeep <- table(CCD6D8@idents)
mat <- CCD6D8@net$weight
selectedLabels <- rownames(mat)[rowSums(mat) > 0]
groupSize <- as.numeric(labelsKeep[selectedLabels])
mat <- mat[selectedLabels, selectedLabels]

png('CCALL.png', width = 500*5, height = 500*3, res = 300)
par(mfrow = c(3,5), mar = c(.1,1,1,.1), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, 
  color.use = matchColors2[rownames(mat)],
  vertex.label.color = matchColors2[rownames(mat)],
  vertex.label.cex = 0.7,
  vertex.weight = groupSize, 
  weight.scale = T, 
  edge.weight.max = 7, 
  alpha.edge = 0.7,
  title.name = rownames(mat)[i])
}
dev.off()

png('CCHM.png', width = 500*5, height = 500*5, res = 300)
netVisual_heatmap(CCD6D8, color.use = matchColors2[levels(CCD6D8@idents)])
dev.off()

sapply(CCD6D8@netP$pathways, function(pathway.show){
    png(paste0('D6D8_', pathway.show,'.png'), width = 1500, height = 1500, res = 300)
    par(bg = NA, mar = c(0,0,0,0), oma = c(0,0,0,0))
    try(netVisual_aggregate(CCD6D8, signaling = pathway.show, layout = "chord", color.use = matchColors2[levels(CCD6D8@idents)]), silent = TRUE)
    dev.off()
})


# Figure 4 Analyses
# D10C <- Read10X('/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler/data/D10C')
# D10N <- Read10X('/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler/data/D10N')
# D13C <- Read10X('/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler/data/D13C')
# D13N <- Read10X('/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler/data/D13N')

# D10C <- CreateSeuratObject(D10C, project = 'D10C')
# D10N <- CreateSeuratObject(D10N, project = 'D10N')
# D13C <- CreateSeuratObject(D13C, project = 'D13C')
# D13N <- CreateSeuratObject(D13N, project = 'D13N')

# CHIR <- merge(D10C, c(D10N, D13C, D13N))
# CHIR <- scQC(CHIR)
# CHIR <- NormalizeData(CHIR)
# CHIR <- FindVariableFeatures(CHIR)
# CHIR <- ScaleData(CHIR)
# CHIR <- RunPCA(CHIR)
# CHIR <- RunUMAP(CHIR, dims = 1:30)

# save(CHIR, file = 'CHIR.RData')
load('CHIR_Annotated.RData')
#CHIR <- CHIR[,CHIR$labelProbability >= 0.8]


OEPDOEP <- CHIR[,Idents(CHIR) %in% c('OEPD', 'OEP')]
OEPDOEP$Treatment <- ifelse(grepl('C', OEPDOEP$orig.ident), 'CHIR+', 'CHIR-')
OEPDOEP$Treatment <- factor(OEPDOEP$Treatment, c('CHIR-', 'CHIR+'))
Idents(OEPDOEP) <- OEPDOEP$Treatment

sapply(c('DKK1', 'LEF1', 'TCF4', 'SP5', 'CCND1', 'NKD1', 'AXIN2', 'JUN', 'AKAP12', 'DLX5', 'APOE', 'ALDH1A1'), function(id){
        png(paste0('chirVln/', id, '.png'), width = 600, height = 1000, res = 300)
        df <- data.frame(E = OEPDOEP@assays$RNA@data[id,], G = OEPDOEP$Treatment)
        P <- ggplot(df, aes(G,E, fill = G)) +
        geom_violin(size = 0.1) + 
        geom_jitter(width = 0.25, size = 0.05) +
        labs(title = id) +
        theme_minimal() +
        theme(panel.grid = element_blank(),
        legend.position = 'None',
        plot.title = element_text(face = 3)) +
        xlab(NULL) +
        ylab('Expression Level') + 
        scale_fill_brewer(palette="Paired")
        print(P)
        dev.off()
})

DE <- FindMarkers(OEPDOEP, 'CHIR+', logfc.threshold = 0)
DE$label <- rownames(DE)
DE$label[DE$p_val_adj >= 0.05] <- NA
DE$label[abs(DE$avg_log2FC) < 1] <- NA
DE$color = 'black'
DE$color[!is.na(DE$label) & DE$avg_log2FC > 1] <- '#045a8d'
DE$color[!is.na(DE$label) & DE$avg_log2FC < -1] <- '#74a9cf'
DE$shape <- ifelse(is.na(DE$label), 16, 8)

png('DE_CHIR.png', width = 1900, height = 1900, res = 300)
ggplot(DE, aes(avg_log2FC, -log10(p_val), label = label)) +
geom_point(color = DE$color, shape = DE$shape) +
geom_text_repel(
        fontface = 3, 
        min.segment.length = 0, 
        segment.size = 0.2,
        bg.color = 'white') +
theme_minimal() +
theme(panel.grid = element_blank()) +
xlab(parse(text = 'log[2]~(Fold-Change)')) +
ylab(parse(text = '-log[10]~(P-value)'))
dev.off()

# library(UCell)
# library(fgsea)
# KEGG <- gmtPathways('/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler/data/KEGG_2021_Human.txt')
# CHIRPathways <- ScoreSignatures_UCell(CHIR@assays$RNA@data, KEGG, name = '', ncores = 20)
# save(CHIRPathways, file = 'UCell_CHIR.RData')
load('UCell_CHIR.RData')
sapply(colnames(CHIRPathways), function(id){
df <- data.frame(E = CHIRPathways[,id], I = CHIR$orig.ident, CT = CHIR$transferedCellType)
df$I <- factor(df$I, c('D10N', 'D13N', 'D10C', 'D13C'))
levels(df$I) <- c('D10-', 'D13-', 'D10+', 'D13+')
png(paste0('PathwayEnrichment/pathwayVln_', gsub('\\/', '', id), '.png'), width = 3500, height = 400, res = 300)
P <- ggplot(df, aes(I, E, fill = CT, color = CT)) +
theme_light() + 
geom_violin() +
facet_grid(cols = vars(CT)) +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
theme(axis.title.x = element_blank()) +
theme(panel.grid = element_blank(), panel.border = element_blank(), legend.position = 'None') +
scale_color_manual(values = matchColors2[levels(df$CT)]) +
scale_fill_manual(values = matchColors2[levels(df$CT)]) +
ylab(id) +
theme(strip.text = element_text(colour = 'black'))
plot(P)
dev.off()
})

load('UCell_CHIR.RData')
sapply(colnames(CHIRPathways), function(id){
df <- data.frame(E = CHIRPathways[,id], I = CHIR$orig.ident, CT = CHIR$transferedCellType)
df$I <- factor(df$I, c('D10N', 'D13N', 'D10C', 'D13C'))
levels(df$I) <- c('-', '-', '+', '+')
png(paste0('PathwayEnrichment/pathwayVln2_', gsub('\\/', '', id), '.png'), width = 3500, height = 400, res = 300)
P <- ggplot(df, aes(I, E, fill = CT, color = CT)) +
theme_light() + 
geom_violin() +
facet_grid(cols = vars(CT)) +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
theme(axis.title.x = element_blank()) +
theme(panel.grid = element_blank(), panel.border = element_blank(), legend.position = 'None') +
scale_color_manual(values = matchColors2[levels(df$CT)]) +
scale_fill_manual(values = matchColors2[levels(df$CT)]) +
ylab(id) +
theme(strip.text = element_text(colour = 'black'))
plot(P)
dev.off()
})

# library(UCell)
# library(fgsea)
# Hallmark <- gmtPathways('/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler/data/MSigDB_Hallmark_2020.txt')
# CHIRPathways <- ScoreSignatures_UCell(CHIR@assays$RNA@data, Hallmark, name = '', ncores = 20)

# K <- lapply(colnames(CHIRPathways), function(id){
# df <- data.frame(E = CHIRPathways[,id], I = CHIR$orig.ident, CT = CHIR$transferedCellType, P = id)
# df$I <- factor(df$I, c('D10N', 'D13N', 'D10C', 'D13C'))
# levels(df$I) <- c('D10-', 'D13-', 'D10+', 'D13+')
# return(df)
# })
# K <- do.call(rbind.data.frame, K)
# K$I <- factor(K$I, c('D10-', 'D13-', 'D10+', 'D13+'))
# png('H.png', width = 3500, height = 5000, res = 300)
# ggplot(K, aes(I, E, color = CT, fill = CT)) +
# theme_light() + 
# geom_violin() +
# facet_grid(cols = vars(CT), rows = vars(P), scales = 'free_y') +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   theme(axis.text.y = element_blank()) +
# theme(strip.text.y = element_text(angle = 0, vjust = 0.5)) +
# scale_color_manual(values = matchColors2[levels(K$CT)]) +
# scale_fill_manual(values = matchColors2[levels(K$CT)]) +
# theme(panel.grid = element_blank(), panel.border = element_blank(), legend.position = 'None') +
# theme(strip.text = element_text(colour = 'black')) +
# ylab('Enrichment Score') +
# xlab('Experiment')
# dev.off()

PA <- plot_density(CHIR[,CHIR$orig.ident %in% 'D10C'], c('DKK1'), 
reduction = 'umap', size = 0.05) + labs(title = 'Day 10 CHIR+')
PB <- plot_density(CHIR[,CHIR$orig.ident %in% 'D13C'], c('DKK1'), 
reduction = 'umap', size = 0.05) + labs(title = 'Day 13 CHIR+')
PC <- plot_density(CHIR[,CHIR$orig.ident %in% 'D10N'], c('DKK1'), 
reduction = 'umap', size = 0.05) + labs(title = 'Day 10 CHIR-')
PD <- plot_density(CHIR[,CHIR$orig.ident %in% 'D13N'], c('DKK1'), 
reduction = 'umap', size = 0.05) + labs(title = 'Day 13 CHIR-')
FP <- PA + PB + PC + PD
FP <- FP & theme_void() &
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
# plot.title = element_text(face = 3),
legend.key.width= unit(1, 'mm')) &
scale_color_gradient2(mid = 'gray90', high = '#012a4a', midpoint = 0, limits = c(-0.05,0.72))
FP <- FP + plot_annotation(title = 'DKK1', theme = theme(plot.title = element_text(face = 3))) 
png('CHIR-DKK1.png', width = 750*2, height = 750*2, res = 300)
print(FP)
dev.off()

PA <- plot_density(CHIR[,CHIR$orig.ident %in% 'D10C'], c('ALDH1A1'), 
reduction = 'umap', size = 0.05) + labs(title = 'Day 10 CHIR+')
PB <- plot_density(CHIR[,CHIR$orig.ident %in% 'D13C'], c('ALDH1A1'), 
reduction = 'umap', size = 0.05) + labs(title = 'Day 13 CHIR+')
PC <- plot_density(CHIR[,CHIR$orig.ident %in% 'D10N'], c('ALDH1A1'), 
reduction = 'umap', size = 0.05) + labs(title = 'Day 10 CHIR-')
PD <- plot_density(CHIR[,CHIR$orig.ident %in% 'D13N'], c('ALDH1A1'), 
reduction = 'umap', size = 0.05) + labs(title = 'Day 13 CHIR-')
FP <- PA + PB + PC + PD
FP <- FP & theme_void() &
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
# plot.title = element_text(face = 3),
legend.key.width= unit(1, 'mm')) &
scale_color_gradient2(mid = 'gray90', high = '#012a4a', midpoint = 0, limits = c(-0.05,0.72))
FP <- FP + plot_annotation(title = 'ALDH1A1', theme = theme(plot.title = element_text(face = 3))) 
png('CHIR-ALDH1A1.png', width = 750*2, height = 750*2, res = 300)
print(FP)
dev.off()

F1CMarkers <- c('SOX2', 'POU5F1', 'NANOG', 
'OTX2', 'RAX', 'GBX2', 'EN1', 'EPCAM',
'FOXD3', 'SNAI1', 'ETS1',
'MEF2C', 'SOX10', 'ERBB3', 'CAMK2B', 
'DLX2', 'FOXC2', 'PITX2', 'FLT1', 'SOX7', 
'TBX2', 'ACTA2', 'PSAP', 'GPX3')
png('t.png', width = 2000, height = 1100, res = 300)
#DotPlot(CHIR[,Idents(CHIR) %in% 'PSC2'], features = F1CMarkers) +
DotPlot(CHIR, features = F1CMarkers, col.min = 0, dot.min = 0.01, scale.by = 'size') +
theme_minimal() +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 3)) +
theme(panel.grid = element_blank()) +
scale_color_gradient2(mid = '#f4f5f6', high = '#012a4a', midpoint = 0) +
xlab('Marker Genes') +
ylab('Cell Types')
dev.off()

BG <- as.data.frame(P147D100@reductions$umap@cell.embeddings)
BG$TP <- 'BG'
BG$CT <- 'BG'
BG$COLOR <- 'gray90'
UG <- as.data.frame(CHIR@reductions$umap@cell.embeddings)
UG$TP <- as.vector(CHIR$orig.ident)
UG$CT <- as.vector(Idents(CHIR))
UG$COLOR <- matchColors2[UG$CT]
ALL <- rbind(BG, UG)

PA <- ggplot(ALL[ALL$TP %in% c('BG', 'D10N'),], aes(UMAP_1, UMAP_2)) +
geom_point(size = 0.01, color = ALL[ALL$TP %in% c('BG', 'D10N'),]$COLOR) +
theme_light() +
labs(title = 'Day 10 CHIR-')

PB <- ggplot(ALL[ALL$TP %in% c('BG', 'D10C'),], aes(UMAP_1, UMAP_2)) +
geom_point(size = 0.01, color = ALL[ALL$TP %in% c('BG', 'D10C'),]$COLOR) +
theme_light() +
labs(title = 'Day 10 CHIR+')

PC <- ggplot(ALL[ALL$TP %in% c('BG', 'D13N'),], aes(UMAP_1, UMAP_2)) +
geom_point(size = 0.01, color = ALL[ALL$TP %in% c('BG', 'D13N'),]$COLOR) +
theme_light()+
labs(title = 'Day 13 CHIR-')

PD <- ggplot(ALL[ALL$TP %in% c('BG', 'D13C'),], aes(UMAP_1, UMAP_2)) +
geom_point(size = 0.01, color = ALL[ALL$TP %in% c('BG', 'D13C'),]$COLOR) +
theme_light() +
labs(title = 'Day 13 CHIR+')

png('P147_UMAPCHIR.png', width = 2000, height = 2000, res = 300)
((PA | PB) / (PC | PD)) & 
theme(panel.border = element_blank(), panel.grid = element_blank()) &
xlab('UMAP 1') &
ylab('UMAP 2')
dev.off()

CP <- table(CHIR$orig.ident, CHIR$transferedCellType)
CP <- round((CP/rowSums(CP))*100,2)
apply(CP,2,function(Y){t.test(Y~c(1,0,1,0))$p.value})
CP <- reshape2::melt(CP)
colnames(CP) <- c('Sample', 'CT', 'Proportion')
CP$Sample <- factor(CP$Sample, levels = c('D13C', 'D13N', 'D10C', 'D10N'))
levels(CP$Sample) <- c('D13 CHIR+', 'D13 CHIR-', 'D10 CHIR+', 'D10 CHIR-')
png('P147_CTT.png', width = 2000, height = 800, res = 300)
ggplot(CP, aes(Proportion, Sample, fill = CT)) +
geom_bar(stat = 'identity') +
theme_light() +
theme(legend.position = 'bottom', 
        panel.grid = element_blank(),
        panel.border = element_blank()) +
        scale_fill_manual(values = matchColors2[levels(CP$CT)]) +
        guides(fill=guide_legend(nrow=3)) +
        labs(fill = 'Cell Type')
dev.off()

CP$G <- ifelse(grepl('\\+', CP$Sample), 'CHIR+', 'CHIR-')

library(dplyr)
CP <- CP %>% group_by(G, CT) %>% summarize(A = round(mean(Proportion),1), L = A - sd(Proportion), U = A + sd(Proportion))
#CP <- CP[!CP$CT %in% c('MEL', 'MYO', 'HC', 'NEP'),]
png('P147_CTS.png', width = 2000, height = 700, res = 300)
ggplot(CP, aes(CT, A, fill = G)) +
theme_light() + 
geom_errorbar(aes(ymin=A*0.7, ymax=U), width=0.4,
                 position=position_dodge(.9)) + 
geom_bar(stat = 'identity', position=position_dodge()) +
# geom_text(aes(label=A), hjust=1, color="black", angle = 90,
        #     position = position_dodge(0.8), size=3)+
scale_fill_brewer(palette="Paired") +
xlab('Cell Type') +
ylab('%') + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
theme(panel.grid = element_blank(), panel.border = element_blank()) +
theme(legend.position = 'bottom', legend.key.height = unit(1,'mm')) +
labs(fill = 'Treatment')
dev.off()

CHIR$IDPerSample <- paste0(CHIR$orig.ident, ' - ', CHIR$transferedCellType)
Idents(CHIR) <- CHIR$IDPerSample

# Pseudobulk analyses
# CHIRPB <- pbsapply(levels(Idents(CHIR)), function(ID){
#         rowSums(CHIR@assays$RNA@counts[,CHIR$IDPerSample %in% ID, drop = FALSE])
# })
# CHIRPBG <- ifelse(grepl('D10N|D13N', colnames(CHIRPB)), 0, 1)

# Communication analyses
library(CellChat)
Idents(CHIR) <- CHIR$transferedCellType
ccD10N <- doCellComHuman(CHIR[,CHIR$orig.ident %in% 'D10N'])
ccD10C <- doCellComHuman(CHIR[,CHIR$orig.ident %in% 'D10C'])
ccD10N <- liftCellChat(ccD10N, levels(Idents(CHIR)))
ccD10C <- liftCellChat(ccD10C, levels(Idents(CHIR)))
ccD10 <- mergeCellChat(list(D10N = ccD10N, D10C = ccD10C), 
                    add.names = c('D10N', 'D10C'), cell.prefix = TRUE)

png('ccD10_1.png', width = 900, height = 1000, res = 300)
compareInteractions(ccD10) + 
compareInteractions(ccD10, measure = 'weight') &
theme_light() +
theme(legend.position = 'None', 
      panel.grid = element_blank(), 
      panel.border = element_blank())
dev.off()

png('ccD10_2.png', width = 1200*1.1, height = 1000*1.1, res = 300)
netVisual_heatmap(ccD10, color.use = matchColors2[levels(Idents(CHIR))])
dev.off()

sapply(sort(levels(Idents(CHIR))), function(cellType){
  png(paste0('/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler/results/CCD10/hr_fSignaling_',gsub('/','_',cellType),'.png'), width = 1500, height = 2000, res = 300)
    try(print(rankNet(ccD10, 
            mode = "comparison", 
            title = cellType,
            sources.use = cellType,  
            stacked = TRUE, 
            do.stat = TRUE, 
            do.flip = TRUE) + 
            theme_light() +
            theme(panel.grid = element_blank(), 
                  panel.border = element_blank()
                  )
    ))
  dev.off()
})

ccD13N <- doCellComHuman(CHIR[,CHIR$orig.ident %in% 'D13N'])
ccD13C <- doCellComHuman(CHIR[,CHIR$orig.ident %in% 'D13C'])
ccD13N <- liftCellChat(ccD13N, levels(Idents(CHIR)))
ccD13C <- liftCellChat(ccD13C, levels(Idents(CHIR)))
ccD13 <- mergeCellChat(list(D13N = ccD13N, D13C = ccD13C), 
                    add.names = c('D13N', 'D13C'), cell.prefix = TRUE)

png('ccD13_1.png', width = 900, height = 1000, res = 300)
compareInteractions(ccD13) + 
compareInteractions(ccD13, measure = 'weight') &
theme_light() +
theme(legend.position = 'None', 
      panel.grid = element_blank(), 
      panel.border = element_blank())
dev.off()

png('ccD13_2.png', width = 1200*1.1, height = 1000*1.1, res = 300)
netVisual_heatmap(ccD13, color.use = matchColors2[levels(Idents(CHIR))])
dev.off()

sapply(sort(levels(Idents(CHIR))), function(cellType){
  png(paste0('/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler/results/CCD13/hr_fSignaling_',gsub('/','_',cellType),'.png'), width = 1500, height = 2000, res = 300)
    try(print(rankNet(ccD13, 
            mode = "comparison", 
            title = cellType,
            sources.use = cellType,  
            stacked = TRUE, 
            do.stat = TRUE, 
            do.flip = TRUE) + 
            theme_light() +
            theme(panel.grid = element_blank(), 
                  panel.border = element_blank()
                  )
    ))
  dev.off()
})


ccN <- doCellComHuman(CHIR[,CHIR$orig.ident %in% c('D10N', 'D13N')])
ccC <- doCellComHuman(CHIR[,CHIR$orig.ident %in% c('D10C', 'D13C')])
ccN <- liftCellChat(ccN, levels(Idents(CHIR)))
ccC <- liftCellChat(ccC, levels(Idents(CHIR)))
CC <- mergeCellChat(list(N = ccN, C = ccC), 
                    add.names = c('CHIR-', 'CHIR+'), cell.prefix = TRUE)

png('CC_1.png', width = 900, height = 1000, res = 300)
compareInteractions(CC) + 
compareInteractions(CC, measure = 'weight') &
theme_light() +
theme(legend.position = 'None', 
      panel.grid = element_blank(), 
      panel.border = element_blank())
dev.off()

png('CC_2.png', width = 1200*1.1, height = 1000*1.1, res = 300)
netVisual_heatmap(CC, color.use = matchColors2[levels(Idents(CHIR))])
dev.off()

sapply(sort(levels(Idents(CHIR))), function(cellType){
  png(paste0('/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler/results/CC/hr_fSignaling_',gsub('/','_',cellType),'.png'), width = 1500, height = 2000, res = 300)
    try(print(rankNet(CC, 
            mode = "comparison", 
            title = cellType,
            sources.use = cellType,  
            stacked = TRUE, 
            do.stat = TRUE, 
            do.flip = TRUE) + 
            theme_light() +
            theme(panel.grid = element_blank(), 
                  panel.border = element_blank()
                  )
    ))
  dev.off()
})

write.csv(netVisual_bubble(ccD10N, return.data = TRUE)[[1]], 'results/allComm_D10N.csv')
write.csv(netVisual_bubble(ccD10C, return.data = TRUE)[[1]], 'results/allComm_D10C.csv')
write.csv(netVisual_bubble(ccD13N, return.data = TRUE)[[1]], 'results/allComm_D13N.csv')
write.csv(netVisual_bubble(ccD13C, return.data = TRUE)[[1]], 'results/allComm_D13C.csv')
write.csv(netVisual_bubble(ccN, return.data = TRUE)[[1]], 'results/allComm_CHIR-.csv')
write.csv(netVisual_bubble(ccC, return.data = TRUE)[[1]], 'results/allComm_CHIR+.csv')


sapply(ccN@netP$pathways, function(pathway.show){
    png(paste0('results/CC/ccN_', pathway.show,'.png'), width = 1500, height = 1500, res = 300)
    par(bg = NA, mar = c(0,0,0,0), oma = c(0,0,0,0))
    try(netVisual_aggregate(ccN, signaling = pathway.show, layout = "chord", color.use = matchColors2[levels(ccN@idents)]), silent = TRUE)
    dev.off()
})

sapply(ccC@netP$pathways, function(pathway.show){
    png(paste0('results/CC/ccC_', pathway.show,'.png'), width = 1500, height = 1500, res = 300)
    par(bg = NA, mar = c(0,0,0,0), oma = c(0,0,0,0))
    try(netVisual_aggregate(ccC, signaling = pathway.show, layout = "chord", color.use = matchColors2[levels(ccC@idents)]), silent = TRUE)
    dev.off()
})


#### Figure 6
F6Data <- P147D100[,Idents(P147D100) %in% c('HC', 'OEP')]

FC <- FoldChange(F6Data, 'OEP')
FC <- FC[order(FC[,1], decreasing = FALSE),]
png('F6-HC.png', width = 4000, height =1200, res = 300)
VlnPlot(F6Data, 
features = c('ATOH1', 'MYO7A', 'OTOF', 'STRC', 'ESPN',
'GPX2', 'PCDH15', 'CDH23', 'USH2A', 'POU4F3', 'CABP2',
'GFI1', 'USH1C', 'RIPOR2', 'MYO6', 'MYO15A', 'CIB2', 
'PCP4', 'CALM2', 'LHFPL5'), 
pt.size = 0, ncol = 10) &
theme_minimal() &
theme(panel.grid = element_blank(), legend.position = 'None') &
xlab(NULL) & 
scale_fill_manual(values = c('#26604F','#445C71')) &
theme(plot.title = element_text(face = 3))
dev.off()
png('F6-HC2.png', width = 4000, height =1200, res = 300)
VlnPlot(F6Data, 
features = rownames(FC)[1:20], 
pt.size = 0, ncol = 10) &
theme_minimal() &
theme(panel.grid = element_blank(), legend.position = 'None') &
xlab(NULL) & 
scale_fill_manual(values = c('#26604F','#445C71')) &
theme(plot.title = element_text(face = 3))
dev.off()

png('F6-OEP.png', width = 4000, height =1200, res = 300)
VlnPlot(F6Data, 
features = c('PAX8', 'PAX2', 'OC90', 'HMX3', 'DLX3', 'DLX5',
'SALL4', 'DUSP6', 'SPRY2', 'SIX1', 'EYA1', 'FOXG1', 'LMX1A',
'OTOA', 'APOE', 'SMOC2', 'SPARCL1', 'FBXO2', 'COL11A1', 'COL9A2'), 
pt.size = 0, ncol = 10) &
theme_minimal() &
theme(panel.grid = element_blank(), legend.position = 'None') &
xlab(NULL) & 
scale_fill_manual(values = c('#26604F','#445C71')) &
theme(plot.title = element_text(face = 3))
dev.off()
FC <- FC[order(FC[,1], decreasing = TRUE),]
png('F6-OEP2.png', width = 4000, height =1200, res = 300)
VlnPlot(F6Data, 
features = rownames(FC)[1:20], 
pt.size = 0, ncol = 10) &
theme_minimal() &
theme(panel.grid = element_blank(), legend.position = 'None') &
xlab(NULL) & 
scale_fill_manual(values = c('#26604F','#445C71')) &
theme(plot.title = element_text(face = 3))
dev.off()

png('F6-DP.png', width = 3000, height = 600, res = 300)
DotPlot(F6Data, features=c('PAX8', 'PAX2', 'OC90', 'HMX3', 'DLX3', 'DLX5',
'SALL4', 'DUSP6', 'SPRY2', 'SIX1', 'EYA1', 'FOXG1', 'LMX1A',
'OTOA', 'APOE', 'SMOC2', 'SPARCL1', 'FBXO2', 'COL11A1', 'COL9A2',
'ATOH1', 'MYO7A', 'OTOF', 'STRC', 'ESPN',
'GPX2', 'PCDH15', 'CDH23', 'USH2A', 'POU4F3', 'CABP2',
'GFI1', 'USH1C', 'RIPOR2', 'MYO6', 'MYO15A', 'CIB2', 
'PCP4', 'CALM2', 'LHFPL5'), dot.min = 0.01, col.min = 0, group.by = 'CT') +
theme_light() + 
theme(panel.grid = element_blank(), panel.border = element_blank()) + 
theme(axis.text.x = element_text(angle = 90, face = 3, vjust = 0.5, hjust=1)) +
theme(legend.position = 'bottom', legend.key.height = unit(1,'mm')) +
scale_color_gradient2(mid = 'gray90', high = '#012a4a', midpoint = 0) +
xlab(NULL) +
ylab(NULL)
dev.off()

CP <- reshape2::melt(table(F6Data$orig.ident, Idents(F6Data)))
CP <- CP[CP[,1] %in% c('D10', 'D13', 'D18', 'D24', 'D30', 'D36'),]
colnames(CP) <- c('Day', 'CT', 'Cells')
CP$Day <- factor(CP$Day, c('D10', 'D13', 'D18', 'D24', 'D30', 'D36'))
png('F6-T.png', width = 500, height = 700, res = 300)
ggplot(CP, aes(Day, Cells, fill = CT)) +
geom_bar(stat = 'identity') +
facet_grid(rows = vars(CT), scales = 'free_y') +
scale_fill_manual(values = c('#26604F','#445C71')) +
theme_minimal() +
theme(panel.grid = element_blank()) +
theme(legend.position = 'None') +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
ylab(NULL) +
theme(strip.text.y = element_text(angle = 0))
dev.off()

D <- P147D100@reductions$umap@cell.embeddings
D <- as.data.frame.array(D)
D$ID <- P147D100$ID
D$Color <- matchColors2[D$ID]
D$TP <- P147D100$orig.ident
D$Color[D$ID != c('HC', 'OEP')] <- 'gray95'
D <- D[order(D$Color, decreasing = TRUE),]
P6 <- ggplot(D, aes(UMAP_1, UMAP_2)) +
geom_point(color = D$Color, size = 0.01) +
theme_light() +
theme(legend.position = 'None', 
        panel.grid = element_blank(),
        panel.border = element_blank()) +
xlab('UMAP 1') +
ylab('UMAP 2') +
scale_color_manual(values = matchColors2[levels(Idents(P147D100[,P147D100$orig.ident == 'D3']))]) +
geom_text_repel(data = LabelPosition[LabelPosition$Label %in% c('HC', 'OEP'),], 
mapping = aes(X,Y, label = Label), 
bg.color = 'white', 
box.padding = 0,   
min.segment.length = 0) 

png('F6-UMAP.png', width = 1500, height = 1500, res = 300)
print(P6)
dev.off()

F6Data <- P147D100[,Idents(P147D100) %in% c('HC', 'OEP')]
F6Data <- F6Data[,F6Data$orig.ident %in% c('D18', 'D24', 'D30', 'D36')]
F6Data$orig.ident <- factor(F6Data$orig.ident, c('D18', 'D24', 'D30', 'D36'))
levels(F6Data$orig.ident) <- c('Day 18', 'Day 24', 'Day 30', 'Day 36')
# F6Data <- FindVariableFeatures(F6Data)
# F6Data <- ScaleData(F6Data)
# # set.seed(1)
# F6Data <- RunPCA(F6Data)
# # set.seed(1)
# F6Data <- RunHarmony(F6Data, 'orig.ident')
# # set.seed(1)
# F6Data <- RunUMAP(F6Data, reduction = 'harmony', dims = 1:30)

png('F6-ISL1-2.png', width = (2200/3), height = 650, res = 300)
plot_density(F6Data, c('ISL1'), reduction = 'umap', size = 0.25) &
theme_void() &
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.width= unit(1, 'mm'), 
legend.position = 'right') &
scale_color_gradient2(mid = 'gray90', high = matchColors2['OEP'], midpoint = 0) &
xlim(c(2,5.4)) &
ylim(c(-3.5,0)) 
dev.off()

png('F6-PAX8-3.png', width = (2200/3), height = 650, res = 300)
plot_density(F6Data, c('PAX8'), reduction = 'umap', size = 0.25) &
theme_void() &
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.width= unit(1, 'mm'), 
legend.position = 'right') &
scale_color_gradient2(mid = 'gray90', high = matchColors2['OEPD'], midpoint = 0) &
xlim(c(2,5.4)) &
ylim(c(-3.5,0)) 
dev.off()

png('F6-R.png', width = (2200/3)*2, height = 650, res = 300)
DimPlot(F6Data, split.by = 'orig.ident', ncol = 4, pt.size = 0.05) + 
theme_light() +
theme(legend.position = 'none', 
panel.grid = element_blank(), 
panel.border = element_blank()) +
scale_color_manual(values = matchColors2[levels(Idents(F6Data))]) +
xlab('UMAP 1') +
ylab('UMAP 2') +
theme(strip.text.x = element_text(color = 'black', hjust=0), 
 strip.background = element_rect(fill = rgb(1,1,1,0))) +
 xlim(c(2,5.4)) +
ylim(c(-3.5,0)) +
theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()

png('F6-M1.png', width = 2200, height = 650, res = 300)
plot_density(F6Data, 
c('SOX2', 'OTOG', 'USH1C'), 
reduction = 'umap', size = 0.25) &
xlab('UMAP 1') &
ylab('UMAP 2') &
theme_void() & 
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
legend.key.width= unit(0.5, 'mm'),
plot.title = element_text(face = 3)) &
scale_color_gradient2(mid = 'gray90', high = '#26604F', midpoint = 0) &
xlim(c(2,5.4)) &
ylim(c(-3.5,0)) 
dev.off()

png('F6-M2.png', width = 2200, height = 650, res = 300)
plot_density(F6Data, 
c('DLX5', 'LMX1A', 'WNT2B'), 
reduction = 'umap', size = 0.25) &
xlab('UMAP 1') &
ylab('UMAP 2') &
theme_void() & 
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
legend.key.width= unit(0.5, 'mm'),
plot.title = element_text(face = 3)) &
scale_color_gradient2(mid = 'gray90', high = '#3D997D', midpoint = 0) & 
xlim(c(2,5.4)) &
ylim(c(-3.5,0)) 
dev.off()

png('F6-M3.png', width = (2200/3)*2, height = 650*2, res = 300)
plot_density(F6Data, 
c('ATOH1', 'POU4F3', 'OCM', 'ANXA4'), 
reduction = 'umap', size = 0.25) &
xlab('UMAP 1') &
ylab('UMAP 2') &
theme_void() & 
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
legend.key.width= unit(0.5, 'mm'),
plot.title = element_text(face = 3)) &
scale_color_gradient2(mid = 'gray90', high = '#445C71', midpoint = 0) & 
xlim(c(2,5.4)) &
ylim(c(-3.5,0)) 
dev.off()


MT <- FindAllMarkers(P147D100, only.pos = TRUE)
write.csv(MT, file = "P147D100_CTMarkers.csv")

png('F6-ISL1-PAX8-1.png', width = 2000, height = 1000, res = 300)
plot_density(P147D100, c('ISL1', 'PAX8'), reduction = 'umap', size = 0.05) &
theme_void() &
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.width= unit(1, 'mm'), 
legend.position = 'right') &
scale_color_gradient2(mid = 'gray90', high = '#012a4a', midpoint = 0)
dev.off()

OEP <- P147D100[,Idents(P147D100) %in% 'OEP']

png('F6-ISL1-1.png', width = 800, height = 700, res = 300)
plot_density(OEP, c('ISL1'), reduction = 'umap', size = 0.25) &
theme_void() &
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.width= unit(1, 'mm'), 
legend.position = 'right') &
scale_color_gradient2(mid = 'gray90', high = matchColors2['OEP'], midpoint = 0) &
xlim(c(2,5)) &
ylim(c(-3,0))
dev.off()

png('F6-PAX8-2.png', width = 800, height = 700, res = 300)
plot_density(OEP, c('PAX8'), reduction = 'umap', size = 0.25) &
theme_void() &
theme(panel.grid = element_blank(), 
panel.border = element_blank(),
plot.title = element_text(face = 3),
legend.key.width= unit(1, 'mm'), 
legend.position = 'right') &
scale_color_gradient2(mid = 'gray90', high = matchColors2['OEPD'], midpoint = 0) &
xlim(c(2,5)) &
ylim(c(-3,0))
dev.off()

newIdents <- as.vector(Idents(OEP))
newIdents[OEP@assays$RNA@counts['ISL1',] != 0] <- 'ISL1'
newIdents[OEP@assays$RNA@counts['PAX8',] != 0] <- 'PAX8'
newIdents[colSums(OEP@assays$RNA@counts[c('ISL1', 'PAX8'),] != 0) == 2] <- 'OEP'


Idents(OEP) <- newIdents
DE <- FindMarkers(OEP, 'ISL1', 'PAX8', logfc.threshold = 0)
DE$label <- rownames(DE)
DE$label[DE$p_val_adj >= 0.05] <- NA
DE$label[abs(DE$avg_log2FC) < 1] <- NA
DE$color = 'black'
DE$color[!is.na(DE$label) & DE$avg_log2FC > 1] <- matchColors2['OEP']
DE$color[!is.na(DE$label) & DE$avg_log2FC < -1] <- matchColors2['OEPD']
DE$shape <- ifelse(is.na(DE$label), 16, 8)

png('DE_ISL1-PAX8.png', width = 1900, height = 1900, res = 300)
ggplot(DE, aes(avg_log2FC, -log10(p_val), label = label)) +
geom_point(color = DE$color, shape = DE$shape) +
geom_text_repel(
        fontface = 3, 
        min.segment.length = 0, 
        segment.size = 0.2,
        bg.color = 'white', max.overlaps = 10, max.iterations = 1e3) +
theme_minimal() +
theme(panel.grid = element_blank()) +
xlab(parse(text = 'log[2]~(Fold-Change)')) +
ylab(parse(text = '-log[10]~(P-value)'))
dev.off()