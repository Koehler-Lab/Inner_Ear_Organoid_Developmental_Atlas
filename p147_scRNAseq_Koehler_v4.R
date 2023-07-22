# Libraries
snapshot <- "/lab-share/RC-Data-Science-e2/Public/DST_software/pipeline_R_snapshot/r4.1-20220803-bch-dst-v1/"
.libPaths(snapshot)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(Nebulosa)
library(patchwork)
library(pbapply)
library(ggrepel)
library(harmony)

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

load('scVelocity_Results.RData')
library(dplyr)
umap_arrows <- arrows$garrows %>%
    as.data.frame() %>%
    mutate(x2 = x0 + (x1 - x0) * 4,
           y2 = y0 + (y1 - y0) * 4)

png('V.png', width = 2000, height = 2000, res = 300)
UMAPPlot(P147D100, label = FALSE) +
theme_light() +
theme(legend.position = 'None', 
        panel.grid = element_blank(),
        panel.border = element_blank()) +
        geom_curve(data = umap_arrows, curvature = 0.2, angle = 45,
                 aes(x = x0, xend = x2, y = y0, yend = y2),
                 size = 1,
                 arrow = arrow(length = unit(4, "points"), type = "closed"),
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
    mutate(x2 = x0 + (x1 - x0) * 4,
           y2 = y0 + (y1 - y0) * 4)
png('umapSS1.png', width = 1500, height = 1500, res = 300)
UMAPPlot(SS1) +
theme_light() +
theme(legend.position = 'bottom', 
        panel.grid = element_blank(),
        panel.border = element_blank()) +
        geom_curve(data = umap_arrows, curvature = 0.2, angle = 45,
                 aes(x = x0, xend = x2, y = y0, yend = y2),
                 size = 1,
                 arrow = arrow(length = unit(4, "points"), type = "closed"),
                 colour = "grey20", alpha = 0.8) +
                 xlab('UMAP 1') +
xlab('UMAP 1') +
ylab('UMAP 2') +
scale_color_manual(values = matchColors2[levels(Idents(SS1))])
dev.off()

load('SS2.RData')
load('scVelocity_SS2.RData')
library(dplyr)
umap_arrows <- arrows$garrows %>%
    as.data.frame() %>%
    mutate(x2 = x0 + (x1 - x0) * 4,
           y2 = y0 + (y1 - y0) * 4)

png('umapSS2.png', width = 1500, height = 1500, res = 300)
UMAPPlot(SS2) +
theme_light() +
theme(legend.position = 'bottom', 
        panel.grid = element_blank(),
        panel.border = element_blank()) +
        geom_curve(data = umap_arrows, curvature = 0.2, angle = 45,
                 aes(x = x0, xend = x2, y = y0, yend = y2),
                 size = 1,
                 arrow = arrow(length = unit(4, "points"), type = "closed"),
                 colour = "grey20", alpha = 0.8) +
                 xlab('UMAP 1') +
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


## Cell Communication Analyses
library('CellChat')
set.seed(1)
D6D8 <- subset(D6D8, idents = levels(CP$CT))

CCD6D8 <- doCellComHuman(D6D8)
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