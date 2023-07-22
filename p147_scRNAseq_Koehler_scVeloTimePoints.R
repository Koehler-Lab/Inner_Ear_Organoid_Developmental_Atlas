# Libraries
snapshot <- "/lab-share/RC-Data-Science-e2/Public/DST_software/pipeline_R_snapshot/r4.1-20220803-bch-dst-v1/"
.libPaths(snapshot)
source('/lab-share/RC-Data-Science-e2/Public/Daniel/rFunctions.R')

library(Seurat)
library(dplyr)
library(ggplot2)
library(velocyto.R)
library(velociraptor)
library(ggrepel)

# Setting WD
setwd('/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler')

# Loading Seurat Object
load('results/D100/P147.RData')
X$orig.bc <- colnames(X)

# Splitting
D0 <- X[,X$orig.ident %in% 'D0']
D0 <- addSUmatrices(D0, '/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler/IEOdev_DSP_d0.loom')

D3 <- X[,X$orig.ident %in% 'D3']
D3 <- addSUmatrices(D3, '/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler/IEOdev_DSP_d3.loom')

D6 <- X[,X$orig.ident %in% 'D6']
D6 <- addSUmatrices(D6, '/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler/IEOdev_DSP_d6.loom')

D8 <- X[,X$orig.ident %in% 'D8']
D8 <- addSUmatrices(D8, '/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler/IEOdev_DSP_d8.loom')

D10 <- X[,X$orig.ident %in% 'D10']
D10 <- addSUmatrices(D10, '/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler/IEOdev_DSP_d10N.loom')

D13 <- X[,X$orig.ident %in% 'D13']
D13 <- addSUmatrices(D13, '/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler/IEOdev_DSP_d13N.loom')

D18 <- X[,X$orig.ident %in% 'D18']
D18 <- addSUmatrices(D18, '/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler/IEOdev_DSP_d18.loom')

D24 <- X[,X$orig.ident %in% 'D24']
D24 <- addSUmatrices(D24, '/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler/IEOdev_DSP_d24.loom')

D30 <- X[,X$orig.ident %in% 'D30']
D30 <- addSUmatrices(D30, '/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler/IEOdev_DSP_d30.loom')

D36 <- X[,X$orig.ident %in% 'D36']
D36 <- addSUmatrices(D36, '/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler/IEOdev_DSP_d36.loom')

# Re-joining
X <- merge(D0, c(D3, D6, D8, D10, D13, D18, D24, D30, D36), merge.dr = c("umap"))
rm(D0, D3, D6, D8, D10, D13, D18, D24, D30, D36)
gc()
X <- RenameCells(X, new.names = X$orig.bc)

# Velocity
tpData <- X[,X$orig.ident %in% c('D0', 'D3')]
tpV <- doVelocity(tpData, 'stochastic')
tpVF <- getVectorField(tpData,tpV, reduction = 'umap', resolution = 50)
save(tpVF, file = 'scVelo_D0D3.RData')

tpData <- X[,X$orig.ident %in% c('D3', 'D6', 'D8')]
tpV <- doVelocity(tpData, 'stochastic')
tpVF <- getVectorField(tpData,tpV, reduction = 'umap', resolution = 50)
save(tpVF, file = 'scVelo_D3D6D8.RData')

tpData <- X[,X$orig.ident %in% c('D8', 'D10', 'D13')]
tpV <- doVelocity(tpData, 'stochastic')
tpVF <- getVectorField(tpData,tpV, reduction = 'umap', resolution = 50)
save(tpVF, file = 'scVelo_D8D10D13.RData')

tpData <- X[,X$orig.ident %in% c('D13', 'D18', 'D24', 'D30', 'D36')]
tpV <- doVelocity(tpData, 'stochastic')
tpVF <- getVectorField(tpData,tpV, reduction = 'umap', resolution = 50)
save(tpVF, file = 'scVelo_D13D18D24D30D36.RData')

load('/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler/P147D100_Annotated.RData')

tpColors <- tenColors <- c('#a9d6e5', '#89c2d9', '#61a5c2', '#468faf', '#2c7da0', '#2a6f97', '#014f86', '#01497c', '#013a63', '#012a4a')
ctColors <- c('#AF0508', '#F1070B', '#F94144', 
'#9A3C09', '#D4520C', '#F3722C', 
'#9B5805', '#D67907', '#F8961E', 
'#B48007', '#F6AF0A', '#F9C74F', 
'#527735', '#70A349', '#90BE6D', 
'#26604F', '#35846C', '#43AA8B', 
'#314352', '#445C71', '#577590')

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

png('V-D0D3.png', width = 1800, height = 1800, res = 300)
load('scVelo_D0D3.RData')
tp <- c('D0', 'D3')
D <- P147D100@reductions$umap@cell.embeddings
D <- as.data.frame.array(D)
D$ID <- P147D100$ID
D$Color <- matchColors2[D$ID]
D$TP <- P147D100$orig.ident
D$Color[!P147D100$orig.ident %in% tp] <- 'gray95'
D <- D[order(D$Color, decreasing = TRUE),]
sIDs <- table(D$ID[!D$Color %in% 'gray95'])
sIDs <- names(sIDs[(sIDs/sum(sIDs)) > 0.01])
P <- ggplot(D, aes(UMAP_1, UMAP_2)) +
  geom_point(color = D$Color, size = 0.01) +
  theme_void() +
  theme(legend.position = 'None', 
          panel.grid = element_blank(),
          panel.border = element_blank()) +
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  scale_color_manual(values = matchColors2[levels(Idents(P147D100[,P147D100$orig.ident == 'D3']))]) +
        geom_segment(data=tpVF, 
        mapping=aes(
            x=start.1, 
            y=start.2, 
            xend=end.1, 
            yend=end.2), 
        size = 0.5,
        lineend = 'round', 
        linejoin = 'round',
        arrow=arrow(length=unit(0.5, "mm")))
print(P)
dev.off()


png('V-D3D6D8.png', width = 1800, height = 1800, res = 300)
load('scVelo_D3D6D8.RData')
tp <- c('D3', 'D6', 'D8')
D <- P147D100@reductions$umap@cell.embeddings
D <- as.data.frame.array(D)
D$ID <- P147D100$ID
D$Color <- matchColors2[D$ID]
D$TP <- P147D100$orig.ident
D$Color[!P147D100$orig.ident %in% tp] <- 'gray95'
D <- D[order(D$Color, decreasing = TRUE),]
sIDs <- table(D$ID[!D$Color %in% 'gray95'])
sIDs <- names(sIDs[(sIDs/sum(sIDs)) > 0.01])
P <- ggplot(D, aes(UMAP_1, UMAP_2)) +
  geom_point(color = D$Color, size = 0.01) +
  theme_void() +
  theme(legend.position = 'None', 
          panel.grid = element_blank(),
          panel.border = element_blank()) +
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  scale_color_manual(values = matchColors2[levels(Idents(P147D100[,P147D100$orig.ident == 'D3']))]) +
        geom_segment(data=tpVF, 
        mapping=aes(
            x=start.1, 
            y=start.2, 
            xend=end.1, 
            yend=end.2), 
        size = 0.5,
        lineend = 'round', 
        linejoin = 'round',
        arrow=arrow(length=unit(0.5, "mm")))
print(P)
dev.off()

png('V-D8D10D13.png', width = 1800, height = 1800, res = 300)
load('scVelo_D8D10D13.RData')
tp <- c('D8', 'D10', 'D13')
D <- P147D100@reductions$umap@cell.embeddings
D <- as.data.frame.array(D)
D$ID <- P147D100$ID
D$Color <- matchColors2[D$ID]
D$TP <- P147D100$orig.ident
D$Color[!P147D100$orig.ident %in% tp] <- 'gray95'
D <- D[order(D$Color, decreasing = TRUE),]
sIDs <- table(D$ID[!D$Color %in% 'gray95'])
sIDs <- names(sIDs[(sIDs/sum(sIDs)) > 0.01])
P <- ggplot(D, aes(UMAP_1, UMAP_2)) +
  geom_point(color = D$Color, size = 0.01) +
  theme_void() +
  theme(legend.position = 'None', 
          panel.grid = element_blank(),
          panel.border = element_blank()) +
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  scale_color_manual(values = matchColors2[levels(Idents(P147D100[,P147D100$orig.ident == 'D3']))]) +
        geom_segment(data=tpVF, 
        mapping=aes(
            x=start.1, 
            y=start.2, 
            xend=end.1, 
            yend=end.2), 
        size = 0.5,
        lineend = 'round', 
        linejoin = 'round',
        arrow=arrow(length=unit(0.5, "mm")))
print(P)
dev.off()

png('V-D13D18D24D30D36.png', width = 1800, height = 1800, res = 300)
load('scVelo_D13D18D24D30D36.RData')
tp <- c('D13', 'D18', 'D24', 'D30', 'D36')
D <- P147D100@reductions$umap@cell.embeddings
D <- as.data.frame.array(D)
D$ID <- P147D100$ID
D$Color <- matchColors2[D$ID]
D$TP <- P147D100$orig.ident
D$Color[!P147D100$orig.ident %in% tp] <- 'gray95'
D <- D[order(D$Color, decreasing = TRUE),]
sIDs <- table(D$ID[!D$Color %in% 'gray95'])
sIDs <- names(sIDs[(sIDs/sum(sIDs)) > 0.01])
P <- ggplot(D, aes(UMAP_1, UMAP_2)) +
  geom_point(color = D$Color, size = 0.01) +
  theme_void() +
  theme(legend.position = 'None', 
          panel.grid = element_blank(),
          panel.border = element_blank()) +
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  scale_color_manual(values = matchColors2[levels(Idents(P147D100[,P147D100$orig.ident == 'D3']))]) +
        geom_segment(data=tpVF, 
        mapping=aes(
            x=start.1, 
            y=start.2, 
            xend=end.1, 
            yend=end.2), 
        size = 0.5,
        lineend = 'round', 
        linejoin = 'round',
        arrow=arrow(length=unit(0.5, "mm")))
print(P)
dev.off()





