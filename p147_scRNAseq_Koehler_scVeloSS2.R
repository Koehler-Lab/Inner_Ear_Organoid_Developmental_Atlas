# Libraries
snapshot <- "/lab-share/RC-Data-Science-e2/Public/DST_software/pipeline_R_snapshot/r4.1-20220803-bch-dst-v1/"
.libPaths(snapshot)
source('/lab-share/RC-Data-Science-e2/Public/Daniel/rFunctions.R')

library(Seurat)
library(dplyr)
library(ggplot2)
library(velocyto.R)
library(velociraptor)

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

# Loading SS1
load('/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler/SS2.RData')
X <- X[,colnames(SS2)]
SS2 <- SS2[,colnames(X)]

# Velocity
V <- doVelocity(X, 'stochastic')
VF <- getVectorField(SS2,V, reduction = 'umap', resolution = 60)

save(VF, file = 'scVelo_SS2.RData')

png('t.png')
UMAPPlot(SS2) + 
theme_light() + 
theme(legend.position = 'None') + 
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
        ylab('UMAP 2')
dev.off()