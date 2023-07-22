# Libraries
snapshot <- "/lab-share/RC-Data-Science-e2/Public/DST_software/pipeline_R_snapshot/r4.1-20220803-bch-dst-v1/"
.libPaths(snapshot)
library(velocyto.R)
library(Seurat)
library(dplyr)
library(ggplot2)

# Setting WD
setwd('/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler')

# Loading Seurat Object
load('SS3.RData')
X <- F5Data
E0 <- X@reductions$umap@cell.embeddings
mdX <- X@meta.data
mdX <- as.data.frame.array(mdX)
mdX$BC <- rownames(mdX)
rownames(mdX) <- paste0(mdX$orig.ident, '_', gsub('-[[:print:]]+$','',rownames(mdX)))

# Reading Matrices
D0 <- read.loom.matrices("IEOdev_DSP_d0.loom")
colnames(D0[[1]]) <- paste0('D0_', gsub('IEOdev_DSP_d0:|x', '', colnames(D0[[1]])))
colnames(D0[[2]]) <- paste0('D0_', gsub('IEOdev_DSP_d0:|x', '', colnames(D0[[2]])))
D0[[1]] <- D0[[1]][VariableFeatures(X),]
D0[[2]] <- D0[[2]][VariableFeatures(X),]

D3 <- read.loom.matrices("IEOdev_DSP_d3.loom")
colnames(D3[[1]]) <- paste0('D3_', gsub('IEOdev_DSP_d3:|x', '', colnames(D3[[1]])))
colnames(D3[[2]]) <- paste0('D3_', gsub('IEOdev_DSP_d3:|x', '', colnames(D3[[2]])))
D3[[1]] <- D3[[1]][VariableFeatures(X),]
D3[[2]] <- D3[[2]][VariableFeatures(X),]

D6 <- read.loom.matrices("IEOdev_DSP_d6.loom")
colnames(D6[[1]]) <- paste0('D6_', gsub('IEOdev_DSP_d6:|x', '', colnames(D6[[1]])))
colnames(D6[[2]]) <- paste0('D6_', gsub('IEOdev_DSP_d6:|x', '', colnames(D6[[2]])))
D6[[1]] <- D6[[1]][VariableFeatures(X),]
D6[[2]] <- D6[[2]][VariableFeatures(X),]

D8 <- read.loom.matrices("IEOdev_DSP_d8.loom")
colnames(D8[[1]]) <- paste0('D8_', gsub('IEOdev_DSP_d8:|x', '', colnames(D8[[1]])))
colnames(D8[[2]]) <- paste0('D8_', gsub('IEOdev_DSP_d8:|x', '', colnames(D8[[2]])))
D8[[1]] <- D8[[1]][VariableFeatures(X),]
D8[[2]] <- D8[[2]][VariableFeatures(X),]

D10 <- read.loom.matrices("IEOdev_DSP_d10N.loom")
colnames(D10[[1]]) <- paste0('D10_', gsub('IEOdev_DSP_d10N:|x', '', colnames(D10[[1]])))
colnames(D10[[2]]) <- paste0('D10_', gsub('IEOdev_DSP_d10N:|x', '', colnames(D10[[2]])))
D10[[1]] <- D10[[1]][VariableFeatures(X),]
D10[[2]] <- D10[[2]][VariableFeatures(X),]

D13 <- read.loom.matrices("IEOdev_DSP_d13N.loom")
colnames(D13[[1]]) <- paste0('D13_', gsub('IEOdev_DSP_d13N:|x', '', colnames(D13[[1]])))
colnames(D13[[2]]) <- paste0('D13_', gsub('IEOdev_DSP_d13N:|x', '', colnames(D13[[2]])))
D13[[1]] <- D13[[1]][VariableFeatures(X),]
D13[[2]] <- D13[[2]][VariableFeatures(X),]

D18 <- read.loom.matrices("IEOdev_DSP_d18.loom")
colnames(D18[[1]]) <- paste0('D18_', gsub('IEOdev_DSP_d18:|x', '', colnames(D18[[1]])))
colnames(D18[[2]]) <- paste0('D18_', gsub('IEOdev_DSP_d18:|x', '', colnames(D18[[2]])))
D18[[1]] <- D18[[1]][VariableFeatures(X),]
D18[[2]] <- D18[[2]][VariableFeatures(X),]

D24 <- read.loom.matrices("IEOdev_DSP_d24.loom")
colnames(D24[[1]]) <- paste0('D24_', gsub('IEOdev_DSP_d24:|x', '', colnames(D24[[1]])))
colnames(D24[[2]]) <- paste0('D24_', gsub('IEOdev_DSP_d24:|x', '', colnames(D24[[2]])))
D24[[1]] <- D24[[1]][VariableFeatures(X),]
D24[[2]] <- D24[[2]][VariableFeatures(X),]

D30 <- read.loom.matrices("IEOdev_DSP_d30.loom")
colnames(D30[[1]]) <- paste0('D30_', gsub('IEOdev_DSP_d30:|x', '', colnames(D30[[1]])))
colnames(D30[[2]]) <- paste0('D30_', gsub('IEOdev_DSP_d30:|x', '', colnames(D30[[2]])))
D30[[1]] <- D30[[1]][VariableFeatures(X),]
D30[[2]] <- D30[[2]][VariableFeatures(X),]

D36 <- read.loom.matrices("IEOdev_DSP_d36.loom")
colnames(D36[[1]]) <- paste0('D36_', gsub('IEOdev_DSP_d36:|x', '', colnames(D36[[1]])))
colnames(D36[[2]]) <- paste0('D36_', gsub('IEOdev_DSP_d36:|x', '', colnames(D36[[2]])))
D36[[1]] <- D36[[1]][VariableFeatures(X),]
D36[[2]] <- D36[[2]][VariableFeatures(X),]

# Merging spliced (S) and unspliced (U) matrices
S <- cbind(D0[[1]], D3[[1]], D6[[1]], D8[[1]], D10[[1]], 
           D13[[1]], D18[[1]], D24[[1]], D30[[1]], D36[[1]])
U <- cbind(D0[[2]], D3[[2]], D6[[2]], D8[[2]], D10[[2]], 
           D13[[2]], D18[[2]], D24[[2]], D30[[2]], D36[[2]])

# Randomly selecting 10K Cells for testing purposes
set.seed(1)
SC <- sample(seq_len(ncol(S)), 10000)
S <- S[,SC]
U <- U[,SC]

# Matching cells with the original Seurat object
mdX <- mdX[intersect(colnames(S), rownames(mdX)),]
X <- X[,mdX$BC]

# Getting the low dimensional representation to be used
E <- X@reductions$umap@cell.embeddings
rownames(E) <- rownames(mdX)
S <- S[,rownames(mdX)]
U <- U[,rownames(mdX)]
D <- X@reductions$pca@cell.embeddings
rownames(D) <- rownames(E)
D <- as.dist(1-cor(t(D)))

# Computing transitions
velocity <- gene.relative.velocity.estimates(S, U, deltaT = 1, kCells = 20, cell.dist = D, fit.quantile = 0.02)

# Computing arrows
arrows <- show.velocity.on.embedding.cor(E,velocity, n = 20, grid.n = 50, show.grid.flow = TRUE, return.details = TRUE)
dev.off()

# Plotting
umap_arrows <- arrows$garrows %>%
    as.data.frame() %>%
    mutate(x2 = x0 + (x1 - x0) * 2,
           y2 = y0 + (y1 - y0) * 2)

png('t.png', width = 2000, height = 2000, res = 300)
E0 <- as.data.frame(E0)
ggplot(E0, aes(UMAP_1, UMAP_2)) +
geom_point(alpha = 0.25, pch = 16, color = 'gray80')+
geom_curve(data = umap_arrows, curvature = 0.2, angle = 45,
                 aes(x = x0, xend = x2, y = y0, yend = y2),
                 size = .1,
                 arrow = arrow(length = unit(1, "points"), type = "closed"),
                 colour = "grey20", alpha = 0.8) +
theme_light() +
theme(panel.grid = element_blank(), panel.border = element_blank()) +
xlab('UMAP 1') + ylab('UMAP 2')
dev.off()


save(velocity, arrows , file = 'scVelocity_SS3.RData')
