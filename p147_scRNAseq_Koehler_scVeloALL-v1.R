# Libraries
snapshot <- "/lab-share/RC-Data-Science-e2/Public/DST_software/pipeline_R_snapshot/r4.1-20220803-bch-dst-v1/"
.libPaths(snapshot)

library(Seurat)
library(dplyr)
library(ggplot2)
library(velocyto.R)
library(velociraptor)

# Setting WD
setwd('/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler')

# # Loading Seurat Object
# load('results/D100/P147.RData')
# E0 <- X@reductions$umap@cell.embeddings
# mdX <- X@meta.data
# mdX <- as.data.frame.array(mdX)
# mdX$BC <- rownames(mdX)
# rownames(mdX) <- paste0(mdX$orig.ident, '_', gsub('-[[:print:]]+$','',rownames(mdX)))

# # Reading Matrices
# D0 <- read.loom.matrices("IEOdev_DSP_d0.loom")
# colnames(D0[[1]]) <- paste0('D0_', gsub('IEOdev_DSP_d0:|x', '', colnames(D0[[1]])))
# colnames(D0[[2]]) <- paste0('D0_', gsub('IEOdev_DSP_d0:|x', '', colnames(D0[[2]])))
# D0[[1]] <- D0[[1]][VariableFeatures(X),]
# D0[[2]] <- D0[[2]][VariableFeatures(X),]

# D3 <- read.loom.matrices("IEOdev_DSP_d3.loom")
# colnames(D3[[1]]) <- paste0('D3_', gsub('IEOdev_DSP_d3:|x', '', colnames(D3[[1]])))
# colnames(D3[[2]]) <- paste0('D3_', gsub('IEOdev_DSP_d3:|x', '', colnames(D3[[2]])))
# D3[[1]] <- D3[[1]][VariableFeatures(X),]
# D3[[2]] <- D3[[2]][VariableFeatures(X),]

# D6 <- read.loom.matrices("IEOdev_DSP_d6.loom")
# colnames(D6[[1]]) <- paste0('D6_', gsub('IEOdev_DSP_d6:|x', '', colnames(D6[[1]])))
# colnames(D6[[2]]) <- paste0('D6_', gsub('IEOdev_DSP_d6:|x', '', colnames(D6[[2]])))
# D6[[1]] <- D6[[1]][VariableFeatures(X),]
# D6[[2]] <- D6[[2]][VariableFeatures(X),]

# D8 <- read.loom.matrices("IEOdev_DSP_d8.loom")
# colnames(D8[[1]]) <- paste0('D8_', gsub('IEOdev_DSP_d8:|x', '', colnames(D8[[1]])))
# colnames(D8[[2]]) <- paste0('D8_', gsub('IEOdev_DSP_d8:|x', '', colnames(D8[[2]])))
# D8[[1]] <- D8[[1]][VariableFeatures(X),]
# D8[[2]] <- D8[[2]][VariableFeatures(X),]

# D10 <- read.loom.matrices("IEOdev_DSP_d10N.loom")
# colnames(D10[[1]]) <- paste0('D10_', gsub('IEOdev_DSP_d10N:|x', '', colnames(D10[[1]])))
# colnames(D10[[2]]) <- paste0('D10_', gsub('IEOdev_DSP_d10N:|x', '', colnames(D10[[2]])))
# D10[[1]] <- D10[[1]][VariableFeatures(X),]
# D10[[2]] <- D10[[2]][VariableFeatures(X),]

# D13 <- read.loom.matrices("IEOdev_DSP_d13N.loom")
# colnames(D13[[1]]) <- paste0('D13_', gsub('IEOdev_DSP_d13N:|x', '', colnames(D13[[1]])))
# colnames(D13[[2]]) <- paste0('D13_', gsub('IEOdev_DSP_d13N:|x', '', colnames(D13[[2]])))
# D13[[1]] <- D13[[1]][VariableFeatures(X),]
# D13[[2]] <- D13[[2]][VariableFeatures(X),]

# D18 <- read.loom.matrices("IEOdev_DSP_d18.loom")
# colnames(D18[[1]]) <- paste0('D18_', gsub('IEOdev_DSP_d18:|x', '', colnames(D18[[1]])))
# colnames(D18[[2]]) <- paste0('D18_', gsub('IEOdev_DSP_d18:|x', '', colnames(D18[[2]])))
# D18[[1]] <- D18[[1]][VariableFeatures(X),]
# D18[[2]] <- D18[[2]][VariableFeatures(X),]

# D24 <- read.loom.matrices("IEOdev_DSP_d24.loom")
# colnames(D24[[1]]) <- paste0('D24_', gsub('IEOdev_DSP_d24:|x', '', colnames(D24[[1]])))
# colnames(D24[[2]]) <- paste0('D24_', gsub('IEOdev_DSP_d24:|x', '', colnames(D24[[2]])))
# D24[[1]] <- D24[[1]][VariableFeatures(X),]
# D24[[2]] <- D24[[2]][VariableFeatures(X),]

# D30 <- read.loom.matrices("IEOdev_DSP_d30.loom")
# colnames(D30[[1]]) <- paste0('D30_', gsub('IEOdev_DSP_d30:|x', '', colnames(D30[[1]])))
# colnames(D30[[2]]) <- paste0('D30_', gsub('IEOdev_DSP_d30:|x', '', colnames(D30[[2]])))
# D30[[1]] <- D30[[1]][VariableFeatures(X),]
# D30[[2]] <- D30[[2]][VariableFeatures(X),]

# D36 <- read.loom.matrices("IEOdev_DSP_d36.loom")
# colnames(D36[[1]]) <- paste0('D36_', gsub('IEOdev_DSP_d36:|x', '', colnames(D36[[1]])))
# colnames(D36[[2]]) <- paste0('D36_', gsub('IEOdev_DSP_d36:|x', '', colnames(D36[[2]])))
# D36[[1]] <- D36[[1]][VariableFeatures(X),]
# D36[[2]] <- D36[[2]][VariableFeatures(X),]

# # Merging spliced (S) and unspliced (U) matrices
# S <- cbind(D0[[1]], D3[[1]], D6[[1]], D8[[1]], D10[[1]], 
#            D13[[1]], D18[[1]], D24[[1]], D30[[1]], D36[[1]])
# U <- cbind(D0[[2]], D3[[2]], D6[[2]], D8[[2]], D10[[2]], 
#            D13[[2]], D18[[2]], D24[[2]], D30[[2]], D36[[2]])

# save(S,U,X, mdX, file = 'P147VelocityData.RData')
load('P147VelocityData.RData')

# Randomly selecting 10K Cells for testing purposes
set.seed(1)
SC <- sample(seq_len(ncol(S)), 20000)
S <- S[,SC]
U <- U[,SC]

# Matching cells with the original Seurat object
mdX <- mdX[intersect(colnames(S), rownames(mdX)),]
X <- X[,mdX$BC]

# Getting the low dimensional representation to be used
E <- X@reductions$umap@cell.embeddings
rownames(E) <- rownames(mdX)
X <- X[,mdX$BC]

# # Velocity
# ad <- import("anndata", convert = FALSE)
# dfobs <- X@meta.data
# rownames(dfobs) <- rownames(mdX)
# dfvar <- data.frame(id = rownames(X[VariableFeatures(X),]))
# rownames(dfvar) <- rownames(X[VariableFeatures(X),])
# adata <- ad$AnnData(
#   X=t(X@assays$RNA@counts[VariableFeatures(X),mdX$BC]),
#   obs=dfobs,
#   var=dfvar,
#   layers=list('spliced'=t(S[,rownames(mdX)]), 'unspliced'=t(U[,rownames(mdX)])),
#   obsm=list('X_umap'=E, 'X_pca'=X@reductions$pca@cell.embeddings[,1:2]) 
#   )

# scv$pp$filter_genes(adata) ## filter
# scv$pp$moments(adata) ## normalize and compute moments
# scv$tl$recover_dynamics(adata) ## model
# scv$tl$velocity(adata, mode='dynamical')
# scv$tl$velocity_graph(adata)
# scv$tl$velocity_pseudotime(adata)

# save(adata, file = 'vTest.RData')
# load('vTest.RData')

# library(velociraptor)
# .make_np_friendly <- function(x) {
#     if (is_sparse(x)) {
#         as(x, "dgCMatrix")
#     } else {
#         as.matrix(x)
#     }
# }

# library(scRNAseq)
# sce <- HermannSpermatogenesisData()
# sce
# sce <- sce[, 1:500]
# library(scuttle)
# sce <- logNormCounts(sce, assay.type=1)

# library(scran)
# dec <- modelGeneVar(sce)
# top.hvgs <- getTopHVGs(dec, n=2000)

# library(velociraptor)
# velo.out <- scvelo(sce, subset.row=top.hvgs, assay.X="spliced")
# velo.out

# library(scater)

# set.seed(100)
# sce <- runPCA(sce, subset_row=top.hvgs)
# sce <- runTSNE(sce, dimred="PCA")

# sce$velocity_pseudotime <- velo.out$velocity_pseudotime
# embedded <- embedVelocity(reducedDim(sce, "TSNE"), velo.out)
# grid.df <- gridVectors(reducedDim(sce, "TSNE"), embedded)
# rownames(embedded) <- colnames(sce)
# colnames(embedded) <- c('tsne1', 'tsne2')


# library(scuttle)
# set.seed(42)
# sce1 <- mockSCE(ncells = 100, ngenes = 500)
# sce2 <- mockSCE(ncells = 100, ngenes = 500)

# datlist <- list(X=counts(sce1), spliced=counts(sce1), unspliced=counts(sce2))

# out <- scvelo(datlist, mode = "dynamical")

# em <- embedVelocity(reducedDim(out, 1), out)[,1:2]

# plotVelocityStream(velo.out, embedded, color.streamlines = TRUE)
# dev.off()
dimnames(E) <- NULL
datList <- list(
    X=X@assays$RNA@counts[VariableFeatures(X),mdX$BC], 
    spliced=S[,rownames(mdX)], 
    unspliced=U[,rownames(mdX)])

P147ALL <- scvelo(datList, mode = "dynamical")
P147EM <- embedVelocity(as.matrix(E), P147ALL)
P147G1 <- gridVectors(E, P147EM, resolution = 50, as.data.frame = FALSE, return.intermediates = TRUE)
P147G2 <- gridVectors(E, P147EM, resolution = 50)

ggplot(as.data.frame(E), aes(V1,V2)) +
geom_segment(data=P147G2, mapping=aes(x=start.1, y=start.2, 
        xend=end.1, yend=end.2), arrow=arrow(length=unit(0.05, "inches")))
dev.off()

save(P147ALL, P147EM, E, P147G1, P147G2, file = 'P147_scVeloResults.RData')


load('P147_scVeloResults.RData')

library(scater)
reducedDim(P147ALL, "UMAP") <- E

png('AAAA1.png', width = 1500, height = 1500, res = 300)
M <- velociraptor::plotVelocityStream(
    P147ALL, 
    P147EM, 
    color_by = matchColors2[levels(Idents(P147D100))],
    use.dimred = 2, 
    grid.resolution = 50,
    scale = TRUE,
    stream.L = 3,
    stream.res = 1,
    stream.width = 1,
    arrow.length = 0.2,
    arrow.angle = 25
    ) +
    geom_point(size = 0.1, color = 'gray80')
M$layers <- list(M$layers[[3]], M$layers[[2]])
print(M)
dev.off()

