# singularity exec --bind /lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler:/data/,/lab-share/RC-Data-Science-e2/Public/DST_software/pipeline_R_snapshot/:/snapshot/ r-monocle3:1.0.0--r42h9f5acd7_3 R

snapshot <- "/snapshot/r4.2.3-20230324-bch-dst-v1"
.libPaths(snapshot)
library(Seurat)
library(monocle3)
library(ggplot2)

# Loading data
load('/data/P147D100_Annotated.RData')

tpColors <- tenColors <- c('#a9d6e5', '#89c2d9', '#61a5c2', '#468faf', '#2c7da0', '#2a6f97', '#014f86', '#01497c', '#013a63', '#012a4a')
ctColors <- c('#AF0508', '#F1070B', '#F94144', 
'#9A3C09', '#D4520C', '#F3722C', 
'#9B5805', '#D67907', '#F8961E', 
'#B48007', '#F6AF0A', '#F9C74F', 
'#527735', '#70A349', '#90BE6D', 
'#26604F', '#35846C', '#43AA8B', 
'#314352', '#445C71', '#577590')

Label1 <- unique(P147D100$CT)
Label2 <- unique(P147D100$ID)
matchColors2 <- ctColors
names(matchColors2) <- Label2

Idents(P147D100) <- factor(Idents(P147D100), Label2)

getTrajectory <- function(X, nDim = 30, transferUMAP = FALSE, joinALL = FALSE){
  geneInfo <- data.frame(gene_short_name = rownames(X))
  rownames(geneInfo) <- rownames(X)
  cellInfo <- X@meta.data
  cds <- new_cell_data_set(X@assays$RNA@counts,
                          cell_metadata = cellInfo,
                          gene_metadata = geneInfo)
  cds <- preprocess_cds(cds, num_dim = nDim)
  cds <- align_cds(cds, alignment_group = "orig.ident")
  cds <- reduce_dimension(cds)
  if(transferUMAP){
    reducedDims(cds)[['UMAP']] <- X@reductions$umap@cell.embeddings
  }
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds, use_partition = !joinALL)
  attr(X, 'monocle3') <- cds
  return(X)
}

orderCells <- function(X, rootNode){
  attr(X, 'monocle3') <- order_cells(attr(X, 'monocle3'), root_pr_nodes = rootNode)
  return(X)
}

P147D100 <- getTrajectory(P147D100, transferUMAP = TRUE, joinALL = TRUE)

png('/data/P147-Monocle3_TransferredUMAP.png', width = 2000*1.4, height = 1500*1.5, res = 300)
P2 <- plot_cells(attr(P147D100, 'monocle3'), label_principal_points = TRUE,  color_cells_by = "ID") +
theme_void() +
guides(color= guide_legend("Cell Type", override.aes = list(size=5), ncol = 2)) +
scale_color_manual(values = matchColors2[levels(Idents(P147D100))])
print(P2)
dev.off()

P147D100 <- getTrajectory(P147D100, transferUMAP = FALSE, joinALL = TRUE)

png('/data/P147-Monocle3.png', width = 2000*1.4, height = 1500*1.5, res = 300)
P2 <- plot_cells(attr(P147D100, 'monocle3'), label_principal_points = TRUE,  color_cells_by = "ID") +
theme_void() +
guides(color= guide_legend("Cell Type", override.aes = list(size=5), ncol = 2)) +
scale_color_manual(values = matchColors2[levels(Idents(P147D100))])
print(P2)
dev.off()
