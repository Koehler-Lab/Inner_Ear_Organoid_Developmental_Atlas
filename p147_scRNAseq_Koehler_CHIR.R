require(symphony)
require(Seurat)
snapshot <- "/lab-share/RC-Data-Science-e2/Public/DST_software/pipeline_R_snapshot/r4.1-20220803-bch-dst-v1/"
.libPaths(snapshot)
source('/lab-share/RC-Data-Science-e2/Public/Daniel/rFunctions.R')

setwd('/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler')

# Loading Data
load('/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler/P147D100_Annotated.RData')
load('/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler/CHIR.RData')

doTransferLabel <- function(X, Y){
# Getting metadata from the reference
MD <- data.frame(B = X@meta.data[,1], L = Idents(X))

# Removing previous instances of the function
if(file.exists('.UMAP')){
    file.remove('.UMAP')
}

# Generating a Symphony Reference
S <- buildReference(X@assays$RNA@counts, MD, 
                        vars = 'B', 
                        do_umap = TRUE, 
                        verbose = TRUE, 
                        d = 30,
                        save_uwot_path = '.UMAP')

# Replacing embeeding by the one in the Seurat object
S$umap$embedding[,1] <- X@reductions$umap@cell.embeddings[,1]
S$umap$embedding[,2] <- X@reductions$umap@cell.embeddings[,2]
U <- uwot::load_uwot('.UMAP')
U$embedding[,1] <- X@reductions$umap@cell.embeddings[,1]
U$embedding[,2] <- X@reductions$umap@cell.embeddings[,2]
file.remove('.UMAP')
U <- uwot::save_uwot(U, file = '.UMAP', verbose = FALSE)

# Querying the new data into the generated reference
set.seed(1)
Q <- mapQuery(Y@assays$RNA@counts, Y@meta.data,
                      ref_obj = S)

# Transfering Labels
Q <- knnPredict(Q, S, S$meta_data$L, k = 5)

# Adding new metadata to the Y object
Y$transferedCellType <- Q$meta_data$cell_type_pred_knn
Y$labelProbability <- Q$meta_data$cell_type_pred_knn_prob
Idents(Y) <- Y$transferedCellType
Y@reductions$umap@cell.embeddings[,1] <- Q$umap[,1]
Y@reductions$umap@cell.embeddings[,2] <- Q$umap[,2]

# Returning annotated object
return(Y)
}

# Applying the function
CHIR <- doTransferLabel(P147D100[,P147D100$orig.ident %in% c('D10', 'D13')], CHIR)
save(CHIR, file = 'CHIR_Annotated.RData')
