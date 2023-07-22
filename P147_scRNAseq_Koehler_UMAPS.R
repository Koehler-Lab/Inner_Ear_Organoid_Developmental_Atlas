snapshot <- "/lab-share/RC-Data-Science-e2/Public/DST_software/pipeline_R_snapshot/r4.1-20220803-bch-dst-v1/"
.libPaths(snapshot)
library(Seurat)
library(ggplot2)
library(ggrepel)

setwd('/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler/')
# Loading data
load('/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler/P147D100_Annotated.RData')

tpColors <- tenColors <- c('#a9d6e5', '#89c2d9', '#61a5c2', '#468faf', '#2c7da0', '#2a6f97', '#014f86', '#01497c', '#013a63', '#012a4a')
ctColors <- c('#AF0508', '#F1070B', '#F94144', 
'#9A3C09', '#D4520C', '#F3722C', 
'#9B5805', '#D67907', '#F8961E', 
'#B48007', '#F6AF0A', '#F9C74F', 
'#527735', '#70A349', '#90BE6D', 
'#26604F', '#35846C', '#43AA8B', 
'#314352', '#445C71', '#577590')

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

for(tp in unique((P147D100$orig.ident))){
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
  labs(title = tp)
 
  png(paste0('UMAP-',tp,'.png'), width = 1500, height = 1500, res = 300)
  print(P)
  dev.off()
}
