library(Seurat)
library(symphony)
snapshot <- "/lab-share/RC-Data-Science-e2/Public/DST_software/pipeline_R_snapshot/r4.1-20220803-bch-dst-v1/"
.libPaths(snapshot)
source('/lab-share/RC-Data-Science-e2/Public/Daniel/rFunctions.R')
library(ggplot2)
library(ggrepel)

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


load('P147D100_Annotated.RData')
Label1 <- unique(P147D100$CT)
Label2 <- unique(P147D100$ID)
matchColors2 <- ctColors
names(matchColors2) <- Label2
Idents(P147D100) <- factor(Idents(P147D100), Label2)

LabelPosition <- lapply(split(as.data.frame(P147D100@reductions$umap@cell.embeddings), Idents(P147D100)), function(X){apply(X,2,median)})
Labels <- names(LabelPosition)
LabelPosition <- do.call(rbind.data.frame, LabelPosition)
colnames(LabelPosition) <- c('X', 'Y')
LabelPosition$Label <- Labels


D6 <- P147D100[,P147D100$orig.ident %in% 'D6']
Y <- readRDS('/lab-share/RC-Data-Science-e2/Public/Daniel/p147-scRNAseq-Koehler/data/d6_User_Comp3.rds')

Y <- doTransferLabel(D6, Y, TRUE)
Y$user <- factor(Y$user, c('Matt', 'Jiyoon (WA25)', 'Jiyoon (DSP)'))

D <- P147D100@reductions$umap@cell.embeddings
D <- as.data.frame.array(D)
D$ID <- P147D100$ID
D$Color <- matchColors2[D$ID]
D$TP <- P147D100$orig.ident
D$Color <- 'gray95'
D <- D[!P147D100$orig.ident %in% 'D6',]
A <- Y@reductions$umap@cell.embeddings
A <- as.data.frame.array(A)
A$ID <- Idents(Y)
A$Color <- matchColors2[as.vector(A$ID)]
A$TP <- Y$user
A <- A[A$TP %in% 'Matt',]
D <- rbind(D,A)
sLabels <- table(D$ID[D$Color != 'gray95'])
sLabels <- names(sLabels[sLabels >= 50])
UMAPD6 <- ggplot(D, aes(UMAP_1, UMAP_2)) +
geom_point(color = D$Color, size = 0.01) +
theme_void() +
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
labs(title = 'Experimenter 1 (DSP)')

png('D6-E1.png', width = 1500, height = 1500, res = 300)
print(UMAPD6)
dev.off()


D <- P147D100@reductions$umap@cell.embeddings
D <- as.data.frame.array(D)
D$ID <- P147D100$ID
D$Color <- matchColors2[D$ID]
D$TP <- P147D100$orig.ident
D$Color <- 'gray95'
D <- D[!P147D100$orig.ident %in% 'D6',]
A <- Y@reductions$umap@cell.embeddings
A <- as.data.frame.array(A)
A$ID <- Idents(Y)
A$Color <- matchColors2[as.vector(A$ID)]
A$TP <- Y$user
A <- A[A$TP %in% 'Jiyoon (DSP)',]
D <- rbind(D,A)
sLabels <- table(D$ID[D$Color != 'gray95'])
sLabels <- names(sLabels[sLabels >= 50])
UMAPD6 <- ggplot(D, aes(UMAP_1, UMAP_2)) +
geom_point(color = D$Color, size = 0.01) +
theme_void() +
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
labs(title = 'Experimenter 2 (DSP)')

png('D6-E2.png', width = 1500, height = 1500, res = 300)
print(UMAPD6)
dev.off()

D <- P147D100@reductions$umap@cell.embeddings
D <- as.data.frame.array(D)
D$ID <- P147D100$ID
D$Color <- matchColors2[D$ID]
D$TP <- P147D100$orig.ident
D$Color <- 'gray95'
D <- D[!P147D100$orig.ident %in% 'D6',]
A <- Y@reductions$umap@cell.embeddings
A <- as.data.frame.array(A)
A$ID <- Idents(Y)
A$Color <- matchColors2[as.vector(A$ID)]
A$TP <- Y$user
A <- A[A$TP %in% 'Jiyoon (WA25)',]
D <- rbind(D,A)
sLabels <- table(D$ID[D$Color != 'gray95'])
sLabels <- names(sLabels[sLabels >= 50])
UMAPD6 <- ggplot(D, aes(UMAP_1, UMAP_2)) +
geom_point(color = D$Color, size = 0.01) +
theme_void() +
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
labs(title = 'Experimenter 2 (WA25)')

png('D6-E3.png', width = 1500, height = 1500, res = 300)
print(UMAPD6)
dev.off()

Y$user <- factor(Y$user, c('Matt', 'Jiyoon (DSP)', 'Jiyoon (WA25)'))

CP <- table(Y$user, Idents(Y))
CP <- CP/rowSums(CP) * 100
CP <- reshape2::melt(CP)
levels(CP[,1]) <- c('Experimenter 1 (DSP)', 'Experimenter 2 (DSP)', 'Experimenter 2 (WA25)')
CP[,1] <- factor(CP[,1], rev(c('Experimenter 1 (DSP)', 'Experimenter 2 (DSP)', 'Experimenter 2 (WA25)')))
colnames(CP) <- c('E', 'CT', 'P')

png('D6-CP.png', width = 3000, height = 750, res = 300)
P <- ggplot(CP, aes(P, E, fill = CT)) +
geom_bar(stat = 'identity') +
scale_fill_manual(values = matchColors2[levels(CP$CT)]) +
theme_light() +
theme(panel.border = element_blank(), 
panel.grid = element_blank(),
legend.position = 'bottom') +
guides(fill=guide_legend(ncol=1)) +
xlab('Proportion (%)') +
ylab(NULL) +
labs(fill = 'Cell Type') +
guides(fill=guide_legend(nrow=2,byrow=TRUE))
print(P)
dev.off()