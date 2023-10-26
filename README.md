##### Project: Single cell analysis of the role of CD47 on innate and adaptive immune cells during infection #########

###########################################################################################################################################
#### Below are codes used to generate main figures and supplementary figures.##############################################################
###########################################################################################################################################
 
#### Figure 1 ############

library(Seurat)
library(SeuratObject)
library(dplyr)
library(pheatmap)
library(BiocManager)
library(multtest)
library(metap)
library(limma)
library(ggplot2)
library(cowplot)
lamin_all <- readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 1/lamin_all.rds")
lamin_all <- ScaleData(object = lamin_all, verbose = FALSE)
lamin_all <- RunPCA(object = lamin_all, npcs = 30, verbose = FALSE)
lamin_all <- RunUMAP(object = lamin_all, reduction = "pca", dims = 1:20) 
lamin_all <- FindNeighbors(object = lamin_all, reduction = "pca", dims = 1:20)
lamin_all <- FindClusters(lamin_all, resolution = 0.3)
DimPlot(lamin_all, reduction = "umap", label = TRUE, repel = TRUE , raster=FALSE, label.size = 8)+  theme(text = element_text(size = 15)) + NoLegend()
new.cluster.ids <-c("CD8 T cells","CD4 T cells", "B cells","Monocytes","NK cells","NK cells","CD4 T cells", "Platelet", "Dendritic cells","CD4 T cells","B cells","Monocytes","B cells")
names(new.cluster.ids) <- levels(lamin_all)
lamin_all <- RenameIdents(lamin_all, new.cluster.ids)
DimPlot(lamin_all, reduction = "umap", label = TRUE, repel = TRUE , raster=FALSE, label.size = 8)+  theme(text = element_text(size = 15)) + NoLegend()
lamin_all <- RunTSNE(object = lamin_all)
DimPlot(object = lamin_all, reduction = "tsne", label = TRUE, repel = TRUE,  label.size = 12, pt.size = 2) + theme(text = element_text(size = 25)) + NoLegend()
#ggsave("pic1.png", height = 13, width = 11)
DefaultAssay(lamin_all) <- "RNA"
features<- c('CD47', "SIRPA" ) 
DotPlot(lamin_all, features = features, cols = c("blue", "red", "green"), dot.scale = 20, group.by = "condition") + RotatedAxis() + theme(axis.text = element_text(size = 20)) + theme(axis.title.x = element_blank())+ theme(axis.title.y = element_blank())
FeaturePlot(lamin_all, features = c("CD47"), raster=FALSE, reduction = "tsne", pt.size = 1, cols = c( "darkgray", "blue", "darkred"))
FeaturePlot(lamin_all, features = c("SIRPA"), raster=FALSE, reduction = "tsne", pt.size = 1, cols = c( "darkgray", "blue", "darkred"))
#ggsave("pic1.png", height = 6, width = 6 )
VlnPlot(lamin_all, features = c("SIRPA"), raster=FALSE, pt.size = 1) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + NoLegend() + theme(axis.title.y = element_text(size = 25)) + theme(title = element_text(size = 20))
VlnPlot(lamin_all, features = c("CD47"), pt.size = 0, raster=FALSE, split.by = "condition")

RidgePlot(lamin_all, features = c( "CD47"), ncol = NULL) + NoLegend() + theme(title = element_text(size = 20)) + theme(axis.text.y = element_text(size = 25) ) + theme(axis.title.y = element_blank()) +theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20) )
ggsave("pic1.png", height = 11, width = 8 )


#### CD47 categorization of Immune cells
DefaultAssay(lamin_all) <- "RNA"
lamin_all_CD47lo <- subset(lamin_all, subset = CD47<0.5)
lamin_all_CD47int <- subset(lamin_all, subset = CD47>0.5 & CD47<1.5)
lamin_all_CD47hi <- subset(lamin_all, subset = CD47>1.5)

saveRDS(lamin_all_CD47lo, file = "C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 1/lamin_all_CD47lo.rds")
saveRDS(lamin_all_CD47int, file = "C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 1/lamin_all_CD47int.rds")
saveRDS(lamin_all_CD47hi, file = "C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 1/lamin_all_CD47hi.rds")

#integration
lamin_all_CD47lo <- readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 1/lamin_all_CD47lo.rds")
lamin_all_CD47lo$CD47_expression <- "CD47 Low"
lamin_all_CD47lo <- NormalizeData(lamin_all_CD47lo, normalization.method = "LogNormalize", scale.factor = 10000)
lamin_all_CD47lo <- FindVariableFeatures(lamin_all_CD47lo, selection.method = "vst", nfeatures = 2000)
lamin_all_CD47lo <- ScaleData(object = lamin_all_CD47lo, verbose = FALSE)

lamin_all_CD47int <- readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 1/lamin_all_CD47int.rds")
lamin_all_CD47int$CD47_expression <- "CD47 Inter."
lamin_all_CD47int <- NormalizeData(lamin_all_CD47int, normalization.method = "LogNormalize", scale.factor = 10000)
lamin_all_CD47int <- FindVariableFeatures(lamin_all_CD47int, selection.method = "vst", nfeatures = 2000)
lamin_all_CD47int <- ScaleData(object = lamin_all_CD47int, verbose = FALSE)

lamin_all_CD47hi <- readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 1/lamin_all_CD47hi.rds")
lamin_all_CD47hi$CD47_expression <- "CD47 High"
lamin_all_CD47hi <- NormalizeData(lamin_all_CD47hi, normalization.method = "LogNormalize", scale.factor = 10000)
lamin_all_CD47hi <- FindVariableFeatures(lamin_all_CD47hi, selection.method = "vst", nfeatures = 2000)
lamin_all_CD47hi <- ScaleData(object = lamin_all_CD47hi, verbose = FALSE)

anchors <- FindIntegrationAnchors(object.list = list(lamin_all_CD47lo, lamin_all_CD47int, lamin_all_CD47hi ), dims = 1:20)
immune_combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = immune_combined) <- "integrated"
saveRDS(immune_combined, file = "C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 1/immune_combined.rds")
immune_combined <- readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 1/immune_combined.rds")
immune_combined <- ScaleData(object = immune_combined, verbose = FALSE)
immune_combined <- RunPCA(object = immune_combined, npcs = 30, verbose = FALSE)
immune_combined <- RunUMAP(object = immune_combined, reduction = "pca", dims = 1:20) 
immune_combined <- FindNeighbors(object = immune_combined, reduction = "pca", dims = 1:20)
immune_combined <- FindClusters(immune_combined, resolution = 0.3)
DimPlot(immune_combined, reduction = "umap", label = TRUE, repel = TRUE , raster=FALSE, label.size = 8)+  theme(text = element_text(size = 15)) + NoLegend()
immune_combined <- RunTSNE(object = immune_combined, check_duplicates = FALSE)
DimPlot(object = immune_combined, reduction = "tsne", label = TRUE, repel = TRUE,  label.size = 10, pt.size = 2) + theme(text = element_text(size = 25)) + NoLegend()

Idents(immune_combined) <- "CD47_expression"
DefaultAssay(immune_combined) <- "RNA"
FeaturePlot(immune_combined, features = c("CD47"), raster=FALSE,label = F, reduction = "tsne", pt.size = 1, cols = c( "darkgray", "blue", "darkred"), split.by = "CD47_expression" )
FeaturePlot(immune_combined, features = c("SIRPA"), raster=FALSE, reduction = "tsne", pt.size = 1, cols = c( "gray96", "red3"))
#ggsave("pic1.png", height = 8, width = 10 )

act_cytx <- c("CD47","TLR7","STAT1", "STAT2", "IRF7", "IRF9", "CD69", "ISG15", "LY6E", "IFI6","MX1", "IFIT2", "IFIT3", "IFIT5", "OASL", "IFNG","GZMB","GZMA","GZMH", "GZMK",  "NKG7", "IL2RB","PRF1", "GNLY", "KLRB1", "PDCD1", "LAG3", "TIGIT")
DotPlot(immune_combined, features = act_cytx, cols = c("blue", "red", "green"), dot.scale = 20) + RotatedAxis() + theme(axis.text = element_text(size = 20)) + theme(axis.title.x = element_blank())+ theme(axis.title.y = element_blank())

ave <- AverageExpression(immune_combined)$RNA
gene <- intersect(rownames(ave), act_cytx)
ave_lamin <- ave[gene,]
colfunc<-colorRampPalette(c("blue3","gray96","red3"))
pheatmap(ave_lamin, col=(colfunc(100)), scale="row",border_color=NA, fontsize = 20,cluster_cols = FALSE)

VlnPlot(immune_combined, features = c("CD69"), pt.size = 0, raster=FALSE)+NoLegend()+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20))
VlnPlot(immune_combined, features = c("LY6E"), pt.size = 0, raster=FALSE)+NoLegend()+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20))
VlnPlot(immune_combined, features = c("GZMB"), pt.size = 0, raster=FALSE)+NoLegend()+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20))
VlnPlot(immune_combined, features = c("KLRB1"), pt.size = 0, raster=FALSE)+NoLegend()+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20))
VlnPlot(immune_combined, features = c("GNLY"), pt.size = 0, raster=FALSE)+NoLegend()+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20))
#ggsave("pic1.png", height = 6, width = 5)


##### Figure 2: Monocytes analysis ############

library(Seurat)
library(SeuratObject)
library(dplyr)
library(pheatmap)
library(BiocManager)
library(multtest)
library(metap)
library(limma)
library(ggplot2)
library(cowplot)
monocyte_dc <- readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 2/monocyte_dc.rds")
monocyte_dc <- ScaleData(object = monocyte_dc, verbose = FALSE)
monocyte_dc <- RunPCA(object = monocyte_dc, npcs = 30, verbose = FALSE)
monocyte_dc <- RunUMAP(object = monocyte_dc, reduction = "pca", dims = 1:20) 
monocyte_dc <- FindNeighbors(object = monocyte_dc, reduction = "pca", dims = 1:20)
monocyte_dc <- FindClusters(monocyte_dc, resolution = 0.1)
DimPlot(monocyte_dc, reduction = "umap", label = TRUE, repel = TRUE , raster=FALSE, label.size = 8)+  theme(text = element_text(size = 15)) + NoLegend()
DefaultAssay(monocyte_dc) <- "RNA"
new.cluster.ids <-c("CD1c DC","CD14 Monocytes","CD14 Monocytes","CD14 Monocytes","CD1c DC","CD1c DC","CD16 Monocytes", "AS DC", "CD1c DC", "CLEC9A DC", "pDC", "CD1c DC")
names(new.cluster.ids) <- levels(monocyte_dc)
monocyte_dc <- RenameIdents(monocyte_dc, new.cluster.ids)
monocyte_dc <- RunTSNE(object = monocyte_dc)
DimPlot(object = monocyte_dc, reduction = "tsne", label = TRUE, repel = TRUE,  label.size = 10, pt.size = 0) + theme(text = element_text(size = 25)) + NoLegend()
#ggsave("pic1.png", height = 7, width = 10)
VlnPlot(monocyte_dc, features = c("CD47"), pt.size = 0, raster=FALSE) + NoLegend() + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20))


#figure 2a
lamin_Monocytes <- readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 2/Monocytes/lamin_Monocytes.rds")
lamin_Monocytes <- FindVariableFeatures(lamin_Monocytes, selection.method = "vst", nfeatures = 2000)
lamin_Monocytes<- ScaleData(object = lamin_Monocytes, verbose = FALSE)
lamin_Monocytes<- RunPCA(object = lamin_Monocytes, npcs = 30, verbose = FALSE)
lamin_Monocytes<- RunUMAP(object = lamin_Monocytes, reduction = "pca", dims = 1:20) 
lamin_Monocytes<- FindNeighbors(object = lamin_Monocytes, reduction = "pca", dims = 1:20)
lamin_Monocytes<- FindClusters(lamin_Monocytes, resolution = 0.3)
DimPlot(lamin_Monocytes, reduction = "umap", label = TRUE, repel = TRUE , raster=FALSE, label.size = 10) +  theme(text = element_text(size = 25)) + NoLegend()

DefaultAssay(lamin_Monocytes) <- "RNA"
VlnPlot(lamin_Monocytes, features = c("CD47"), pt.size = 2, raster=FALSE) + NoLegend()
lamin_Monocytes_CD47lo <- subset(lamin_Monocytes, subset = CD47<0.5)
lamin_Monocytes_CD47int <- subset(lamin_Monocytes, subset = CD47>0.5 & CD47<1.5)
lamin_Monocytes_CD47hi <- subset(lamin_Monocytes, subset = CD47>1.5)

anchors <- FindIntegrationAnchors(object.list = list(lamin_Monocytes_CD47lo, lamin_Monocytes_CD47int, lamin_Monocytes_CD47hi ), dims = 1:20)
Monocytes_combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = Monocytes_combined) <- "integrated"
saveRDS(Monocytes_combined, file = "C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 2/Monocytes/Monocytes_combined.rds")
Monocytes_combined <- readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 2/Monocytes/Monocytes_combined.rds")
Monocytes_combined <- ScaleData(object = Monocytes_combined, verbose = FALSE)
Monocytes_combined <- RunPCA(object = Monocytes_combined, npcs = 30, verbose = FALSE)
Monocytes_combined <- RunUMAP(object = Monocytes_combined, reduction = "pca", dims = 1:20) 
Monocytes_combined <- FindNeighbors(object = Monocytes_combined, reduction = "pca", dims = 1:20)
Monocytes_combined <- FindClusters(Monocytes_combined, resolution = 0.3)
DimPlot(Monocytes_combined, reduction = "umap", label = TRUE, repel = TRUE , raster=FALSE, label.size = 8)+  theme(text = element_text(size = 15)) + NoLegend()

Idents(Monocytes_combined) <- "CD47_expression"
DefaultAssay(Monocytes_combined) <- "RNA"

### pathway analysis for figure 2b
library(Signac)
library(monocle3)
library(Matrix)
library(patchwork)
set.seed(1234)
library(SCORPIUS)
library(patchwork)
library(tidyverse)
library(magrittr)
library(ReactomePA)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(biomaRt)
library(org.Hs.eg.db)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(monocle)
library(ggrepel)
library(EnhancedVolcano)
library(Nebulosa)
library(BiocFileCache)
library(RColorBrewer)
library(RColorBrewer)
library(Nebulosa)
library(fgsea)
library(msigdbr)
library(ReactomeGSA)
library(pathview)
geneSets_H = msigdbr(species = "Homo sapiens", category = "H")
head(geneSets_H)
geneSets = geneSets_H %>% split(x = .$gene_symbol, f = .$gs_name)
Monocytes_combined <-  AddModuleScore(Monocytes_combined, list(geneSets[[1]]), pool = NULL, nbin = 24, ctrl = 100, k = FALSE, assay = NULL, name = names(geneSets[1]), seed = 1, search = FALSE)
for (i in 1:50){
  Monocytes_combined <-  AddModuleScore(Monocytes_combined, list(geneSets[[i]]), pool = NULL, nbin = 24, ctrl = 100, k = FALSE, assay = NULL, name = names(geneSets[i]), seed = 1, search = FALSE)
  
}
head(Monocytes_combined[])
colnames(Monocytes_combined@meta.data)
score <- Monocytes_combined@meta.data[,c(40,47,51,63,64,65,66,67,77,80,84,85)]
ave <- aggregate(score,list(score$CD47_expression), mean)
rownames(ave) <- ave$Group.1
ave <- ave[,-c(2)]
ave <- ave[,-c(1)]
ave1 <- t(ave)
colfunc<-colorRampPalette(c("blue3","gray96","red3"))
heatmap(ave1,col=(colfunc(100)),Colv=NA, scale = "row",cluster_cols = FALSE, density.info = "none", trace = "none", border_color=NA, fontsize = 4, margins = c(5, 35))
VlnPlot(object = Monocytes_combined, "HALLMARK_INFLAMMATORY_RESPONSE1", raster=FALSE, pt.size = 0 ) + stat_summary(fun = median, geom='point', size = 15, colour = "black", shape = 45) + NoLegend() + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 25) ) + theme(axis.text.y = element_text(size = 25)) + theme(title = element_text(size = 20))
c("HALLMARK_INTERFERON_ALPHA_RESPONSE1", "HALLMARK_COMPLEMENT1", "HALLMARK_INFLAMMATORY_RESPONSE1" , "HALLMARK_INTERFERON_GAMMA_RESPONSE1" )


##### Figure 2b

monocytes <- readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 2/monocyte/Monocytes.rds")
monocytes <- FindVariableFeatures(monocytes, selection.method = "vst", nfeatures = 2000)
monocytes <- ScaleData(object = monocytes, verbose = FALSE)
monocytes <- RunPCA(object = monocytes, npcs = 30, verbose = FALSE)
monocytes <- RunUMAP(object = monocytes, reduction = "pca", dims = 1:20) 
monocytes <- FindNeighbors(object = monocytes, reduction = "pca", dims = 1:20)
monocytes <- FindClusters(monocytes, resolution = 0.5)
DimPlot(monocytes, reduction = "umap", label = TRUE, repel = TRUE , raster=FALSE, label.size = 8)+  theme(text = element_text(size = 15)) + NoLegend()


VlnPlot(monocytes, features = c("CD47"), pt.size = 1, raster=FALSE) + NoLegend()
monocytes_CD47lo <- subset(monocytes, subset = CD47<0.5)
monocytes_CD47int <- subset(monocytes, subset = CD47>0.5 & CD47<1.5)
monocytes_CD47hi <- subset(monocytes, subset = CD47>1.5)

#integration
monocytes_CD47lo$CD47_expression <- "CD47 Low"
monocytes_CD47lo <- NormalizeData(monocytes_CD47lo, normalization.method = "LogNormalize", scale.factor = 10000)
monocytes_CD47lo <- FindVariableFeatures(monocytes_CD47lo, selection.method = "vst", nfeatures = 2000)
monocytes_CD47lo <- ScaleData(object = monocytes_CD47lo, verbose = FALSE)
monocytes_CD47int$CD47_expression <- "CD47 Inter."
monocytes_CD47int <- NormalizeData(monocytes_CD47int, normalization.method = "LogNormalize", scale.factor = 10000)
monocytes_CD47int <- FindVariableFeatures(monocytes_CD47int, selection.method = "vst", nfeatures = 2000)
monocytes_CD47int <- ScaleData(object = monocytes_CD47int, verbose = FALSE)
monocytes_CD47hi$CD47_expression <- "CD47 High"
monocytes_CD47hi <- NormalizeData(monocytes_CD47hi, normalization.method = "LogNormalize", scale.factor = 10000)
monocytes_CD47hi <- FindVariableFeatures(monocytes_CD47hi, selection.method = "vst", nfeatures = 2000)
monocytes_CD47hi <- ScaleData(object = monocytes_CD47hi, verbose = FALSE)
anchors <- FindIntegrationAnchors(object.list = list(monocytes_CD47lo, monocytes_CD47int, monocytes_CD47hi ), dims = 1:20)
monocytes_combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = monocytes_combined) <- "integrated"
saveRDS(monocytes_combined, file = "C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 2/monocyte/monocytes_combined.rds")
monocytes_combined <- readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 2/monocyte/monocytes_combined.rds")
monocytes_combined <- ScaleData(object = monocytes_combined, verbose = FALSE)
monocytes_combined <- RunPCA(object = monocytes_combined, npcs = 30, verbose = FALSE)
monocytes_combined <- RunUMAP(object = monocytes_combined, reduction = "pca", dims = 1:20) 
monocytes_combined <- FindNeighbors(object = monocytes_combined, reduction = "pca", dims = 1:20)
monocytes_combined <- FindClusters(monocytes_combined, resolution = 0.5)
DimPlot(monocytes_combined, reduction = "umap", label = TRUE, repel = TRUE , raster=FALSE, label.size = 8)+  theme(text = element_text(size = 15)) + NoLegend()
Idents(monocytes_combined) <- "CD47_expression"
DefaultAssay(monocytes_combined) <- "RNA"
act_cytx <- c("CD47","TLR7","STAT1", "STAT2", "IRF7", "IRF9", "CD69", "ISG15", "LY6E", "IFI6","MX1", "IFIT2", "IFIT3", "IFIT5", "OASL", "CD86","GZMB","GZMA","GZMH", "GZMK",  "NKG7", "IL2RB","PRF1", "GNLY", "KLRB1", "PDCD1", "LAG3", "TIGIT")
ave <- AverageExpression(monocytes_combined)$RNA
gene <- intersect(rownames(ave), act_cytx)
ave_lamin <- ave[gene,]
colfunc<-colorRampPalette(c("blue3","gray96","red3"))
pheatmap(ave_lamin, col=(colfunc(100)), scale="row",border_color=NA, fontsize = 20,cluster_cols = FALSE)
DotPlot(monocytes_combined, features = act_cytx, cols = c("blue", "red", "green"), dot.scale = 20) + RotatedAxis() + theme(axis.text = element_text(size = 20)) + theme(axis.title.x = element_blank())+ theme(axis.title.y = element_blank())

DefaultAssay(monocytes) <- "RNA"

#IRF7
VlnPlot(monocytes, features = c("IRF7"), pt.size = 1, raster=FALSE) + NoLegend()
IRF7_low <- subset(monocytes, subset = IRF7<0.5)
IRF7_hi <- subset(monocytes, subset = IRF7>0.5)
IRF7_low$IRF7_expression <- "IRF7 low"
IRF7_hi$IRF7_expression <- "IRF7 high"
anchors <- FindIntegrationAnchors(object.list = list(IRF7_low, IRF7_hi), dims = 1:20)
IRF7_combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = IRF7_combined) <- "integrated"
IRF7_combined<-readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 2/monocyte/IRF7_combined.rds")
Idents(IRF7_combined) <- "IRF7_expression"
DefaultAssay(IRF7_combined) <- "RNA"
VlnPlot(IRF7_combined, features = c("CD47"), pt.size = 0, raster=FALSE)+NoLegend()+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) + theme(axis.title = element_text(size = 30))


#MX1
VlnPlot(monocytes, features = c("MX1"), pt.size = 1, raster=FALSE) + NoLegend()
MX1_low <- subset(monocytes, subset = MX1<0.5)
MX1_hi <- subset(monocytes, subset = MX1>0.5)
MX1_low$MX1_expression <- "MX1 low"
MX1_hi$MX1_expression <- "MX1 high"
anchors <- FindIntegrationAnchors(object.list = list(MX1_low, MX1_hi), dims = 1:20)
MX1_combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = MX1_combined) <- "integrated"
MX1_combined<-readRDS( "C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 2/monocyte/MX1_combined.rds")
Idents(MX1_combined) <- "MX1_expression"
DefaultAssay(MX1_combined) <- "RNA"
VlnPlot(MX1_combined, features = c("CD47"), pt.size = 0, raster=FALSE)+NoLegend()+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20))


#ISG15
VlnPlot(monocytes, features = c("ISG15"), pt.size = 1, raster=FALSE) + NoLegend()
ISG15_low <- subset(monocytes, subset = ISG15<0.5)
ISG15_hi <- subset(monocytes, subset = ISG15>0.5)
ISG15_low$ISG15_expression <- "ISG15 low"
ISG15_hi$ISG15_expression <- "ISG15 high"
anchors <- FindIntegrationAnchors(object.list = list(ISG15_low, ISG15_hi), dims = 1:20)
ISG15_combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = ISG15_combined) <- "integrated"
ISG15_combined<-readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 2/monocyte/ISG15_combined.rds")
Idents(ISG15_combined) <- "ISG15_expression"
DefaultAssay(ISG15_combined) <- "RNA"
VlnPlot(ISG15_combined, features = c("CD47"), pt.size = 0, raster=FALSE)+NoLegend()+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20))

#IL1B
VlnPlot(monocytes, features = c("IL1B"), pt.size = 1, raster=FALSE) + NoLegend()
IL1B_low <- subset(monocytes, subset = IL1B<0.5)
IL1B_hi <- subset(monocytes, subset = IL1B>0.5)
IL1B_low$IL1B_expression <- "IL1B low"
IL1B_hi$IL1B_expression <- "IL1B high"
anchors <- FindIntegrationAnchors(object.list = list(IL1B_low, IL1B_hi), dims = 1:20)
IL1B_combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = IL1B_combined) <- "integrated"
IL1B_combined<-readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 2/monocyte/IL1B_combined.rds")
Idents(IL1B_combined) <- "IL1B_expression"
DefaultAssay(IL1B_combined) <- "RNA"
VlnPlot(IL1B_combined, features = c("CD47"), pt.size = 0, raster=FALSE)+NoLegend()+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20))

#CD86
VlnPlot(monocytes, features = c("CD86"), pt.size = 1, raster=FALSE) + NoLegend()
CD86_low <- subset(monocytes, subset = CD86<0.5)
CD86_hi <- subset(monocytes, subset = CD86>0.5)
CD86_low$CD86_expression <- "CD86 low"
CD86_hi$CD86_expression <- "CD86 high"
anchors <- FindIntegrationAnchors(object.list = list(CD86_low, CD86_hi), dims = 1:20)
CD86_combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = CD86_combined) <- "integrated"
CD86_combined<-readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 2/monocyte/CD86_combined.rds")
Idents(CD86_combined) <- "CD86_expression"
DefaultAssay(CD86_combined) <- "RNA"
VlnPlot(CD86_combined, features = c("CD47"), pt.size = 0, raster=FALSE)+NoLegend()+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20))

#IFNAR1
VlnPlot(monocytes, features = c("IFNAR1"), pt.size = 1, raster=FALSE) + NoLegend()
IFNAR1_low <- subset(monocytes, subset = IFNAR1<0.5)
IFNAR1_hi <- subset(monocytes, subset = IFNAR1>0.5)
IFNAR1_low$IFNAR1_expression <- "IFNAR1 low"
IFNAR1_hi$IFNAR1_expression <- "IFNAR1 high"
anchors <- FindIntegrationAnchors(object.list = list(IFNAR1_low, IFNAR1_hi), dims = 1:20)
IFNAR1_combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = IFNAR1_combined) <- "integrated"
IFNAR1_combined<-readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 2/monocyte/IFNAR1_combined.rds")
Idents(IFNAR1_combined) <- "IFNAR1_expression"
DefaultAssay(IFNAR1_combined) <- "RNA"
VlnPlot(IFNAR1_combined, features = c("CD47"), pt.size = 0, raster=FALSE)+NoLegend()+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20))

#NKG7
VlnPlot(monocytes, features = c("NKG7"), pt.size = 1, raster=FALSE) + NoLegend()
NKG7_low <- subset(monocytes, subset = NKG7<0.5)
NKG7_hi <- subset(monocytes, subset = NKG7>0.5)
NKG7_low$NKG7_expression <- "NKG7 low"
NKG7_hi$NKG7_expression <- "NKG7 high"
anchors <- FindIntegrationAnchors(object.list = list(NKG7_low, NKG7_hi), dims = 1:20)
NKG7_combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = NKG7_combined) <- "integrated"
NKG7_combined<- readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 2/monocyte/NKG7_combined.rds")
Idents(NKG7_combined) <- "NKG7_expression"
DefaultAssay(NKG7_combined) <- "RNA"
VlnPlot(NKG7_combined, features = c("CD47"), pt.size = 0, raster=FALSE)+NoLegend()+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20))

#TNF
VlnPlot(monocytes, features = c("TNF"), pt.size = 1, raster=FALSE) + NoLegend()
TNF_low <- subset(monocytes, subset = TNF<0.5)
TNF_hi <- subset(monocytes, subset = TNF>0.5)
TNF_low$TNF_expression <- "TNF low"
TNF_hi$TNF_expression <- "TNF high"
anchors <- FindIntegrationAnchors(object.list = list(TNF_low, TNF_hi), dims = 1:20)
TNF_combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = TNF_combined) <- "integrated"
TNF_combined<- readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 2/monocyte/TNF_combined.rds")
Idents(TNF_combined) <- "TNF_expression"
DefaultAssay(TNF_combined) <- "RNA"
VlnPlot(TNF_combined, features = c("CD47"), pt.size = 0, raster=FALSE)+NoLegend()+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20))


####### Suppl. Figure 1 D,E
VlnPlot(monocytes, features = c("CD47"), pt.size = 1, raster=FALSE) + NoLegend()
CD47_low <- subset(monocytes, subset = CD47<0.5)
CD47_hi <- subset(monocytes, subset = CD47>0.5)
CD47_low$CD47_expression <- "CD47 low"
CD47_hi$CD47_expression <- "CD47 high"
anchors <- FindIntegrationAnchors(object.list = list(CD47_low, CD47_hi), dims = 1:20)
CD47_combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = CD47_combined) <- "integrated"
CD47_combined<-readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 2/monocyte/CD47_combined.rds")
Idents(CD47_combined) <- "CD47_expression"
DefaultAssay(CD47_combined) <- "RNA"
VlnPlot(CD47_combined, features = c("CD47"), pt.size = 0, raster=FALSE)+NoLegend()+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20))
CD47_combined <- ScaleData(object = CD47_combined, features = rownames(CD47_combined), verbose = FALSE)
CD47_hi_low <- FindMarkers(CD47_combined, ident.1 = "CD47 high", ident.2 = "CD47 low", min.pct = 0, logfc.threshold = 0)
saveRDS(CD47_hi_low, file = "C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 2/monocyte/CD47_hi_low.rds")
high_vs_loww <- FindMarkers(CD47_combined, ident.1 = "CD47 low", ident.2 = "CD47 high", min.pct = 0.1, logfc.threshold = 0.1)
EnhancedVolcano(high_vs_loww, lab = rownames(high_vs_loww), x = "avg_log2FC", y = "p_val_adj", xlab = bquote(~avg_Log[2]~ 'fold change'),  title = "Monocytes: CD47 high vs CD47 low cells", subtitle = bquote(italic("avg_log2FC")), labSize = 6, pointSize = 2, axisLabSize=30, titleLabSize=25, subtitleLabSize=20, captionLabSize=15, legendPosition = "none", pCutoff = 10e-10,  FCcutoff = 0.05, col=c('grey', 'grey', 'grey', 'red3'), colAlpha = 1, xlim = c(-1, 1),ylim = c(0,310))

EnhancedVolcano(CD47_hi_low, lab = rownames(CD47_hi_low), x = "avg_log2FC", y = "p_val_adj", xlab = bquote(~avg_Log[2]~ 'fold change'),  title = "Monocytes: CD47 high vs CD47 low cells", subtitle = bquote(italic("avg_log2FC")), labSize = 6, pointSize = 2, axisLabSize=30, titleLabSize=25, subtitleLabSize=20, captionLabSize=15, legendPosition = "none", pCutoff = 10e-100,  FCcutoff = 0.05, col=c('grey', 'grey', 'grey', 'red3'), colAlpha = 1, xlim = c(-1, 1),ylim = c(0,310))
ggsave("pic1.png", height = 12, width = 10 )
lamin <- rownames(CD47_hi_low)
go_enrich <- enrichGO(gene =lamin , OrgDb = "org.Hs.eg.db",  keyType = 'SYMBOL', readable = T, ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
dotplot(go_enrich, showCategory = 10, title = "GO: Monocytes: CD47 high vs CD47 low cells")
top100 <- slice_max(CD47_hi_low, n = 100, order_by = avg_log2FC)
lamin <- rownames(top100)

go_enrich@result
go_enrich@result %>% head()
write.csv(go_enrich,"GO.csv")
go_enrich@result %>%
  filter(Description== "response to virus" |
           Description=="positive regulation of leukocyte activation" | 
           Description=="positive regulation of cytokine production"|
           Description=="chemokine-mediated signaling pathway" |
           Description=="natural killer cell mediated cytotoxicity" |
           Description=="type I interferon production" |
           Description=="positive regulation of lymphocyte activation" |
           Description=="cellular response to type I interferon" |
           Description=="regulation of adaptive immune response" |
           Description=="regulation of T cell activation" |
           Description=="complement activation, classical pathway"  )%>%
  arrange(desc(p.adjust)) %>%
  ggplot(aes(y=factor(Description), x= Count, fill= p.adjust)) +
  geom_bar(stat = "identity", width = 0.7, size=20)  + 
  theme(axis.text.y = element_text(size=20),
        axis.title.y = element_blank()) +
  scale_fill_gradient (low= "blue", high= "red") +
  ggtitle( "Monocytes: CD47 high vs CD47 low cells") + theme(title = element_text(size = 25))
#ggsave("pic1.png", height = 6, width = 6 )


######### Figure 2: pDC analysis

lamin_pDC <- readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 2/Figure 2b/lamin_pDC.rds")
lamin_pDC <- FindVariableFeatures(lamin_pDC, selection.method = "vst", nfeatures = 2000)
lamin_pDC<- ScaleData(object = lamin_pDC, verbose = FALSE)
lamin_pDC<- RunPCA(object = lamin_pDC, npcs = 30, verbose = FALSE)
lamin_pDC<- RunUMAP(object = lamin_pDC, reduction = "pca", dims = 1:20) 
lamin_pDC<- FindNeighbors(object = lamin_pDC, reduction = "pca", dims = 1:20)
lamin_pDC<- FindClusters(lamin_pDC, resolution = 0.3)
DimPlot(lamin_pDC, reduction = "umap", label = TRUE, repel = TRUE , raster=FALSE, label.size = 10) +  theme(text = element_text(size = 25)) + NoLegend()
DefaultAssay(lamin_pDC) <- "RNA"
VlnPlot(lamin_pDC, features = c("CD47"), pt.size = 2, raster=FALSE) + NoLegend()
lamin_pDC_CD47lo <- subset(lamin_pDC, subset = CD47<0.5)
lamin_pDC_CD47int <- subset(lamin_pDC, subset = CD47>0.5 & CD47<1.6)
lamin_pDC_CD47hi <- subset(lamin_pDC, subset = CD47>1.6)

lamin_pDC_CD47lo$CD47_expression <- "CD47 Low"
lamin_pDC_CD47lo <- NormalizeData(lamin_pDC_CD47lo, normalization.method = "LogNormalize", scale.factor = 10000)
lamin_pDC_CD47lo <- FindVariableFeatures(lamin_pDC_CD47lo, selection.method = "vst", nfeatures = 2000)
lamin_pDC_CD47lo <- ScaleData(object = lamin_pDC_CD47lo, verbose = FALSE)
lamin_pDC_CD47int$CD47_expression <- "CD47 Inter."
lamin_pDC_CD47int <- NormalizeData(lamin_pDC_CD47int, normalization.method = "LogNormalize", scale.factor = 10000)
lamin_pDC_CD47int <- FindVariableFeatures(lamin_pDC_CD47int, selection.method = "vst", nfeatures = 2000)
lamin_pDC_CD47int <- ScaleData(object = lamin_pDC_CD47int, verbose = FALSE)
lamin_pDC_CD47hi$CD47_expression <- "CD47 High"
lamin_pDC_CD47hi <- NormalizeData(lamin_pDC_CD47hi, normalization.method = "LogNormalize", scale.factor = 10000)
lamin_pDC_CD47hi <- FindVariableFeatures(lamin_pDC_CD47hi, selection.method = "vst", nfeatures = 2000)
lamin_pDC_CD47hi <- ScaleData(object = lamin_pDC_CD47hi, verbose = FALSE)

anchors <- FindIntegrationAnchors(object.list = list(lamin_pDC_CD47lo, lamin_pDC_CD47int, lamin_pDC_CD47hi ), dims = 1:20)
pDC_combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = pDC_combined) <- "integrated"
saveRDS(pDC_combined, file = "C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 2/Figure 2b/pDC_combined.rds")
pDC_combined <- readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 2/Figure 2b/pDC_combined.rds")
pDC_combined <- ScaleData(object = pDC_combined, verbose = FALSE)
pDC_combined <- RunPCA(object = pDC_combined, npcs = 30, verbose = FALSE)
pDC_combined <- RunUMAP(object = pDC_combined, reduction = "pca", dims = 1:20) 
pDC_combined <- FindNeighbors(object = pDC_combined, reduction = "pca", dims = 1:20)
pDC_combined <- FindClusters(pDC_combined, resolution = 0.3)
DimPlot(pDC_combined, reduction = "umap", label = TRUE, repel = TRUE , raster=FALSE, label.size = 8)+  theme(text = element_text(size = 15)) + NoLegend()
Idents(pDC_combined) <- "CD47_expression"
DefaultAssay(pDC_combined) <- "RNA"
VlnPlot(pDC_combined, features = c("CD47"), raster=FALSE, pt.size = 0) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + NoLegend() + theme(axis.title.y = element_text(size = 25)) + theme(title = element_text(size = 20))

allcluster<- c("HLA-E", "HLA-DB5", "HLA-DB1", "HLA-DRA", "HLA-DQA1", "HLA-DPA1", "HLA-DMA", "CD5", "CD81", "BST2", "MT-CO1", "MT-CO3", "MT-CO3", "MT-ATP6", "NRP1", "S100A8", "MX1", "PRF1", "JCHAIN", "IRF8", "BCL11A", "DERL3", "GPR183", "PTCRA", "FAM129C", "SPIB", "IRF7", "MAP3K", "GZM", "IGJ", "AK128525", "SERPINF1", "ITM2C", "TCF4", "BCL11A", "LILRA4","TNFAIP2", "CLEC4C", "IL3RA", "MZB1", "RPL23", "PLEK", "RPS11", "SERBP1", "SET", "TPM3", "RPS17", "RPS20", "MTPN", "PTGDS", "FCER1G", "FTL", "RPS5", "NACA", "RPL19", "RPL29", "CLIC3", "PLCG2", "MT-ND1", "MT-ND5", "MT-ATP6", "MT-ND2", "MT-ND4", "MT-CO1", "MALAT1", "OGT", "ISG15", "LY6E", "IFI6", "ISG20", "IFI44L", "OAS1", "OAS3", "PARP12", "NKG7", "GNLY","PRF1","TNFSF10", "KLRB1", "GZMK")
all <- c( "TLR9", "TLR7", "TNF", "ZFP36", "SPARC", "TNFSF4", "PF4V1", "CXCL3", "CXCL5", "PPBP", "CCL5", "CCL3", "CCL2", "CCRL2", "CCL4", "CXCR2", "FFAR2", "TNFSF14", "IL1B", "IL1RN", "SERPINE1", "PDE4D", "TNFAIP3", "NLRP3", "TRIB1", "NFKBIA", "JUND", "JUN", "FOS", "F2R", "PTGS2", "EIF1", "DUSP1", "AREG", "BTG1", "KLF6", "LAPTM5",  "MYD88", "CD36", "LIPA", "RETN", "ANXA2","MAPK1", "TNFRSF21", "TNFAIP2", "IL10RB", "CD40", "CD38", "CD81", "HLA-DRA", "HLA-DMA", "HLA-C", "HLA-A", "HLA-E", "HLA-B", "HLA-F", "B2M", "USP8", "USP15", "USP24", "USP48", "USP3", "SLAMF7", "SOX4", "UFM1", "ARMCX3", "TPR", "AHI1", "UBTF", "LINC02812", "CLINT1", "TMEM123", "NEK8", "OXR1", "TMX3" )
lamin<- c("LILRA4","TNFAIP2", "CLEC4C", "IL3RA", "MZB1","TLR9", "TLR7","STAT1", "STAT2", "IRF7","CD47", "ISG15", "LY6E", "IFI6", "ISG20", "MX1", "OAS1", "OAS2", "OAS3", "PARP12", "NRP1", "S100A8", "IFIT2", "IFIT3", "IFIT5",  "RBBP6", "SHFL", "TRIM22", "TRIM25", "TRIM32", "TRIM69",  "IRF1","JCHAIN", "IRF8", "BCL11A", "DERL3", "GPR183", "PTCRA", "FAM129C", "SPIB", "IRF7", "MAP3K", "GZM", "IGJ", "AK128525", "SERPINF1", "ITM2C", "TCF4", "BCL11A", "IRF3", "IRF9", "CGAS", "USP18","SET", "TPM3", "RPS17", "RPS20", "MTPN", "PTGDS", "FCER1G", "SOCS1", "SOCS3", "IFITM2", "IFNL1", "IFNLR1", "IFNAR1", "IFNAR2", "IFNGR1", "IL10RB", "CD40", "CD38", "CD81", "HLA-DRA", "HLA-DMA", "HLA-C", "HLA-A", "HLA-E", "HLA-B", "HLA-F", "B2M", "USP8", "USP15", "USP24", "USP48", "USP3", "SLAMF7", "SOX4",  "ARMCX3", "TPR", "AHI1", "UBTF", "LINC02812", "CLINT1", "TMEM123", "NEK8", "OXR1", "TMX3")
all <- c( "NKG7", "GNLY","PRF1","TNFSF10", "KLRB1", "GZMK", "GZMB", "TNF", "CXCL3", "CCL5", "CCL3", "CCL2", "CCL4", "CXCR2", "TNFSF14", "IL1B", "IL1RN",  "TNFAIP3", "TNFAIP2", "IL10RB", "CD40", "CD38", "CD81", "HLA-DRA", "HLA-DMA", "HLA-C", "HLA-A", "HLA-E", "HLA-B", "HLA-F", "B2M",  "SLAMF7",  "IRF7", "ISG15", "LY6E", "IFI6", "ISG20", "MX1", "OAS1", "OAS2", "OAS3", "PARP12",  "IFIT2", "IFIT3", "IFIT5", "TRIM22", "TRIM25", "TRIM32", "TRIM69",  "IRF8" )
lamin<- c("TLR9", "TLR7","STAT1", "STAT2", "IRF7", "ISG15", "LY6E", "IFI6", "ISG20", "MX1", "OAS1", "OAS2", "OAS3", "PARP12",  "IFIT2", "IFIT3", "IFIT5",  "RBBP6", "SHFL", "TRIM22", "TRIM25", "TRIM32", "TRIM69",  "IRF1", "IRF3", "IRF9", "CGAS", "USP18", "SOCS1", "SOCS3", "IFITM2", "IFNL1", "IFNLR1", "IFNAR1", "IFNAR2", "IFNGR1", "IL10RB", "CD40", "CD38", "CD81", "HLA-DRA", "HLA-DMA", "HLA-C", "HLA-A", "HLA-E", "HLA-B", "HLA-F", "B2M", "USP8", "USP15", "USP24", "USP48", "USP3", "SLAMF7", "SOX4",  "ARMCX3", "TPR", "AHI1", "UBTF", "LINC02812", "CLINT1", "TMEM123", "NEK8", "OXR1", "TMX3")
ISG <- c( "STAT1", "STAT2", "ISG15", "LY6E", "IFI6", "ISG20", "MX1", "CD47", "OAS1", "OAS2", "OAS3", "PARP12",  "IFIT2", "IFIT3", "IFIT5",  "RBBP6", "SHFL", "TRIM22", "TRIM25", "TRIM32", "TRIM69",  "IRF1", "IRF3", "IRF9", "CGAS", "USP18", "SOCS1", "SOCS3", "IFITM2", "IFNL1", "IFNLR1", "IFNAR1", "IFNAR2", "IFNGR1")
proinflammation <- c("TLR9", "TNF", "ZFP36", "SPARC", "TNFSF4", "PF4V1",  "CXCL3", "CXCL5", "PPBP", "CCL5", "CCL3", "CXCL3", "CCL2", "CCL4", "CXCR2", "FFAR2", "TNFSF14", "IL1B", "IL1RN", "SERPINE1", "PDE4D", "IL6", "TNFAIP3", "NLRP3", "TRIB1", "NFKBIA", "JUND", "JUN", "FOS", "F2R", "PTGS2", "EIF1", "DUSP1", "AREG", "BTG1", "KLF6", "LAPTM5", "MYD88", "CD36", "LIPA", "RETN", "ANXA2", "TMEM176A", "TMEM176A")
all1 <- c("NKG7", "GNLY","PRF1","TNFSF10", "KLRB1", "GZMK", "GZMB", "TNF", "CXCL3", "CCL5", "CCL3", "CCL2","HLA-DRA", "HLA-DMA", "HLA-C", "HLA-A", "HLA-E", "HLA-B", "HLA-F", "B2M",  "SLAMF7", "ISG15", "LY6E", "IFI6", "ISG20", "MX1", "OAS1", "OAS2", "OAS3", "PARP12",  "IFIT2", "IFIT3", "IFIT5", "TRIM22", "TRIM25", "TRIM32", "TRIM69",  "IRF8")

IFNgenes <-c("TLR7","STAT1", "STAT2", "IRF7", "IRF9", "IRF8", "ISG15", "LY6E", "IFI6", "ISG20", "MX1", "OAS1", "OAS2", "OAS3",  "IFIT2", "IFIT3", "IFIT5",  "TRIM22", "TRIM25", "TRIM32", "TRIM69", "IFITM2", "IFNLR1", "IFNAR1", "IFNGR1")
antigenpretn <-c("HLA-C", "HLA-A", "HLA-E", "HLA-B", "HLA-F","HLA-DB5", "HLA-DB1", "HLA-DRA", "HLA-DQA1", "HLA-DPA1", "HLA-DMA", "MICA", "MICB", "LMP2",  "CIITA", "NLRC5" ,"CD40", "CD80", "B2M", " ICAM1", " ICAM2", " ICAM3")

ave <- AverageExpression(pDC_combined)$RNA
gene <- intersect(rownames(ave), antigenpretn)
ave_lamin <- ave[gene,]
colfunc<-colorRampPalette(c("blue3","gray96","red3"))
pheatmap(ave_lamin, col=(colfunc(100)), scale="row",border_color=NA, fontsize = 20,cluster_cols = FALSE)

### Gene ontology analysis

CD47hi_low <- FindMarkers(pDC_combined, ident.1 = "CD47 High", ident.2 = "CD47 Low", min.pct = 0.1, logfc.threshold = 0.1)
saveRDS(CD47hi_low, file = "C:/Users/au672897/Desktop/CD47 project/files/pDC/CD47hi_low.rds" )
top100_HIV_TLR9_pDC <- slice_max(CD47hi_low, n = 100, order_by = avg_log2FC)
top100gene_HIV_TLR9_pDC <- rownames(top100_HIV_TLR9_pDC)
top100gene_HIV_TLR9_pDC
go_enrich <- enrichGO(gene =top100gene_HIV_TLR9_pDC , OrgDb = "org.Hs.eg.db",  keyType = 'SYMBOL', readable = T, ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05)

GO_analysis <- enrichGO(gene, OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE, pool = TRUE)
GO_analysis

dotplot(go_enrich, showCategory = 20) + 
  theme(axis.text.y = element_text(size=20),
        axis.title.y = element_blank()) +
  scale_fill_gradient (low= "blue", high= "red") +
  ggtitle("pDC: CD47 High vs CD47 Low cells") + theme(title = element_text(size = 25))

write.csv(go_enrich,"GO.csv")
go_enrich@result
go_enrich@result %>% head()

go_enrich@result %>%
  filter(Description== "regulation of viral process" |
           Description=="defense response to virus" | 
           Description=="response to interferon-beta"|
           Description=="response to interferon-alpha"|
           Description== "type I interferon signaling pathway" |
           Description=="cytokine-mediated signaling pathway" |
           Description=="regulation of innate immune response" |
           Description=="defense response to symbiont" |
           Description=="positive regulation of innate immune response" |
           Description=="negative regulation of viral genome replication" |
           Description=="negative regulation of viral process" |
           Description=="defense response to symbiont" |
           Description=="complement activation, classical pathway" |
           Description=="cellular response to type I interferon" )%>%
  arrange(desc(p.adjust)) %>%
  ggplot(aes(y=factor(Description), x= Count, fill= p.adjust)) +
  geom_bar(stat = "identity", width = 0.7, size=20)  + 
  theme(axis.text.y = element_text(size=20),
        axis.title.y = element_blank()) +
  scale_fill_gradient (low= "blue", high= "red") +
  ggtitle("pDC: CD47 High vs CD47 Low cells") + theme(title = element_text(size = 25))


#### Figure 3 ##########
lamin_NKcells <- readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 3/lamin_NKcells.rds")
lamin_NKcells <- FindVariableFeatures(lamin_NKcells, selection.method = "vst", nfeatures = 2000)
lamin_NKcells<- ScaleData(object = lamin_NKcells, verbose = FALSE)
lamin_NKcells<- RunPCA(object = lamin_NKcells, npcs = 30, verbose = FALSE)
lamin_NKcells<- RunUMAP(object = lamin_NKcells, reduction = "pca", dims = 1:20) 
lamin_NKcells<- FindNeighbors(object = lamin_NKcells, reduction = "pca", dims = 1:20)
lamin_NKcells<- FindClusters(lamin_NKcells, resolution = 0.1)
DimPlot(lamin_NKcells, reduction = "umap", label = TRUE, repel = TRUE , raster=FALSE, label.size = 10) +  theme(text = element_text(size = 25)) + NoLegend()
lamin_NKcells<-RenameIdents(lamin_NKcells, "0"= "CD56dim NK cells", "1"= "CD56bright NK cells", "2"= "CD16neg NK cells", "3"= "NK/T cells")

DefaultAssay(lamin_NKcells) <- "RNA"
VlnPlot(lamin_NKcells, features = c("CD47"), pt.size = 0, raster=FALSE)+NoLegend()+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20))
lamin_NKcells_CD47lo <- subset(lamin_NKcells, subset = CD47<0.5)
lamin_NKcells_CD47int <- subset(lamin_NKcells, subset = CD47>0.5 & CD47<2.5)
lamin_NKcells_CD47hi <- subset(lamin_NKcells, subset = CD47>2.5)

lamin_NKcells_CD47lo$CD47_expression <- "CD47 Low"
lamin_NKcells_CD47lo <- NormalizeData(lamin_NKcells_CD47lo, normalization.method = "LogNormalize", scale.factor = 10000)
lamin_NKcells_CD47lo <- FindVariableFeatures(lamin_NKcells_CD47lo, selection.method = "vst", nfeatures = 2000)
lamin_NKcells_CD47lo <- ScaleData(object = lamin_NKcells_CD47lo, verbose = FALSE)
lamin_NKcells_CD47int$CD47_expression <- "CD47 Inter."
lamin_NKcells_CD47int <- NormalizeData(lamin_NKcells_CD47int, normalization.method = "LogNormalize", scale.factor = 10000)
lamin_NKcells_CD47int <- FindVariableFeatures(lamin_NKcells_CD47int, selection.method = "vst", nfeatures = 2000)
lamin_NKcells_CD47int <- ScaleData(object = lamin_NKcells_CD47int, verbose = FALSE)
lamin_NKcells_CD47hi$CD47_expression <- "CD47 High"
lamin_NKcells_CD47hi <- NormalizeData(lamin_NKcells_CD47hi, normalization.method = "LogNormalize", scale.factor = 10000)
lamin_NKcells_CD47hi <- FindVariableFeatures(lamin_NKcells_CD47hi, selection.method = "vst", nfeatures = 2000)
lamin_NKcells_CD47hi <- ScaleData(object = lamin_NKcells_CD47hi, verbose = FALSE)

anchors <- FindIntegrationAnchors(object.list = list(lamin_NKcells_CD47lo, lamin_NKcells_CD47int, lamin_NKcells_CD47hi ), dims = 1:20)
NK_combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = NK_combined) <- "integrated"
saveRDS(NK_combined, file = "C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 3/NK_combined.rds")
NK_combined <- readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 3/NK_combined.rds")
NK_combined <- ScaleData(object = NK_combined, verbose = FALSE)
NK_combined <- RunPCA(object = NK_combined, npcs = 30, verbose = FALSE)
NK_combined <- RunUMAP(object = NK_combined, reduction = "pca", dims = 1:20) 
NK_combined <- FindNeighbors(object = NK_combined, reduction = "pca", dims = 1:20)
NK_combined <- FindClusters(NK_combined, resolution = 0.3)
DimPlot(NK_combined, reduction = "umap", label = TRUE, repel = TRUE , raster=FALSE, label.size = 8)+  theme(text = element_text(size = 15)) + NoLegend()
Idents(NK_combined) <- "CD47_expression"
DefaultAssay(NK_combined) <- "RNA"
NKgenes <- c("KLRK1", "KLRC1", "NCR3", "KLRB1","ITGB2", "KLRD1", "FCGR3A","CD47","MX1","ISG15", "LY6E", "IFI6", "ISG20","TRIM22","IFIT3",  "NKG7", "GNLY","PRF1", "GZMB", "CCL5", "IFNG")
DotPlot(NK_combined, features = NKgenes, cols = c("blue", "red", "green"), dot.scale = 15, group.by= "CD47_expression") + RotatedAxis()+ theme(axis.text = element_text(size = 25)) + theme(axis.title.x = element_blank())+theme(axis.title.y = element_blank())

all <- c( "NKG7", "GNLY","PRF1","TNFSF10", "KLRB1", "GZMK", "GZMB", "TNF", "CXCL3", "CCL5", "CCL3", "CCL2", "CCL4", "CXCR2", "TNFSF14", "IL1B", "IL1RN",  "TNFAIP3", "TNFAIP2", "IL10RB", "CD40", "CD38", "CD81", "HLA-DRA", "HLA-DMA", "HLA-C", "HLA-A", "HLA-E", "HLA-B", "HLA-F", "B2M",  "SLAMF7",  "IRF7", "ISG15", "LY6E", "IFI6", "ISG20", "MX1", "OAS1", "OAS2", "OAS3", "PARP12",  "IFIT2", "IFIT3", "IFIT5", "TRIM22", "TRIM25", "TRIM32", "TRIM69",  "IRF8" )
lamin<- c("TLR9", "TLR7","STAT1", "STAT2", "IRF7", "ISG15", "LY6E", "IFI6", "ISG20", "MX1", "OAS1", "OAS2", "OAS3", "PARP12",  "IFIT2", "IFIT3", "IFIT5",  "RBBP6", "SHFL", "TRIM22", "TRIM25", "TRIM32", "TRIM69",  "IRF1", "IRF3", "IRF9", "CGAS", "USP18", "SOCS1", "SOCS3", "IFITM2", "IFNL1", "IFNLR1", "IFNAR1", "IFNAR2", "IFNGR1", "IL10RB", "CD40", "CD38", "CD81", "HLA-DRA", "HLA-DMA", "HLA-C", "HLA-A", "HLA-E", "HLA-B", "HLA-F", "B2M", "USP8", "USP15", "USP24", "USP48", "USP3", "SLAMF7", "SOX4",  "ARMCX3", "TPR", "AHI1", "UBTF", "LINC02812", "CLINT1", "TMEM123", "NEK8", "OXR1", "TMX3")
ISG <- c( "STAT1", "STAT2", "ISG15", "LY6E", "IFI6", "ISG20", "MX1", "CD47", "OAS1", "OAS2", "OAS3", "PARP12",  "IFIT2", "IFIT3", "IFIT5",  "RBBP6", "SHFL", "TRIM22", "TRIM25", "TRIM32", "TRIM69",  "IRF1", "IRF3", "IRF9", "CGAS", "USP18", "SOCS1", "SOCS3", "IFITM2", "IFNL1", "IFNLR1", "IFNAR1", "IFNAR2", "IFNGR1")
proinflammation <- c("TLR9", "TNF", "ZFP36", "SPARC", "TNFSF4", "PF4V1",  "CXCL3", "CXCL5", "PPBP", "CCL5", "CCL3", "CXCL3", "CCL2", "CCL4", "CXCR2", "FFAR2", "TNFSF14", "IL1B", "IL1RN", "SERPINE1", "PDE4D", "IL6", "TNFAIP3", "NLRP3", "TRIB1", "NFKBIA", "JUND", "JUN", "FOS", "F2R", "PTGS2", "EIF1", "DUSP1", "AREG", "BTG1", "KLF6", "LAPTM5", "MYD88", "CD36", "LIPA", "RETN", "ANXA2", "TMEM176A", "TMEM176A")
all1 <- c("NKG7", "GNLY","PRF1","TNFSF10", "KLRB1", "GZMK", "GZMB", "TNF", "CXCL3", "CCL5", "CCL3", "CCL2","HLA-DRA", "HLA-DMA", "HLA-C", "HLA-A", "HLA-E", "HLA-B", "HLA-F", "B2M",  "SLAMF7", "ISG15", "LY6E", "IFI6", "ISG20", "MX1", "OAS1", "OAS2", "OAS3", "PARP12",  "IFIT2", "IFIT3", "IFIT5", "TRIM22", "TRIM25", "TRIM32", "TRIM69",  "IRF8")
ave <- AverageExpression(NK_combined)$RNA
gene <- intersect(rownames(ave), all1)
ave_lamin <- ave[gene,]
colfunc<-colorRampPalette(c("blue3","gray96","red3"))
pheatmap(ave_lamin, col=(colfunc(100)), scale="row",border_color=NA, fontsize = 20,cluster_cols = FALSE)


#MX1
VlnPlot(lamin_NKcells, features = c("MX1"), pt.size = 1, raster=FALSE) + NoLegend()
MX1_low <- subset(lamin_NKcells, subset = MX1<0.5)
MX1_hi <- subset(lamin_NKcells, subset = MX1>0.5)
MX1_low$MX1_expression <- "MX1 low"
MX1_hi$MX1_expression <- "MX1 high"
anchors <- FindIntegrationAnchors(object.list = list(MX1_low, MX1_hi), dims = 1:20)
MX1_combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = MX1_combined) <- "integrated"
MX1_combined<-readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 3/MX1_combined.rds")
DefaultAssay(MX1_combined) <- "RNA"
Idents(MX1_combined) <- "MX1_expression"
DefaultAssay(MX1_combined) <- "RNA"
VlnPlot(MX1_combined, features = c("CD47"), pt.size = 0, raster=FALSE)+NoLegend()+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) + theme(axis.title = element_text(size = 30))

#ISG15
VlnPlot(lamin_NKcells, features = c("ISG15"), pt.size = 1, raster=FALSE) + NoLegend()
ISG15_low <- subset(lamin_NKcells, subset = ISG15<0.5)
ISG15_hi <- subset(lamin_NKcells, subset = ISG15>0.5)
ISG15_low$ISG15_expression <- "ISG15 low"
ISG15_hi$ISG15_expression <- "ISG15 high"
anchors <- FindIntegrationAnchors(object.list = list(ISG15_low, ISG15_hi), dims = 1:20)
ISG15_combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = ISG15_combined) <- "integrated"
ISG15_combined<-readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 3/ISG15_combined.rds")
DefaultAssay(ISG15_combined) <- "RNA"
Idents(ISG15_combined) <- "ISG15_expression"
DefaultAssay(ISG15_combined) <- "RNA"
VlnPlot(ISG15_combined, features = c("CD47"), pt.size = 0, raster=FALSE)+NoLegend()+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) + theme(axis.title = element_text(size = 30))

#NKG7
VlnPlot(lamin_NKcells, features = c("NKG7"), pt.size = 1, raster=FALSE) + NoLegend()
NKG7_low <- subset(lamin_NKcells, subset = NKG7<0.5)
NKG7_hi <- subset(lamin_NKcells, subset = NKG7>0.5)
NKG7_low$NKG7_expression <- "NKG7 low"
NKG7_hi$NKG7_expression <- "NKG7 high"
anchors <- FindIntegrationAnchors(object.list = list(NKG7_low, NKG7_hi), dims = 1:20)
NKG7_combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = NKG7_combined) <- "integrated"
NKG7_combined<- readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 3/NKG7_combined.rds")
DefaultAssay(NKG7_combined) <- "RNA"
Idents(NKG7_combined) <- "NKG7_expression"
DefaultAssay(NKG7_combined) <- "RNA"
VlnPlot(NKG7_combined, features = c("CD47"), pt.size = 0, raster=FALSE)+NoLegend()+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) + theme(axis.title = element_text(size = 30))


#PRF1
VlnPlot(lamin_NKcells, features = c("PRF1"), pt.size = 1, raster=FALSE) + NoLegend()
PRF1_low <- subset(lamin_NKcells, subset = PRF1<0.5)
PRF1_hi <- subset(lamin_NKcells, subset = PRF1>0.5)
PRF1_low$PRF1_expression <- "PRF1 low"
PRF1_hi$PRF1_expression <- "PRF1 high"
anchors <- FindIntegrationAnchors(object.list = list(PRF1_low, PRF1_hi), dims = 1:20)
PRF1_combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = PRF1_combined) <- "integrated"
PRF1_combined<- readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 3/PRF1_combined.rds")
DefaultAssay(PRF1_combined) <- "RNA"
Idents(PRF1_combined) <- "PRF1_expression"
DefaultAssay(PRF1_combined) <- "RNA"
VlnPlot(PRF1_combined, features = c("CD47"), pt.size = 0, raster=FALSE)+NoLegend()+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) + theme(axis.title = element_text(size = 30))

#GZMB
VlnPlot(lamin_NKcells, features = c("GZMB"), pt.size = 1, raster=FALSE) + NoLegend()
GZMB_low <- subset(lamin_NKcells, subset = GZMB<0.5)
GZMB_hi <- subset(lamin_NKcells, subset = GZMB>0.5)
GZMB_low$GZMB_expression <- "GZMB low"
GZMB_hi$GZMB_expression <- "GZMB high"
anchors <- FindIntegrationAnchors(object.list = list(GZMB_low, GZMB_hi), dims = 1:20)
GZMB_combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = GZMB_combined) <- "integrated"
GZMB_combined<- readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 3/GZMB_combined.rds")
DefaultAssay(GZMB_combined) <- "RNA"
Idents(GZMB_combined) <- "GZMB_expression"
DefaultAssay(GZMB_combined) <- "RNA"
VlnPlot(GZMB_combined, features = c("CD47"), pt.size = 0, raster=FALSE)+NoLegend()+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) + theme(axis.title = element_text(size = 30))

#GNLY
VlnPlot(lamin_NKcells, features = c("GNLY"), pt.size = 1, raster=FALSE) + NoLegend()
GNLY_low <- subset(lamin_NKcells, subset = GNLY<0.5)
GNLY_hi <- subset(lamin_NKcells, subset = GNLY>0.5)
GNLY_low$GNLY_expression <- "GNLY low"
GNLY_hi$GNLY_expression <- "GNLY high"
anchors <- FindIntegrationAnchors(object.list = list(GNLY_low, GNLY_hi), dims = 1:20)
GNLY_combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = GNLY_combined) <- "integrated"
GNLY_combined<- readRDS( "C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 3/GNLY_combined.rds")
DefaultAssay(GNLY_combined) <- "RNA"
Idents(GNLY_combined) <- "GNLY_expression"
DefaultAssay(GNLY_combined) <- "RNA"
VlnPlot(GNLY_combined, features = c("CD47"), pt.size = 0, raster=FALSE)+NoLegend()+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) + theme(axis.title = element_text(size = 30))


##### Figure 4 ########

library(Seurat)
library(SeuratObject)
library(dplyr)
library(pheatmap)
library(BiocManager)
library(multtest)
library(metap)
library(limma)
library(ggplot2)
library(cowplot)
lamin_CD8Tcells<- readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 4/lamin_CD8Tcells.rds")
lamin_CD8Tcells <- FindVariableFeatures(lamin_CD8Tcells, selection.method = "vst", nfeatures = 2000)
lamin_CD8Tcells<- ScaleData(object = lamin_CD8Tcells, verbose = FALSE)
lamin_CD8Tcells<- RunPCA(object = lamin_CD8Tcells, npcs = 30, verbose = FALSE)
lamin_CD8Tcells<- RunUMAP(object = lamin_CD8Tcells, reduction = "pca", dims = 1:20) 
lamin_CD8Tcells<- FindNeighbors(object = lamin_CD8Tcells, reduction = "pca", dims = 1:20)
lamin_CD8Tcells<- FindClusters(lamin_CD8Tcells, resolution = 0.2)
DimPlot(lamin_CD8Tcells, reduction = "umap", label = TRUE, repel = TRUE , raster=FALSE, label.size = 10) +  theme(text = element_text(size = 25)) + NoLegend()
DefaultAssay(lamin_CD8Tcells) <- "RNA"
new.cluster.ids <-c("Naive CD8+ T","Effector Mem. CD8+ T","Central Mem. CD8+ T","Effector CD8+ T","Effector CD8+ T")
names(new.cluster.ids) <- levels(lamin_CD8Tcells)
lamin_CD8Tcells <- RenameIdents(lamin_CD8Tcells, new.cluster.ids)
lamin_CD8Tcells <- RunTSNE(object = lamin_CD8Tcells)
DimPlot(object = lamin_CD8Tcells, reduction = "tsne", label = TRUE, repel = TRUE,  label.size = 10, pt.size = 2) + theme(text = element_text(size = 25)) + NoLegend()
VlnPlot(lamin_CD8Tcells, features = c("CD47"), pt.size = 0, raster=FALSE)+NoLegend()+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20))
RidgePlot(lamin_CD8Tcells, features = c( "CD47"), ncol = NULL) + NoLegend() + theme(title = element_text(size = 20)) + theme(axis.text.y = element_text(size = 25) ) + theme(axis.title.y = element_blank()) +theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20) )
RidgePlot(lamin_CD8Tcells, features = c( "PRF1"), ncol = NULL) + NoLegend() + theme(title = element_text(size = 20)) + theme(axis.text.y = element_text(size = 25) ) + theme(axis.title.y = element_blank()) +theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20) )

lamin_CD8Tcells_CD47lo <- subset(lamin_CD8Tcells, subset = CD47<0.5)
lamin_CD8Tcells_CD47int <- subset(lamin_CD8Tcells, subset = CD47>0.5 & CD47<2.5)
lamin_CD8Tcells_CD47hi <- subset(lamin_CD8Tcells, subset = CD47>2.5)


#integration
lamin_CD8Tcells_CD47lo$CD47_expression <- "CD47 Low"
lamin_CD8Tcells_CD47lo <- NormalizeData(lamin_CD8Tcells_CD47lo, normalization.method = "LogNormalize", scale.factor = 10000)
lamin_CD8Tcells_CD47lo <- FindVariableFeatures(lamin_CD8Tcells_CD47lo, selection.method = "vst", nfeatures = 2000)
lamin_CD8Tcells_CD47lo <- ScaleData(object = lamin_CD8Tcells_CD47lo, verbose = FALSE)
lamin_CD8Tcells_CD47int$CD47_expression <- "CD47 Inter."
lamin_CD8Tcells_CD47int <- NormalizeData(lamin_CD8Tcells_CD47int, normalization.method = "LogNormalize", scale.factor = 10000)
lamin_CD8Tcells_CD47int <- FindVariableFeatures(lamin_CD8Tcells_CD47int, selection.method = "vst", nfeatures = 2000)
lamin_CD8Tcells_CD47int <- ScaleData(object = lamin_CD8Tcells_CD47int, verbose = FALSE)
lamin_CD8Tcells_CD47hi$CD47_expression <- "CD47 High"
lamin_CD8Tcells_CD47hi <- NormalizeData(lamin_CD8Tcells_CD47hi, normalization.method = "LogNormalize", scale.factor = 10000)
lamin_CD8Tcells_CD47hi <- FindVariableFeatures(lamin_CD8Tcells_CD47hi, selection.method = "vst", nfeatures = 2000)
lamin_CD8Tcells_CD47hi <- ScaleData(object = lamin_CD8Tcells_CD47hi, verbose = FALSE)

anchors <- FindIntegrationAnchors(object.list = list(lamin_CD8Tcells_CD47lo, lamin_CD8Tcells_CD47int, lamin_CD8Tcells_CD47hi ), dims = 1:20)
CD8T_combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = CD8T_combined) <- "integrated"

CD8T_combined <- readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 4/CD8T_combined.rds")
CD8T_combined <- ScaleData(object = CD8T_combined, verbose = FALSE)
CD8T_combined <- RunPCA(object = CD8T_combined, npcs = 30, verbose = FALSE)
CD8T_combined <- RunUMAP(object = CD8T_combined, reduction = "pca", dims = 1:20) 
CD8T_combined <- FindNeighbors(object = CD8T_combined, reduction = "pca", dims = 1:20)
CD8T_combined <- FindClusters(CD8T_combined, resolution = 0.3)
DimPlot(CD8T_combined, reduction = "umap", label = TRUE, repel = TRUE , raster=FALSE, label.size = 8)+  theme(text = element_text(size = 15)) + NoLegend()

Idents(CD8T_combined) <- "CD47_expression"
DefaultAssay(CD8T_combined) <- "RNA"
cytotoxic_exhaustion <-c("IFNG","GZMB", "CCR7",  "NKG7", "IL2RB","PTPRC",  "PRF1", "HMGB1", "KLRB1",  "PDCD1", "LAG3", "TNFSF14", "KLRG1" )
act_cytx <- c("CD47", "CD69","TNFRSF9", "CD28", "ISG15","ISG20",  "OAS1", "LY6E", "IFI6","MX1", "OAS1", "OASL", "IFNG","GZMB","GZMA","GZMH", "GZMK",  "NKG7", "IL2RB","PRF1", "ISG15", "KLRB1", "PRF1", "PDCD1", "LAG3", "TIGIT")
ave <- AverageExpression(CD8T_combined)$RNA
gene <- intersect(rownames(ave), cytotoxic_exhaustion)
ave_lamin <- ave[gene,]
colfunc<-colorRampPalette(c("blue3","gray96","red3"))
pheatmap(ave_lamin, col=(colfunc(100)), scale="row",border_color=NA, fontsize = 20,cluster_cols = FALSE)
RidgePlot(CD8T_combined, features = c( "CD47"), ncol = NULL) + NoLegend() + theme(title = element_text(size = 20)) + theme(axis.text.y = element_text(size = 25) ) + theme(axis.title.y = element_blank()) +theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20) )
ggsave("pic1.png", height = 9, width = 4 )
VlnPlot(CD8T_combined, features = c("CD47"), raster=FALSE, pt.size = 0) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + NoLegend() + theme(axis.title.y = element_text(size = 25)) + theme(title = element_text(size = 20))
ggsave("pic1.png", height = 9, width = 4 )

#NKG7
VlnPlot(lamin_CD8Tcells, features = c("NKG7"), pt.size = 1, raster=FALSE) + NoLegend()
NKG7_low <- subset(lamin_CD8Tcells, subset = NKG7<0.5)
NKG7_hi <- subset(lamin_CD8Tcells, subset = NKG7>0.5)
NKG7_low$NKG7_expression <- "NKG7 low"
NKG7_hi$NKG7_expression <- "NKG7 high"
anchors <- FindIntegrationAnchors(object.list = list(NKG7_low, NKG7_hi), dims = 1:20)
NKG7_combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = NKG7_combined) <- "integrated"
DefaultAssay(NKG7_combined) <- "RNA"
Idents(NKG7_combined) <- "NKG7_expression"
DefaultAssay(NKG7_combined) <- "RNA"
VlnPlot(NKG7_combined, features = c("CD47"), pt.size = 0, raster=FALSE)+NoLegend()+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) + theme(axis.title = element_text(size = 30))

#PRF1
VlnPlot(lamin_CD8Tcells, features = c("PRF1"), pt.size = 1, raster=FALSE) + NoLegend()
PRF1_low <- subset(lamin_CD8Tcells, subset = PRF1<0.5)
PRF1_hi <- subset(lamin_CD8Tcells, subset = PRF1>0.5)
PRF1_low$PRF1_expression <- "PRF1 low"
PRF1_hi$PRF1_expression <- "PRF1 high"
anchors <- FindIntegrationAnchors(object.list = list(PRF1_low, PRF1_hi), dims = 1:20)
PRF1_combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = PRF1_combined) <- "integrated"
DefaultAssay(PRF1_combined) <- "RNA"
Idents(PRF1_combined) <- "PRF1_expression"
DefaultAssay(PRF1_combined) <- "RNA"
VlnPlot(PRF1_combined, features = c("CD47"), pt.size = 0, raster=FALSE)+NoLegend()+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) + theme(axis.title = element_text(size = 30))

#GZMB
VlnPlot(lamin_CD8Tcells, features = c("GZMB"), pt.size = 1, raster=FALSE) + NoLegend()
GZMB_low <- subset(lamin_CD8Tcells, subset = GZMB<0.5)
GZMB_hi <- subset(lamin_CD8Tcells, subset = GZMB>0.5)
GZMB_low$GZMB_expression <- "GZMB low"
GZMB_hi$GZMB_expression <- "GZMB high"
anchors <- FindIntegrationAnchors(object.list = list(GZMB_low, GZMB_hi), dims = 1:20)
GZMB_combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = GZMB_combined) <- "integrated"
DefaultAssay(GZMB_combined) <- "RNA"
Idents(GZMB_combined) <- "GZMB_expression"
DefaultAssay(GZMB_combined) <- "RNA"
VlnPlot(GZMB_combined, features = c("CD47"), pt.size = 0, raster=FALSE)+NoLegend()+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) + theme(axis.title = element_text(size = 30))

#ISG15
VlnPlot(lamin_CD8Tcells, features = c("ISG15"), pt.size = 1, raster=FALSE) + NoLegend()
ISG15_low <- subset(lamin_CD8Tcells, subset = ISG15<0.5)
ISG15_hi <- subset(lamin_CD8Tcells, subset = ISG15>0.5)
ISG15_low$ISG15_expression <- "ISG15 low"
ISG15_hi$ISG15_expression <- "ISG15 high"
anchors <- FindIntegrationAnchors(object.list = list(ISG15_low, ISG15_hi), dims = 1:20)
ISG15_combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = ISG15_combined) <- "integrated"
DefaultAssay(ISG15_combined) <- "RNA"
Idents(ISG15_combined) <- "ISG15_expression"
DefaultAssay(ISG15_combined) <- "RNA"
VlnPlot(ISG15_combined, features = c("CD47"), pt.size = 0, raster=FALSE)+NoLegend()+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) + theme(axis.title = element_text(size = 30))

### suppl.3 CD4 ########

CD4T<- readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 4/CD4 T cells/CD4T.rds")
CD4T <- FindVariableFeatures(CD4T, selection.method = "vst", nfeatures = 2000)
CD4T<- ScaleData(object = CD4T, verbose = FALSE)
CD4T<- RunPCA(object = CD4T, npcs = 30, verbose = FALSE)
CD4T<- RunUMAP(object = CD4T, reduction = "pca", dims = 1:20) 
CD4T<- FindNeighbors(object = CD4T, reduction = "pca", dims = 1:20)
CD4T<- FindClusters(CD4T, resolution = 0.3)
DimPlot(CD4T, reduction = "umap", label = TRUE, repel = TRUE , raster=FALSE, label.size = 10) +  theme(text = element_text(size = 25)) + NoLegend()
DefaultAssay(CD4T) <- "RNA"
CD4T <- RunTSNE(object = CD4T, check_duplicates = FALSE)
DimPlot(object = CD4T, reduction = "tsne", label = TRUE, repel = TRUE,  label.size = 10, pt.size = 2) + theme(text = element_text(size = 25)) + NoLegend()
new.cluster.ids <-c("Naive CD4+ T","Effector CD4+ T","Central Mem. CD4+ T","Effector Mem. CD4+ T", "Naive CD4+ T", "Cytotoxic CD4+ T", "Cytotoxic CD4+ T")
names(new.cluster.ids) <- levels(CD4T)
CD4T <- RenameIdents(CD4T, new.cluster.ids)
CD4T <- RunTSNE(object = CD4T)
DimPlot(object = CD4T, reduction = "tsne", label = TRUE, repel = TRUE,  label.size = 10, pt.size = 2) + theme(text = element_text(size = 25)) + NoLegend()
ggsave("pic1.png", height = 7, width = 10)
VlnPlot(CD4T, features = c("CD47"), pt.size = 0, raster=FALSE)+NoLegend()+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20))
RidgePlot(CD4T, features = c( "CD3D"), ncol = NULL) + NoLegend() + theme(title = element_text(size = 20)) + theme(axis.text.y = element_text(size = 25) ) + theme(axis.title.y = element_blank()) +theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20) )

###### Figure 5 ##########

library(Seurat)
library(SeuratObject)
library(dplyr)
library(pheatmap)
library(BiocManager)
library(multtest)
library(metap)
library(limma)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(Signac)
library(monocle3)
library(Matrix)
library(patchwork)
B_cells <- readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 5/B_cells.rds")
B_cells <- ScaleData(object = B_cells, verbose = FALSE)
B_cells <- RunPCA(object = B_cells, npcs = 30, verbose = FALSE)
B_cells <- RunUMAP(object = B_cells, reduction = "pca", dims = 1:20) 
B_cells <- FindNeighbors(object = B_cells, reduction = "pca", dims = 1:20)
B_cells <- FindClusters(B_cells, resolution = 0.3)
DimPlot(B_cells, reduction = "umap", label = TRUE, repel = TRUE , raster=FALSE, label.size = 8)+  theme(text = element_text(size = 15)) + NoLegend()
new.cluster.ids <-c("IgM Mem. B cells","Classical Mem. B cells","Naive B cells","Transitional B cells","Classical Mem. B cells", "Double neg. B cells", "Naive B cells")
names(new.cluster.ids) <- levels(B_cells)
B_cells <- RenameIdents(B_cells, new.cluster.ids)
#B_cells <- RunTSNE(object = B_cells)
DimPlot(object = B_cells, reduction = "tsne", label = TRUE, repel = TRUE,  label.size = 10, pt.size = 2) + theme(text = element_text(size = 25)) + NoLegend()
DefaultAssay(B_cells) <- "RNA"
VlnPlot(B_cells, features = c("CD47"), pt.size = 0, raster=FALSE) + NoLegend() + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20))
RidgePlot(B_cells, features = c( "IGHM"), ncol = NULL) + NoLegend() + theme(title = element_text(size = 20)) + theme(axis.text.y = element_text(size = 25) ) + theme(axis.title.y = element_blank()) +theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size = 20) )
DefaultAssay(B_cells) <- "RNA"
VlnPlot(B_cells, features = c("CD47"), pt.size = 2, raster=FALSE) + NoLegend()
B_cells_CD47lo <- subset(B_cells, subset = CD47<0.5)
B_cells_CD47int <- subset(B_cells, subset = CD47>0.5 & CD47<2)
B_cells_CD47hi <- subset(B_cells, subset = CD47>2)
B_cells_CD47lo$CD47_expression <- "CD47 Low"
B_cells_CD47lo <- NormalizeData(B_cells_CD47lo, normalization.method = "LogNormalize", scale.factor = 10000)
B_cells_CD47lo <- FindVariableFeatures(B_cells_CD47lo, selection.method = "vst", nfeatures = 2000)
B_cells_CD47lo <- ScaleData(object = B_cells_CD47lo, verbose = FALSE)
B_cells_CD47int$CD47_expression <- "CD47 Inter."
B_cells_CD47int <- NormalizeData(B_cells_CD47int, normalization.method = "LogNormalize", scale.factor = 10000)
B_cells_CD47int <- FindVariableFeatures(B_cells_CD47int, selection.method = "vst", nfeatures = 2000)
B_cells_CD47int <- ScaleData(object = B_cells_CD47int, verbose = FALSE)
B_cells_CD47hi$CD47_expression <- "CD47 High"
B_cells_CD47hi <- NormalizeData(B_cells_CD47hi, normalization.method = "LogNormalize", scale.factor = 10000)
B_cells_CD47hi <- FindVariableFeatures(B_cells_CD47hi, selection.method = "vst", nfeatures = 2000)
B_cells_CD47hi <- ScaleData(object = B_cells_CD47hi, verbose = FALSE)

anchors <- FindIntegrationAnchors(object.list = list(B_cells_CD47lo, B_cells_CD47int, B_cells_CD47hi ), dims = 1:20)
B_combinedb <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = B_combinedb) <- "integrated"
#B_combined <- readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 5/Figure 5b/B_combinedb.rds")
B_combinedb <- ScaleData(object = B_combinedb, verbose = FALSE)
B_combinedb <- RunPCA(object = B_combinedb, npcs = 30, verbose = FALSE)
B_combinedb <- RunUMAP(object = B_combinedb, reduction = "pca", dims = 1:20) 
B_combinedb <- FindNeighbors(object = B_combinedb, reduction = "pca", dims = 1:20)
B_combinedb <- FindClusters(B_combinedb, resolution = 0.3)
DimPlot(B_combinedb, reduction = "umap", label = TRUE, repel = TRUE , raster=FALSE, label.size = 8)+  theme(text = element_text(size = 15)) + NoLegend()

Idents(B_combinedb) <- "CD47_expression"
DefaultAssay(B_combinedb) <- "RNA"
Bcell_genes <- c("MS4A1", "TNFRSF17","CD274",  "PIGS", "DERL3", "CD19", "CD22", "TCL1A", "CD83", "CD37", "CD79B", "MZB1",   "CD40", "CD80",  "IGHM", "IGHD", "TLR9","STAT1", "STAT2", "IRF7", "ISG15", "LY6E", "IFI6", "ISG20", "MX1", "OAS1", "OAS2", "HLA-DRA", "HLA-DMA", "HLA-C", "HLA-E", "HLA-B", "B2M",  "SLAMF7")
DotPlot(B_combinedb, features = Bcell_genes, cols = c("blue", "red", "green"), dot.scale = 15, group.by= "CD47_expression") + RotatedAxis()+ theme(axis.text = element_text(size = 25)) + theme(axis.title.x = element_blank())+theme(axis.title.y = element_blank())
ave <- AverageExpression(B_combinedb)$RNA
gene <- intersect(rownames(ave), Bcell_genes)
ave_lamin <- ave[gene,]
colfunc<-colorRampPalette(c("blue3","gray96","red3"))
pheatmap(ave_lamin, col=(colfunc(100)), scale="row",border_color=NA, fontsize = 20,cluster_cols = FALSE)

lamin_Bcells <- readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 5/Figure 5b/lamin_Bcells.rds")
lamin_Bcells <- FindVariableFeatures(lamin_Bcells, selection.method = "vst", nfeatures = 2000)
lamin_Bcells<- ScaleData(object = lamin_Bcells, verbose = FALSE)
lamin_Bcells<- RunPCA(object = lamin_Bcells, npcs = 30, verbose = FALSE)
lamin_Bcells<- RunUMAP(object = lamin_Bcells, reduction = "pca", dims = 1:20) 
lamin_Bcells<- FindNeighbors(object = lamin_Bcells, reduction = "pca", dims = 1:20)
lamin_Bcells<- FindClusters(lamin_Bcells, resolution = 0.3)
DimPlot(lamin_Bcells, reduction = "umap", label = TRUE, repel = TRUE , raster=FALSE, label.size = 10) +  theme(text = element_text(size = 25)) + NoLegend()
DefaultAssay(lamin_Bcells) <- "RNA"
VlnPlot(lamin_Bcells, features = c("CD47"), pt.size = 2, raster=FALSE) + NoLegend()
lamin_Bcells_CD47lo <- subset(lamin_Bcells, subset = CD47<0.5)
lamin_Bcells_CD47int <- subset(lamin_Bcells, subset = CD47>0.5 & CD47<2)
lamin_Bcells_CD47hi <- subset(lamin_Bcells, subset = CD47>2)

#integration
lamin_Bcells_CD47lo$CD47_expression <- "CD47 Low"
lamin_Bcells_CD47lo <- NormalizeData(lamin_Bcells_CD47lo, normalization.method = "LogNormalize", scale.factor = 10000)
lamin_Bcells_CD47lo <- FindVariableFeatures(lamin_Bcells_CD47lo, selection.method = "vst", nfeatures = 2000)
lamin_Bcells_CD47lo <- ScaleData(object = lamin_Bcells_CD47lo, verbose = FALSE)
lamin_Bcells_CD47int$CD47_expression <- "CD47 Inter."
lamin_Bcells_CD47int <- NormalizeData(lamin_Bcells_CD47int, normalization.method = "LogNormalize", scale.factor = 10000)
lamin_Bcells_CD47int <- FindVariableFeatures(lamin_Bcells_CD47int, selection.method = "vst", nfeatures = 2000)
lamin_Bcells_CD47int <- ScaleData(object = lamin_Bcells_CD47int, verbose = FALSE)
lamin_Bcells_CD47hi$CD47_expression <- "CD47 High"
lamin_Bcells_CD47hi <- NormalizeData(lamin_Bcells_CD47hi, normalization.method = "LogNormalize", scale.factor = 10000)
lamin_Bcells_CD47hi <- FindVariableFeatures(lamin_Bcells_CD47hi, selection.method = "vst", nfeatures = 2000)
lamin_Bcells_CD47hi <- ScaleData(object = lamin_Bcells_CD47hi, verbose = FALSE)

anchors <- FindIntegrationAnchors(object.list = list(lamin_Bcells_CD47lo, lamin_Bcells_CD47int, lamin_Bcells_CD47hi ), dims = 1:20)
B_combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = B_combined) <- "integrated"
B_combined <- readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 5/Figure 5b/B_combined.rds")
B_combined <- ScaleData(object = B_combined, verbose = FALSE)
B_combined <- RunPCA(object = B_combined, npcs = 30, verbose = FALSE)
B_combined <- RunUMAP(object = B_combined, reduction = "pca", dims = 1:20) 
B_combined <- FindNeighbors(object = B_combined, reduction = "pca", dims = 1:20)
B_combined <- FindClusters(B_combined, resolution = 0.3)
DimPlot(B_combined, reduction = "umap", label = TRUE, repel = TRUE , raster=FALSE, label.size = 8)+  theme(text = element_text(size = 15)) + NoLegend()

Idents(B_combined) <- "CD47_expression"
DefaultAssay(B_combined) <- "RNA"
B_combined.markers <- FindAllMarkers(B_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
DotPlot(B_combined, features = Bcell_genes, cols = c("blue", "red", "green"), dot.scale = 15, group.by= "CD47_expression") + RotatedAxis()+ theme(axis.text = element_text(size = 25)) + theme(axis.title.x = element_blank())+theme(axis.title.y = element_blank())

Bcell_genes <- c("MS4A1", "TNFRSF17","CD274",  "PIGS", "DERL3", "CD19", "CD22", "TCL1A", "CD83", "CD37", "CD79B", "MZB1",   "CD40", "CD80",  "IGHM", "IGHD", "TLR9","STAT1", "STAT2", "IRF7", "ISG15", "LY6E", "IFI6", "ISG20", "MX1", "OAS1", "OAS2", "HLA-DRA", "HLA-DMA", "HLA-C", "HLA-E", "HLA-B", "B2M",  "SLAMF7")
ave <- AverageExpression(B_combined)$RNA
gene <- intersect(rownames(ave), Bcell_genes)
ave_lamin <- ave[gene,]
colfunc<-colorRampPalette(c("blue3","gray96","red3"))
pheatmap(ave_lamin, col=(colfunc(100)), scale="row",border_color=NA, fontsize = 20,cluster_cols = FALSE)

#### DEGs and ontology #####

CD47lo_CD47hi <- FindMarkers(B_combined, ident.1 = "CD47 High", ident.2 = "CD47 Low", min.pct = 0, logfc.threshold = 0)
saveRDS(CD47lo_CD47hi, file = "C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 5/Figure 5b/CD47lo_CD47hi.rds" )
CD47lo_CD47hi <- readRDS("C:/Users/au672897/Desktop/CD47 project/Figure files/Figure 5/Figure 5b/CD47lo_CD47hi.rds")
EnhancedVolcano(CD47lo_CD47hi, lab = rownames(CD47lo_CD47hi), x = "avg_log2FC", y = "p_val_adj", xlab = bquote(~avg_Log[2]~ 'fold change'),  title = "B cells CD47 high vs CD47 low", subtitle = bquote(italic("avg_log2FC")), labSize = 7, pointSize = 2, axisLabSize=30, titleLabSize=25, subtitleLabSize=20, captionLabSize=15, legendPosition = "none", pCutoff = 10e-10,  FCcutoff = 0.05, col=c('grey', 'grey', 'grey', 'red3'), colAlpha = 0.9, xlim = c(-1, 1),ylim = c(0,310))

top100_CD47lo_CD47hi <- slice_max(CD47lo_CD47hi, n = 100, order_by = avg_log2FC)
top100gene_CD47lo_CD47hi <- rownames(top100_CD47lo_CD47hi)
top100gene_CD47lo_CD47hi
go_enrich <- enrichGO(gene =top100gene_CD47lo_CD47hi , OrgDb = "org.Hs.eg.db",  keyType = 'SYMBOL', readable = T, ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
barplot(go_enrich, showCategory = 20, color = "p.adjust", title = "CD47 High ") 

write.csv(go_enrich,"GO.csv")
go_enrich@result
go_enrich@result %>% head()

go_enrich@result %>%
  filter(Description== "cellular response to type I interferon" |
           Description=="B cell receptor signaling pathway" | 
           Description=="positive regulation of B cell activation"|
           Description== "antigen receptor-mediated signaling pathway" |
           Description=="humoral immune response" |
           Description=="immunoglobulin mediated immune response" |
           Description=="antimicrobial humoral response" |
           Description=="regulation of viral life cycle" |
           Description=="complement activation, classical pathway" |
           Description=="cellular response to type I interferon" |
           Description=="lymphocyte mediated immunity"|
           Description=="positive regulation of cell activation" |
           Description=="cellular response to type I interferon" |
           Description=="lymphocyte mediated immunity"|
           Description=="response to interferon-beta"|
           Description=="response to interferon-beta"|
           Description=="activation of immune response"|
           Description=="response to chemokine"|
           Description=="cellular response to interferon-gamma" |
           Description=="programmed necrotic cell death" |
           Description==" activation of immune response " |
           Description=="chemokine-mediated signaling pathway"   )%>%
  arrange(desc(p.adjust)) %>%
  ggplot(aes(y=factor(Description), x= Count, fill= p.adjust)) +
  geom_bar(stat = "identity", width = 0.7, size=20)  + 
  theme(axis.text.y = element_text(size=20),
        axis.title.y = element_blank()) +
  scale_fill_gradient (low= "blue", high= "red") +
  ggtitle("B cells CD47 high vs CD47 low") + theme(title = element_text(size = 25))


