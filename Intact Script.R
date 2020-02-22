library(Seurat)
library(dplyr)
library(ggplot2)
DefaultAssay(Intact2) <- "integrated"

#Set Current idents
Idents(object = CtrlvARKOvIntact.combined) <- "seurat_clusters"

#Setup workspace to make file calling & saving easy

#setwd("W:C:\Users\cbaker\Documents\CB")


#Intact

Intact.data <- Read10X(data.dir = "./filtered_feature_bc_matrix")
Intact <- CreateSeuratObject(counts = Intact.data, project = "IntactSC", min.cells = 3, min.features = 200)

Intact[["percent.mt"]] <- PercentageFeatureSet(Intact, pattern = "^mt-")
VlnPlot(Intact, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(Intact, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Intact, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
Intact <- subset(Intact, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 15)
VlnPlot(Intact, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Intact <- NormalizeData(Intact, normalization.method = "LogNormalize", scale.factor = 10000)
Intact <- NormalizeData(Intact)
Intact <- FindVariableFeatures(Intact, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(Intact), 10)
plot3 <- VariableFeaturePlot(Intact)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
CombinePlots(plots = list(plot3, plot4))
all.genes <- rownames(Intact)
Intact <- ScaleData(Intact, features = all.genes)
Intact <- RunPCA(Intact, features = VariableFeatures(object = Intact))
print(Intact[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Intact, dims = 1:2, reduction = "pca")
DimPlot(Intact, reduction = "pca")
ElbowPlot(Intact, ndims = 30)
Intact <- JackStraw(Intact, num.replicate = 100)
Intact <- ScoreJackStraw(Intact, dims = 1:20)
JackStrawPlot(Intact, dims = 1:14)
Intact <- FindNeighbors(Intact, dims = 1:24)
Intact <- FindClusters(Intact, resolution = 0.5)
head(Idents(Intact), 5)


Intact.markers <- FindAllMarkers(Intact, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Intact.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- Intact.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
Intact <- RunTSNE(Intact, dims = 1:24)
DimPlot(Intact, reduction = "tsne", label = TRUE)
FeaturePlot(Intact, features = c("Krt8", "Vim", "Myh11", "Fbln", "Krt5", "Pecam", "Tyrobp", "EGFP", "Ar"))
RidgePlot(Intact, features = c("Ar", "EGFP"), ncol = 2)
FeaturePlot(Intact, features = c("Ar", "EGFP"), blend = TRUE)
FeaturePlot(Intact, pt.size = 1.5, features = c("Ar", "EGFP"), blend = TRUE)
FeaturePlot(Intact, pt.size = 1, features = c("Ar", "EGFP"), blend = TRUE)
top5 <- Intact.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Intact, features = top5$gene) + NoLegend()
plotA <- FeatureScatter(Intact, feature1 = "Vim", feature2 = "Acta2")
plotB <- DimPlot(Intact)
CombinePlots(plots = list(plotA, plotB))
write.csv(Intact.markers, file = "Intactmarkers.csv")
DoHeatmap(Intact, features = IntactTop10$gene) + NoLegend() + scale_fill_gradientn(colors = c("blue", "white", "red"))

FeaturePlot(Intact, features = c("EGFP", "Ar", "Gli1", "MYC-transgene", "Myc", "Epcam", "Vim", "Krt8", "Krt5", "Acta2", "Fbln1", "Pecam1"), cols = c("lightgrey", "red"), min.cutoff = 0)
DotPlot(Intact, features = c("Ccl5", "Cd3d",	"Cd3e",	"Cd28",	"Cxcr6", "Aqp1", "Cldn5", "Pecam1", "Cdh5", "Plvap",	"Tyrobp", "Spi1", "C1qc", "C1qa", "C1qb", "Actg2", "Myh11", "Tagln", "Acta2", "Tpm2", "Lum", "Pdgfra", "Ptgs2", "Apod", "Fbln1", "Aqp3", "Col17a1", "Krt5", "Krt15", "Krt14", "Stard10", "Krt19", "Krt18", "Krt8", "Cldn3", "Ar", "Gli1", "EGFP", "MYC-transgene", "Myc"), cols = c("light grey", "red"), col.min = 0) + RotatedAxis()
row.names(Intact)

IntactandCtrlMyc.combined <- RenameIdents(object = IntactandCtrlMyc.combined, '0' = "0", '1' = "E", '2' = "E", '3' = "3", '4' = "E", '5' = "E", '6' = "6", '7' = "7", '8' = "8", '9' = "E", '10' = "E", '11' = "E", '12' = "12", '13' = "13")
Epi2 <- subset(x = IntactandCtrlMyc.combined, idents = "E")

DefaultAssay(IntactandCtrlMyc.combined) <- "RNA"
Idents(object = IntactandCtrlMyc.combined) <- "seurat_clusters"
Idents(object = IntactandCtrlMyc.combined) <- "MycExp"
MycPosIntactandCtrlMyc.combined <- subset(x=IntactandCtrlMyc.combined, subset = `MYC-transgene` > 0.01)
MycNegIntactandCtrlMyc.combined <- subset(x=IntactandCtrlMyc.combined, subset = `MYC-transgene` < 0.01)
Idents(object = MycPosIntactandCtrlMyc.combined) <- "Mycpos"
Idents(object = MycNegIntactandCtrlMyc.combined) <- "Mycneg"
MycPosIntactandCtrlMyc.combined[["MycExp"]] <- Idents(object = MycPosIntactandCtrlMyc.combined)
MycNegIntactandCtrlMyc.combined[["MycExp"]] <- Idents(object = MycNegIntactandCtrlMyc.combined)
IntactandCtrlMyc.combinedMyc <- merge(x = MycPosIntactandCtrlMyc.combined, y = MycNegIntactandCtrlMyc.combined)
IntactandCtrlMyc.combined[["MycExp"]] <- Idents(object = IntactandCtrlMyc.combinedMyc)
DimPlot(IntactandCtrlMyc.combined, reduction = "tsne", group.by = "MycExp")
IntactandCtrlMyc.combinedMyc.markers <- FindAllMarkers(IntactandCtrlMyc.combined, logfc.threshold = 0)
IntactandCtrlMyc.combinedMyc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
IntactandCtrlMyc.combinedMyctop10 <- IntactandCtrlMyc.combinedMyc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(IntactandCtrlMyc.combined, features = IntactandCtrlMyc.combinedMyctop10$gene) + NoLegend() + scale_fill_gradientn(colors = c("blue", "white", "red"))
write.csv(IntactandCtrlMyc.combinedMyc.markers, file = "IntactandCtrlMyc.combinedMyc.csv")

#Merge
#Stash old idents
MycCtrl[["orig.clusters"]] <- Idents(object = MycCtrl)
Intact[["orig.clusters"]] <- Idents(object = Intact)

#Set Current idents


MycCtrl$stim <- "Control"
Intact$stim <- "Intact"

IntactandCtrlMyc.anchors <- FindIntegrationAnchors(object.list = list(MycCtrl, Intact), dims = 1:20)

IntactandCtrlMyc.combined <- IntegrateData(anchorset = IntactandCtrlMyc.anchors, dims = 1:20)

DefaultAssay(IntactandCtrlMyc.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
IntactandCtrlMyc.combined <- ScaleData(IntactandCtrlMyc.combined, verbose = FALSE)
IntactandCtrlMyc.combined <- RunPCA(IntactandCtrlMyc.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
IntactandCtrlMyc.combined <- FindNeighbors(IntactandCtrlMyc.combined, reduction = "pca", dims = 1:20)
IntactandCtrlMyc.combined <- FindClusters(IntactandCtrlMyc.combined, resolution = 0.5)
IntactandCtrlMyc.combined <- RunTSNE(IntactandCtrlMyc.combined, reduction = "pca", dims = 1:20)
DimPlot(IntactandCtrlMyc.combined, reduction = "tsne", label = TRUE)
DimPlot(IntactandCtrlMyc.combined, reduction = "tsne", group.by = "stim")

DefaultAssay(IntactandCtrlMyc.combined) <- "RNA"
all.genes <- rownames(IntactandCtrlMyc.combined)
IntactandCtrlMyc.combined <- ScaleData(IntactandCtrlMyc.combined, features = all.genes)

ARKOvCtrlMyc.combined.markers <- FindAllMarkers(ARKOvCtrlMyc.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ARKOvCtrlMyc.combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(ARKOvCtrlMyc.combined.markers, file = "ARKOvCtrlMyc.combined.csv")

DotPlot(IntactandCtrlMyc.combined, features = c("Ccl5", "Cd3d",	"Cd3e",	"Cd28",	"Cxcr6", "Aqp1", "Cldn5", "Pecam1", "Cdh5", "Plvap",	"Tyrobp", "Spi1", "C1qc", "C1qa", "C1qb", "Actg2", "Myh11", "Tagln", "Acta2", "Tpm2", "Lum", "Pdgfra", "Ptgs2", "Apod", "Fbln1", "Aqp3", "Col17a1", "Krt5", "Krt15", "Krt14", "Stard10", "Krt19", "Krt18", "Krt8", "Cldn3", "Ar", "Gli1", "EGFP", "MYC-transgene", "Myc"), cols = c("light grey", "red"), col.min = 0) + RotatedAxis()
row.names(IntactandCtrlMyc.combined)
Idents(object = IntactandCtrlMyc.combined) <- "stim"
IntactandCtrlMyc.combined$stim.seurat_clusters <- paste(Idents(IntactandCtrlMyc.combined), IntactandCtrlMyc.combined$seurat_clusters, sep = "_")
IntactandCtrlMyc.combined$stim <- Idents(IntactandCtrlMyc.combined)

Idents(object = IntactandCtrlMyc.combined) <- "stim.seurat_clusters.MycExp"
IntactandCtrlMyc.combined$stim.seurat_clusters.MycExp <- paste(Idents(IntactandCtrlMyc.combined), IntactandCtrlMyc.combined$MycExp, sep = "_")
IntactandCtrlMyc.combined$stim.seurat_clusters <- Idents(IntactandCtrlMyc.combined)


table(IntactandCtrlMyc.combined$stim.seurat_clusters.MycExp)

Epi <- RenameIdents(object = Epi, '0' = "0", '1' = "E", '2' = "E", '3' = "3", '4' = "E", '5' = "E", '6' = "6", '7' = "7", '8' = "Normal", '9' = "E", '10' = "E", '11' = "E")
NormalLum <- subset(x = Epi, idents = "Normal")

IntactandCtrlMyc.combined <- RenameIdents(object = IntactandCtrlMyc.combined, 'Intact_1_Mycpos' = "HGPIN")

IntactMyc1 <- subset(x = IntactandCtrlMyc.combined, idents = "HGPIN")

#Merge
#Stash old idents
MycARKO[["orig.clusters"]] <- Idents(object = MycARKO)
Intact[["orig.clusters"]] <- Idents(object = Intact)

#Set Current idents

Idents(object = IntactandARKOMyc.combined) <- "seurat_clusters"

MycCtrl$stim <- "ARKO"
Intact$stim <- "Intact"

IntactandARKOMyc.anchors <- FindIntegrationAnchors(object.list = list(MycARKO, Intact), dims = 1:20)

IntactandARKOMyc.combined <- IntegrateData(anchorset = IntactandARKOMyc.anchors, dims = 1:20)

DefaultAssay(IntactandARKOMyc.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
IntactandARKOMyc.combined <- ScaleData(IntactandARKOMyc.combined, verbose = FALSE)
IntactandARKOMyc.combined <- RunPCA(IntactandARKOMyc.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
IntactandARKOMyc.combined <- FindNeighbors(IntactandARKOMyc.combined, reduction = "pca", dims = 1:20)
IntactandARKOMyc.combined <- FindClusters(IntactandARKOMyc.combined, resolution = 0.5)
IntactandARKOMyc.combined <- RunTSNE(IntactandARKOMyc.combined, reduction = "pca", dims = 1:20)
DimPlot(IntactandARKOMyc.combined, reduction = "tsne", label = TRUE)
DimPlot(IntactandARKOMyc.combined, reduction = "tsne", group.by = "stim")

DefaultAssay(IntactandARKOMyc.combined) <- "RNA"
all.genes <- rownames(IntactandARKOMyc.combined)
IntactandARKOMyc.combined <- ScaleData(IntactandARKOMyc.combined, features = all.genes)
                     
IntactandARKOMyc.combined <- RenameIdents(object = IntactandARKOMyc.combined, '0' = "Fb", '1' = "Fb", '2' = "Fb", '3' = "L", '4' = "E", '5' = "L", '6' = "6", '7' = "7", '8' = "Fb", '9' = "9", '10' = "L", '11' = "11", '12' = "L", '13' = "13", '14' = "14", '15' = "Fb", '16' = "16")
Fib2 <- subset(x = IntactandARKOMyc.combined, idents = "Fb")

DefaultAssay(Lum3) <- "RNA"
Idents(object = Lum3) <- "seurat_clusters"
Lum3[["percent.mt"]] <- PercentageFeatureSet(Lum3, pattern = "^mt-")
VlnPlot(Lum3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(Lum3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Lum3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
Lum3 <- NormalizeData(Lum3, normalization.method = "LogNormalize", scale.factor = 10000)
Lum3 <- NormalizeData(Lum3)
Lum3 <- FindVariableFeatures(Lum3, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(Lum3), 10)
plot3 <- VariableFeaturePlot(Lum3)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
CombinePlots(plots = list(plot3, plot4))
all.genes <- rownames(Lum3)
Lum3 <- ScaleData(Lum3, features = all.genes)
Lum3 <- RunPCA(Lum3, features = VariableFeatures(object = Lum3))
print(Lum3[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Lum3, dims = 1:2, reduction = "pca")
DimPlot(Lum3, reduction = "pca")
ElbowPlot(Lum3)
Lum3 <- JackStraw(Lum3, num.replicate = 100)
Lum3 <- ScoreJackStraw(Lum3, dims = 1:20)
JackStrawPlot(intact, dims = 1:14)
intact <- FindNeighbors(intact, dims = 1:20)
intact <- FindClusters(intact, resolution = 0.5)
head(Idents(intact), 5)

#Analysis of expression profiles with whole data set

intact.markers <- FindAllMarkers(intact, min.pct = 0.25, logfc.threshold = 0.25)
intact.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

DotPlot(Lum3, features = c("Ccl5", "Cd3d",	"Cd3e",	"Cd28",	"Cxcr6", "Aqp1", "Cldn5", "Pecam1", "Cdh5", "Plvap",	"Tyrobp", "Spi1", "C1qc", "C1qa", "C1qb", "Actg2", "Myh11", "Tagln", "Acta2", "Tpm2", "Lum", "Pdgfra", "Ptgs2", "Apod", "Fbln1", "Aqp3", "Col17a1", "Krt5", "Krt15", "Krt14", "Stard10", "Krt19", "Krt18", "Krt8", "Cldn3", "Ar", "Gli1", "EGFP", "MYC-transgene", "Myc"), cols = c("light grey", "red"), col.min = 0) + RotatedAxis()

#Idents(object = Lum3) <- "MycExp"
#MycPosLum3 <- subset(x=Lum3, subset = `MYC-transgene` > 0.01)
#MycNegLum3 <- subset(x=Lum3, subset = `MYC-transgene` < 0.01)
#Idents(object = MycPosLum3) <- "Mycpos"
#Idents(object = MycNegLum3) <- "Mycneg"
#MycPosLum3[["MycExp"]] <- Idents(object = MycPosLum3)
#MycNegLum3[["MycExp"]] <- Idents(object = MycNegLum3)
#Lum3Myc <- merge(x = MycPosLum3, y = MycNegLum3)
Lum3[["MycExp"]] <- Idents(object = Lum3Myc)
DotPlot(Lum3, features = c("Ccl5", "Cd3d",	"Cd3e",	"Cd28",	"Cxcr6", "Aqp1", "Cldn5", "Pecam1", "Cdh5", "Plvap",	"Tyrobp", "Spi1", "C1qc", "C1qa", "C1qb", "Actg2", "Myh11", "Tagln", "Acta2", "Tpm2", "Lum", "Pdgfra", "Ptgs2", "Apod", "Fbln1", "Aqp3", "Col17a1", "Krt5", "Krt15", "Krt14", "Stard10", "Krt19", "Krt18", "Krt8", "Cldn3", "Ar", "Gli1", "EGFP", "MYC-transgene", "Myc"), cols = c("light grey", "red"), col.min = 0) + RotatedAxis()

DimPlot(Lum3, reduction = "tsne", group.by = "MycExp")
Lum2Myc.markers <- FindAllMarkers(Lum2, logfc.threshold = 0)
Lum2Myc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
Lum2Myctop10 <- Lum2Myc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(Lum2, features = Lum2Myctop10$gene) + NoLegend() + scale_fill_gradientn(colors = c("blue", "white", "red"))
write.csv(Lum2Myc.markers, file = "Lum2Myc.csv")

Idents(object = Lum3) <- "stim.seurat_clusters"
Lum2$stim.MycExp <- paste(Idents(Lum2), Lum2$MycExp, sep = "_")
Lum2$stim <- Idents(Lum2)

IntactMycposvsARKOMycPos.markers <- FindMarkers(Lum2, ident.1 = "Intact_Mycpos", ident.2 = "ARKO_Mycpos", logfc.threshold = 0.1, min.pct = 0.1)
write.csv(IntactMycposvsARKOMycPos.markers, file = "IntactMycposvsARKOMycPosmarkers.csv")

table(NormalLum$seurat_clusters)


IntactPINvNormal <- merge(x = NormalLum, y = IntactMyc1)

IntactPINvNormal.markers <- FindAllMarkers(IntactPINvNormal, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(IntactPINvNormal.markers, file = "IntactPINvNormalmarkers.csv")

Idents(object = CtrlvARKOvIntact.combined) <- "stim.seurat_clusters"
VlnPlot(Lum2, features = c("Timp1", "Timp2"), pt.size = 0, idents = c("Intact_Mycpos", "ARKO_Mycpos"))


#Stash old idents
MycCtrl[["orig.clusters"]] <- Idents(object = MycCtrl)
MycARKO[["orig.clusters"]] <- Idents(object = MycARKO)
Intact[["orig.clusters"]] <- Idents(object = Intact)


#Set Current idents


MycCtrl$stim <- "Control"
MycARKO$stim <- "ARKO"
Intact$stim <- "Intact"


CtrlvARKOvIntact.anchors <- FindIntegrationAnchors(object.list = list(MycCtrl, MycARKO, Intact), dims = 1:20)

CtrlvARKOvIntact.combined <- IntegrateData(anchorset = CtrlvARKOvIntact.anchors, dims = 1:20)

DefaultAssay(CtrlvARKOvIntact.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
CtrlvARKOvIntact.combined <- ScaleData(CtrlvARKOvIntact.combined, verbose = FALSE)
CtrlvARKOvIntact.combined <- RunPCA(CtrlvARKOvIntact.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
CtrlvARKOvIntact.combined <- FindNeighbors(CtrlvARKOvIntact.combined, reduction = "pca", dims = 1:20)
CtrlvARKOvIntact.combined <- FindClusters(CtrlvARKOvIntact.combined, resolution = 0.5)
CtrlvARKOvIntact.combined <- RunTSNE(CtrlvARKOvIntact.combined, reduction = "pca", dims = 1:20)
DimPlot(CtrlvARKOvIntact.combined, reduction = "tsne", label = TRUE)
DimPlot(CtrlvARKOvIntact.combined, reduction = "tsne", group.by = "stim")

DefaultAssay(CtrlvARKOvIntact.combined) <- "RNA"
all.genes <- rownames(CtrlvARKOvIntact.combined)
CtrlvARKOvIntact.combined <- ScaleData(CtrlvARKOvIntact.combined, features = all.genes)


DotPlot(CtrlvARKOvIntact.combined, features = c("Ccl5", "Cd3d",	"Cd3e",	"Cd28",	"Cxcr6", "Aqp1", "Cldn5", "Pecam1", "Cdh5", "Plvap",	"Tyrobp", "Spi1", "C1qc", "C1qa", "C1qb", "Actg2", "Myh11", "Tagln", "Acta2", "Tpm2", "Lum", "Pdgfra", "Ptgs2", "Apod", "Fbln1", "Aqp3", "Col17a1", "Krt5", "Krt15", "Krt14", "Stard10", "Krt19", "Krt18", "Krt8", "Cldn3", "Ar", "Gli1", "EGFP"), cols = c("light grey", "red")) + RotatedAxis()
FeaturePlot(CtrlvARKOvIntact.combined, features = c("EGFP", "Ar", "Gli1", "Epcam", "Vim", "Krt8", "Krt5", "Acta2", "Fbln1"), cols = c("lightgrey", "red"), min.cutoff = 0)

Idents(object = Fib2) <- "seurat_clusters"
CtrlvARKOvIntact.combined$stim.seurat_clusters <- paste(Idents(CtrlvARKOvIntact.combined), CtrlvARKOvIntact.combined$seurat_clusters, sep = "_")

Idents(object = CtrlvARKOvIntact.combined) <- "stim"
Lum3$stim.MycExp <- paste(Idents(Lum3), Lum3$MycExp, sep = "_")

table(CtrlvARKOvIntact.combined$stim.seurat_clusters)

CtrlvARKOvIntact.combined <- RenameIdents(object = CtrlvARKOvIntact.combined, '0' = "Fb", '1' = "Fb", '2' = "L", '3' = "3", '4' = "Fb", '5' = "5", '6' = "L", '7' = "L", '8' = "8", '9' = "9", '10' = "L", '11' = "L", '12' = "L", '13' = "L", '14' = "14", '15' = "15", '16' = "16", '17' = "17")
Lum3 <- subset(x = CtrlvARKOvIntact.combined, idents = "L")

Fib2.markers <- FindAllMarkers(Fib2, logfc.threshold = 0.1, min.pct = 0.1)
write.csv(Fib2.markers, file = "IntactvARKOFb.csv")


Clusters2and6v7and10.markers <- FindMarkers(CtrlvARKOvIntact.combined, ident.1 = "HG", ident.2 = "LG", logfc.threshold = 0.1, min.pct = 0.1)
write.csv(Clusters2and6v7and10.markers, file = "Clusters2and6v7and10markers.csv")

VlnPlot(Lum3, features = c("Timp1"), pt.size = 0, idents = c("Control_Mycpos", "ARKO_Mycpos", "Intact_Mycpos"))


DefaultAssay(CtrlvARKOvIntact.combined) <- "RNA"
Idents(object = CtrlvARKOvIntact.combined) <- "stim"
Idents(object = CtrlvARKOvIntact.combined) <- "MycExp"
MycPosCtrlvARKOvIntact.combined <- subset(x=CtrlvARKOvIntact.combined, subset = `MYC-transgene` > 0.01)
MycNegCtrlvARKOvIntact.combined <- subset(x=CtrlvARKOvIntact.combined, subset = `MYC-transgene` < 0.01)
Idents(object = MycPosCtrlvARKOvIntact.combined) <- "Mycpos"
Idents(object = MycNegCtrlvARKOvIntact.combined) <- "Mycneg"
MycPosCtrlvARKOvIntact.combined[["MycExp"]] <- Idents(object = MycPosCtrlvARKOvIntact.combined)
MycNegCtrlvARKOvIntact.combined[["MycExp"]] <- Idents(object = MycNegCtrlvARKOvIntact.combined)
CtrlvARKOvIntact.combinedMyc <- merge(x = MycPosCtrlvARKOvIntact.combined, y = MycNegCtrlvARKOvIntact.combined)
CtrlvARKOvIntact.combined[["MycExp"]] <- Idents(object = CtrlvARKOvIntact.combinedMyc)
DimPlot(CtrlvARKOvIntact.combined, reduction = "tsne", group.by = "MycExp")
CtrlvARKOvIntact.combinedMyc.markers <- FindAllMarkers(CtrlvARKOvIntact.combined, logfc.threshold = 0)
CtrlvARKOvIntact.combinedMyc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
CtrlvARKOvIntact.combinedMyctop10 <- CtrlvARKOvIntact.combinedMyc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(CtrlvARKOvIntact.combined, features = CtrlvARKOvIntact.combinedMyctop10$gene) + NoLegend() + scale_fill_gradientn(colors = c("blue", "white", "red"))
write.csv(CtrlvARKOvIntact.combinedMyc.markers, file = "CtrlvARKOvIntact.combinedMyc.csv")


