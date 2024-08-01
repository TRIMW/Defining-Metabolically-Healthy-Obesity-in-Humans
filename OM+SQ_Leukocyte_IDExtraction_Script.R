MS##OM Immune Identification and Sectioning
OM_Immune <- subset(x = OM_merge, idents = c("T cells", "NK cells", "Myeloid cells", "B cells"))
OM_Immune <- ScaleData(OM_Immune, features = VariableFeatures(OM_Immune))
OM_Immune <- RunPCA(OM_Immune, features = VariableFeatures(OM_Immune), npcs = 100)
OM_Immune <- RunUMAP(OM_Immune, reduction = "pca", dims = 1:50)
OM_Immune <- FindNeighbors(OM_Immune, dims = 1:10)
OM_Immune <- FindClusters(OM_Immune, resolution = 0.5)
OM_Immune <- RunUMAP(OM_Immune, dims = 1:10)
DimPlot(OM_Immune)
DimPlot(OM_Immune, split.by = "Sample") + theme(aspect.ratio = 1)
OM_Immune <- subset(x = OM_Immune, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"))

OM_Immune_IDs <- c("CD4 T", "Mac", "CD4 T", "CD8 T", "Granulo", "GD T", "NK Cell", "B Cell", "NKT-like", "Mono", "?", "?")

OM_Immune_markers <- FindAllMarkers(OM_Immune, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
OM_Immune_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

OM_Immune_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> OM_Immune_top10
DoHeatmap(OM_Immune, features = OM_Immune_top10$gene) + NoLegend()


##SQ Immune Identification and Sectioning
SQ_Immune <- subset(x = SQ_merge, idents = c("3", "4", "5", "9", "10", "11"))
SQ_Immune <- ScaleData(SQ_Immune, features = VariableFeatures(SQ_Immune))
SQ_Immune <- RunPCA(SQ_Immune, features = VariableFeatures(SQ_Immune), npcs = 100)
SQ_Immune <- RunUMAP(SQ_Immune, reduction = "pca", dims = 1:50)
ElbowPlot(SQ_Immune, ndims=50)

SQ_Immune <- FindNeighbors(SQ_Immune, dims = 1:12)
SQ_Immune <- FindClusters(SQ_Immune, resolution = 0.5)
SQ_Immune <- RunUMAP(SQ_Immune, dims = 1:12)

SQ_Immune_markers <- FindAllMarkers(SQ_Immune, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SQ_Immune_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

SQ_Immune_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> OM_Immune_top10
DoHeatmap(SQ_Immune, features = OM_Immune_top10$gene) + NoLegend()