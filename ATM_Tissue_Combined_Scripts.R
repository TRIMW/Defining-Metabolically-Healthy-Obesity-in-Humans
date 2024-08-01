
##Combine extrated OM macrophage populations into one seurat object
merged_OM_ATMs <- merge(LAM, PVM, add.cell.ids = c("LAM", "PVM"))
merged_SQ_ATMs <- merge(SQ_LAM, SQ_PVM, add.cell.ids = c("LAM", "PVM"))

VlnPlot(merged_ATMs, features = c('TNF', 'MRC1', 'LYVE1', 'TREM2', 'CD9', 'S100A9'), split.by = 'Sample', pt.size = 0, cols = c('OM_MHO' = 'gray', 'OM_MUO' = 'gold'), combine = TRUE, stack = TRUE, flip = TRUE) +
  NoLegend() + theme(axis.text.x = element_text())
VlnPlot(merged_SQ_ATMs, features = c('TNF', 'MRC1', 'LYVE1', 'TREM2', 'CD9', 'S100A9'), split.by = 'Sample', pt.size = 0, cols = c('OM_MHO' = 'gray', 'OM_MUO' = 'gold'), combine = TRUE, stack = TRUE, flip = TRUE) +
  NoLegend() + theme(axis.text.x = element_text())


Idents(merged_OM_ATMs) <- "Sample"

merged_OM_MHO_ATMs <- subset(merged_OM_ATMs, idents = c("OM_MHO"))
Idents(merged_OM_ATMs) <- "seurat_clusters"

merged_OM_MUO_ATMs <- subset(merged_OM_ATMs, idents = c("OM_MUO"))
Idents(merged_OM_ATMs) <- "seurat_clusters"

Macrophages <- FindVariableFeatures(Macrophages, selection.method = 'vst', nfeatures = 3000)
Macrophages <- ScaleData(Macrophages, features = VariableFeatures(Macrophages))
Macrophages <- RunPCA(Macrophages, features = VariableFeatures(Macrophages), npcs = 100)
Macrophages <- RunUMAP(Macrophages, reduction = "pca", dims = 1:50)
ElbowPlot(Macrophages, ndims=50)

Macrophages <- FindNeighbors(Macrophages, dims = 1:8)
Macrophages <- FindClusters(Macrophages, resolution = 0.5)
Macrophages <- RunUMAP(Macrophages, dims = 1:8)
DimPlot(Macrophages)