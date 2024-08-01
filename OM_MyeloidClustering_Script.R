OM_Myeloid <- subset(x = OM_merge, idents = c('Myeloid cells'))

OM_Myeloid <- ScaleData(OM_Myeloid, features = VariableFeatures(OM_Myeloid))
OM_Myeloid <- RunPCA(OM_Myeloid, features = VariableFeatures(OM_Myeloid), npcs = 100)
OM_Myeloid <- RunUMAP(OM_Myeloid, reduction = "pca", dims = 1:50)
OM_Myeloid <- FindNeighbors(OM_Myeloid, dims = 1:10)
OM_Myeloid <- FindClusters(OM_Myeloid, resolution = 0.5)
OM_Myeloid <- RunUMAP(OM_Myeloid, dims = 1:10)
DimPlot(OM_Myeloid)
DimPlot(OM_Myeloid, split.by = "Sample") + theme(aspect.ratio = 1)


OM_Myeloid_markers <- FindAllMarkers(OM_Myeloid, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
OM_Myeloid_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

OM_Myeloid_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> OM_Myeloid_top10
DoHeatmap(OM_Myeloid, features = OM_Myeloid_top10$gene) + NoLegend()

OM_Myeloid <- subset(x = OM_Myeloid, idents = c("0", "1", "2", "3", "4", "5", "6"))

OM_Myeloid <- subset(x = OM_Myeloid, idents = c("0", "2", "3", "4", "5", "7", "8"))
OM_Myeloid <- ScaleData(OM_Myeloid, features = VariableFeatures(OM_Myeloid))
OM_Myeloid <- ScaleData(OM_Myeloid, features = VariableFeatures(OM_Myeloid))
OM_Myeloid <- RunPCA(OM_Myeloid, features = VariableFeatures(OM_Myeloid), npcs = 100)
OM_Myeloid <- RunUMAP(OM_Myeloid, reduction = "pca", dims = 1:50)
OM_Myeloid <- FindNeighbors(OM_Myeloid, dims = 1:10)
OM_Myeloid <- FindClusters(OM_Myeloid, resolution = 0.3)
OM_Myeloid <- RunUMAP(OM_Myeloid, dims = 1:10)
DimPlot(OM_Myeloid)
DimPlot(OM_Myeloid, split.by = "Sample") + theme(aspect.ratio = 1)


OM_Myeloid_markers <- FindAllMarkers(OM_Myeloid, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
OM_Myeloid_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

OM_Myeloid_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> OM_Myeloid_top10
DoHeatmap(OM_Myeloid, features = OM_Myeloid_top10$gene) + NoLegend()

#Naming clusters
OM_Myeloid_IDs <- c("THBS1+", "LAMs", "cDC2Bs", "Monocytes", "PVMs", "FBLN1+", "cDC2As")
names(OM_Myeloid_IDs) <- levels(OM_Myeloid)
OM_Myeloid <- RenameIdents(OM_Myeloid, OM_Myeloid_IDs)
DimPlot(OM_Myeloid, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

OM_Myeloid_levels <- c('THBS1+', 'FBLN1+', 'Monocytes', 'LAMs', 'PVMs', 'cDC2As', 'cDC2Bs')

Idents(OM_Myeloid) <- factor(Idents(OM_Myeloid), levels= OM_Myeloid_levels)



OM_Myeloid_colors <- c('FBLN1+' = '#005f73', 'LAMs' = 'deeppink', 'cDC2As' = 'peru', 'cDC2Bs' = 'chartreuse3', 'PVMs' = 'skyblue2', 'THBS1+' = 'gold', 'Monocytes' = 'purple')

SCpubr::do_DimPlot(sample = OM_Myeloid,
                   colors.use = OM_Myeloid_colors, plot_cell_borders = TRUE, border.size = 2, border.color = 'black', pt.size = 0.6, label = TRUE, repel = TRUE)

vln_plot <- VlnPlot(OM_Myeloid, 
                    features = c("CLEC10A", "CLEC9A", "LYVE1", "TREM2", "FCGR3A", "FBLN1", "THBS1", "PTPRC"),
                    stack = TRUE, 
                    flip = TRUE, 
                    cols = c("CLEC9A" = "peru", "FBLN1" = "#005f73", "LYVE1" = "skyblue2", "FCGR3A" = "purple", "CLEC10A" = "chartreuse3", "TREM2" = "deeppink", "THBS1" = "gold", "PTPRC" = "gainsboro"))

vln_plot <- vln_plot + theme(legend.position = "none")
print(vln_plot)

##Split Data sets
OM_Myeloid_Copy <- OM_Myeloid
Idents(OM_Myeloid_Copy) <- "Sample"
DimPlot(OM_Myeloid_Copy)

##Isolate Healthy Obese Group
OM_Myeloid_MHO_Copy <- subset(OM_Myeloid_Copy, idents = c("OM_MHO"))
Idents(OM_Myeloid_MHO_Copy) <- "seurat_clusters"
DimPlot(OM_Myeloid_MHO_Copy)

##Isolate unhealthy obese group
OM_Myeloid_MUO_Copy <- subset(OM_Myeloid_Copy, idents = c("OM_MUO"))
Idents(OM_Myeloid_MUO_Copy) <- "seurat_clusters"
DimPlot(OM_Myeloid_MUO_Copy, label = T)

##Relabel Clusters MHO then MUO (Same code)
OM_MHOMUO_Myeloid_IDs <- c("THBS1+", "LAMs", "THBS1+", "cDC2Bs",  "Monocytes", "PVMs", "FBLN1+", "cDC2As")

names(OM_MHOMUO_Myeloid_IDs) <- levels(OM_Myeloid_MHO_Copy)
names(OM_MHOMUO_Myeloid_IDs) <- levels(OM_Myeloid_MUO_Copy)

OM_Myeloid_MHO_Copy <- RenameIdents(OM_Myeloid_MHO_Copy, OM_MHOMUO_Myeloid_IDs)
DimPlot(OM_Myeloid_MHO_Copy, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

OM_Myeloid_MUO_Copy <- RenameIdents(OM_Myeloid_MUO_Copy, OM_MHOMUO_Myeloid_IDs)
DimPlot(OM_Myeloid_MUO_Copy, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


##Re-order clusters
OM_Myeloid_levels <- c('THBS1+', 'FBLN1+', 'Monocytes', 'LAMs', 'PVMs', 'cDC2As', 'cDC2Bs')

Idents(OM_Myeloid_MHO_Copy) <- factor(Idents(OM_Myeloid_MHO_Copy), levels= OM_Myeloid_levels)
Idents(OM_Myeloid_MUO_Copy) <- factor(Idents(OM_Myeloid_MUO_Copy), levels= OM_Myeloid_levels)

SCpubr::do_DimPlot(sample = OM_Myeloid_MHO_Copy,
                   colors.use = OM_Myeloid_colors, plot_cell_borders = TRUE, border.size = 2, border.color = 'black', pt.size = 0.8, label = TRUE, repel = TRUE)
SCpubr::do_DimPlot(sample = OM_Myeloid_MUO_Copy,
                   colors.use = OM_Myeloid_colors, plot_cell_borders = TRUE, border.size = 2, border.color = 'black', pt.size = 0.8, label = TRUE, repel = TRUE)

##Get cell counts
OM_Myeloid_MHO_Cells <- OM_Myeloid_MHO_Copy@active.ident %>% as.data.table
View(OM_Myeloid_MHO_Cells)
write.xlsx2(OM_Myeloid_MHO_Cells, 'OM_Myeloid_MHO_Counts.xlsx', sheetName = "Sheet1",
            col.names = TRUE, row.names = TRUE, append = FALSE)

OM_Myeloid_MUO_Cells <- OM_Myeloid_MUO_Copy@active.ident %>% as.data.table
View(OM_Myeloid_MUO_Cells)
write.xlsx2(OM_Myeloid_MUO_Cells, 'OM_Myeloid_MUO_Counts.xlsx', sheetName = "Sheet1",
            col.names = TRUE, row.names = TRUE, append = FALSE)

## DEG analysis
OM_Mono <- subset(OM_Myeloid, idents = c("Mono"))
DimPlot(OM_Mono, label = T)
OM_Mono_DEGs <- FindMarkers(OM_Mono, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_Mono_DEGs <- subset(OM_Mono_DEGs, p_val_adj < 0.05)
write.table(OM_Mono_DEGs, file="OM_Mono_DEGs", sep=",")


## DEG analysis
OM_THBS1 <- subset(OM_Myeloid, idents = c("THBS1+"))
DimPlot(OM_THBS1, label = T)
OM_THBS1_DEGs <- FindMarkers(OM_THBS1, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_THBS1_DEGs <- subset(OM_THBS1_DEGs, p_val_adj < 0.05)
write.table(OM_THBS1_DEGs, file="OM_THBS1_DEGs", sep=",")

OM_Neu <- ScaleData(OM_Neu, features = VariableFeatures(OM_Neu))
OM_Neu <- RunPCA(OM_Neu, features = VariableFeatures(OM_Neu), npcs = 100)
OM_Neu <- RunUMAP(OM_Neu, reduction = "pca", dims = 1:50)
ElbowPlot(OM_Neu, ndims=50)
OM_Neu <- FindNeighbors(OM_Neu, dims = 1:6)
OM_Neu <- FindClusters(OM_Neu, resolution = 0.2)
OM_Neu <- RunUMAP(OM_Neu, dims = 1:6)
DimPlot(OM_Neu)
DimPlot(OM_Neu, split.by = "Sample")

OM_Neu_markers <- FindAllMarkers(OM_Neu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
OM_Neu_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

OM_Neu_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> OM_Neu_top10
DoHeatmap(OM_Neu, features = OM_Neu_top10$gene) + NoLegend()


OM_PVM <- subset(OM_Myeloid, idents = c("PVM"))
DimPlot(OM_PVM, label = T)
OM_PVM_DEGs <- FindMarkers(OM_PVM, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_PVM_DEGs <- subset(OM_PVM_DEGs, p_val_adj < 0.05)
write.table(OM_PVM_DEGs, file="OM_PVM_DEGs", sep=",")

OM_LAM <- subset(OM_Myeloid, idents = c("LAM"))
DimPlot(OM_LAM, label = T)
OM_LAM_DEGs <- FindMarkers(OM_LAM, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_LAM_DEGs <- subset(OM_LAM_DEGs, p_val_adj < 0.05)
write.table(OM_LAM_DEGs, file="OM_LAM_DEGs", sep=",")

OM_DCs <- subset(OM_Myeloid, idents = c("cDC2B", "cDC2A"))
DimPlot(OM_DCs, label = T)
OM_DCs_DEGs <- FindMarkers(OM_DCs, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_DCs_DEGs <- subset(OM_DCs_DEGs, p_val_adj < 0.05)
write.table(OM_DCs_DEGs, file="OM_DCs_DEGs", sep=",")

OM_DCs <- ScaleData(OM_DCs, features = VariableFeatures(OM_DCs))
OM_DCs <- RunPCA(OM_DCs, features = VariableFeatures(OM_DCs), npcs = 100)
OM_DCs <- RunUMAP(OM_DCs, reduction = "pca", dims = 1:50)
ElbowPlot(OM_Neu, ndims=50)
OM_DCs <- FindNeighbors(OM_DCs, dims = 1:6)
OM_DCs <- FindClusters(OM_DCs, resolution = 0.2)
OM_DCs <- RunUMAP(OM_DCs, dims = 1:6)
DimPlot(OM_DCs)
DimPlot(OM_DCs, split.by = "Sample")

OM_DCs_markers <- FindAllMarkers(OM_DCs, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
OM_DCs_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

OM_DCs_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> OM_DCs_top10
DoHeatmap(OM_DCs, features = OM_DCs_top10$gene) + NoLegend()

OM_cDC2B <- subset(OM_Myeloid, idents = c("cDC2B"))
DimPlot(OM_cDC2B, label = T)
OM_cDC2B_DEGs <- FindMarkers(OM_cDC2B, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_cDC2B_DEGs <- subset(OM_cDC2B_DEGs, p_val_adj < 0.05)
write.table(OM_cDC2B_DEGs, file="OM_cDC2B_DEGs", sep=",")

OM_cDC2A <- subset(OM_Myeloid, idents = c("cDC2A"))
DimPlot(OM_cDC2A, label = T)
OM_cDC2A_DEGs <- FindMarkers(OM_cDC2A, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_cDC2A_DEGs <- subset(OM_cDC2A_DEGs, p_val_adj < 0.05)
write.table(OM_cDC2A_DEGs, file="OM_cDC2A_DEGs", sep=",")