##Extract SQ Immune cells
SQ_Myeloid <- subset(x = SQ_merge, idents = c('Myeloid cells'))

SQ_Myeloid <- ScaleData(SQ_Myeloid, features = VariableFeatures(SQ_Myeloid))
SQ_Myeloid <- RunPCA(SQ_Myeloid, features = VariableFeatures(SQ_Myeloid), npcs = 100)
SQ_Myeloid <- RunUMAP(SQ_Myeloid, reduction = "pca", dims = 1:50)
ElbowPlot(SQ_Myeloid, ndims=50)

SQ_Myeloid <- FindNeighbors(SQ_Myeloid, dims = 1:8)
SQ_Myeloid <- FindClusters(SQ_Myeloid, resolution = 0.3)
SQ_Myeloid <- RunUMAP(SQ_Myeloid, dims = 1:8)
DimPlot(SQ_Myeloid)

SQ_Myeloid_markers <- FindAllMarkers(SQ_Myeloid, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SQ_Myeloid_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

SQ_Myeloid_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> SQ_Myeloid_top10
DoHeatmap(SQ_Myeloid, features = SQ_Myeloid_top10$gene) + NoLegend()


SQ_Myeloid_IDs <- c("LAMs", "cDC2Bs", "Monocytes", "PVMs", "cDC1s", "LAMs", "cDC2As")
names(SQ_Myeloid_IDs) <- levels(SQ_Myeloid)
SQ_Myeloid <- RenameIdents(SQ_Myeloid, SQ_Myeloid_IDs)
DimPlot(SQ_Myeloid, reduction = "umap", pt.size = 0.5)
DimPlot(SQ_Myeloid, reduction = "umap", pt.size = 0.5) + NoLegend()

SQ_Myeloid_levels <- c('Monocytes', 'LAMs', 'PVMs', 'cDC1s', 'cDC2As', 'cDC2Bs')

Idents(SQ_Myeloid) <- factor(Idents(SQ_Myeloid), levels= SQ_Myeloid_levels)

SQ_Myeloid_markers <- FindAllMarkers(SQ_Myeloid, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SQ_Myeloid_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

SQ_Myeloid_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> SQ_Myeloid_top10
DoHeatmap(SQ_Myeloid, features = SQ_Myeloid_top10$gene) + NoLegend()

write.table(SQ_Myeloid_markers, file="SQ_Myeloid_markers.csv", sep=",")

##Split Data sets
SQ_Myeloid_Copy <- SQ_Myeloid
Idents(SQ_Myeloid_Copy) <- "Sample"
DimPlot(SQ_Myeloid_Copy)

##Isolate Healthy Obese Group
SQ_Myeloid_MHO_Copy <- subset(SQ_Myeloid_Copy, idents = c("SQ_MHO"))
Idents(SQ_Myeloid_MHO_Copy) <- "seurat_clusters"
DimPlot(SQ_Myeloid_MHO_Copy)

##Isolate unhealthy obese group
SQ_Myeloid_MUO_Copy <- subset(SQ_Myeloid_Copy, idents = c("SQ_MUO"))
Idents(SQ_Myeloid_MUO_Copy) <- "seurat_clusters"
DimPlot(SQ_Myeloid_MUO_Copy, label = T)

##Relabel Clusters MHO then MUO (Same code)
SQ_MHOMUO_Myeloid_IDs <- c("LAMs", "cDC2Bs", "Monocytes", "PVMs", "cDC1s", "LAMs", "cDC2As")

names(SQ_MHOMUO_Myeloid_IDs) <- levels(SQ_Myeloid_MHO_Copy)
names(SQ_MHOMUO_Myeloid_IDs) <- levels(SQ_Myeloid_MUO_Copy)

SQ_Myeloid_MHO_Copy <- RenameIdents(SQ_Myeloid_MHO_Copy, SQ_MHOMUO_Myeloid_IDs)
DimPlot(SQ_Myeloid_MHO_Copy, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

SQ_Myeloid_MUO_Copy <- RenameIdents(SQ_Myeloid_MUO_Copy, SQ_MHOMUO_Myeloid_IDs)
DimPlot(SQ_Myeloid_MUO_Copy, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


SQ_Myeloid_colors <- c('LAMs' = 'deeppink', 'cDC2As' = 'peru', 'cDC2Bs' = 'chartreuse3', 'PVMs' = 'skyblue2', 'Monocytes' = 'purple', 'cDC1s' = 'dodgerblue')

SCpubr::do_DimPlot(sample = SQ_Myeloid,
                   colors.use = SQ_Myeloid_colors, plot_cell_borders = TRUE, border.size = 2, border.color = 'black', pt.size = 0.8, label = TRUE, repel = TRUE)

SCpubr::do_DimPlot(sample = SQ_Myeloid_MHO_Copy,
                   colors.use = SQ_Myeloid_colors, plot_cell_borders = TRUE, border.size = 2, border.color = 'black', pt.size = 0.8, label = TRUE, repel = TRUE)

SCpubr::do_DimPlot(sample = SQ_Myeloid_MUO_Copy,
                   colors.use = SQ_Myeloid_colors, plot_cell_borders = TRUE, border.size = 2, border.color = 'black', pt.size = 0.8, label = TRUE, repel = TRUE)

vln_plot <- VlnPlot(SQ_Myeloid, 
                    features = c("FCER1A", "LAMP3", "TACSTD2", "LYVE1", "TREM2", "FCGR3A", "PTPRC"),
                    stack = TRUE, 
                    flip = TRUE, 
                    cols = c("TREM2" = "deeppink", "LAMP3" = "peru", "FCER1A" = "chartreuse3", "LYVE1" = "skyblue2", "FCGR3A" = "purple", "TACSTD2" = "dodgerblue", "PTPRC" = "gainsboro"))

# Remove the legend
vln_plot <- vln_plot + theme(legend.position = "none")
print(vln_plot)

##Obtain cell counts
SQ_Myeloid_MHO_Cells <- SQ_Myeloid_MHO_Copy@active.ident %>% as.data.table
View(SQ_Myeloid_MHO_Cells)
write.xlsx2(SQ_Myeloid_MHO_Cells, 'SQ_Myeloid_MHO_Counts.xlsx', sheetName = "Sheet1",
            col.names = TRUE, row.names = TRUE, append = FALSE)

SQ_Myeloid_MUO_Cells <- SQ_Myeloid_MUO_Copy@active.ident %>% as.data.table
View(SQ_Myeloid_MUO_Cells)
write.xlsx2(SQ_Myeloid_MUO_Cells, 'SQ_Myeloid_MUO_Counts.xlsx', sheetName = "Sheet1",
            col.names = TRUE, row.names = TRUE, append = FALSE)

## DEG analysis
SQ_Mono_Neu <- subset(SQ_Myeloid, idents = c("Mono_Neu"))
DimPlot(SQ_Mono_Neu, label = T)

SQ_Mono_Neu <- FindNeighbors(SQ_Mono_Neu, dims = 1:6)
SQ_Mono_Neu <- FindClusters(SQ_Mono_Neu, resolution = 0.2)
SQ_Mono_Neu <- RunUMAP(SQ_Mono_Neu, dims = 1:6)
DimPlot(SQ_Mono_Neu)
DimPlot(SQ_Mono_Neu, split.by = "Sample")

SQ_Mono_Neu_IDs <- c("Neu", "Mono")
names(SQ_Mono_Neu_IDs ) <- levels(SQ_Mono_Neu)
SQ_Mono_Neu <- RenameIdents(SQ_Mono_Neu, SQ_Mono_Neu_IDs )
DimPlot(SQ_Mono_Neu, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(SQ_Mono_Neu, reduction = "umap", pt.size = 0.5) + NoLegend()

SQ_Mono <- subset(SQ_Mono_Neu, idents = c("Mono"))
DimPlot(SQ_Mono, label = T)

SQ_Neu <- subset(SQ_Mono_Neu, idents = c("Neu"))
DimPlot(SQ_Neu, label = T)


SQ_Mono_DEGs <- FindMarkers(SQ_Mono, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 =  "SQ_MUO")
SQ_Mono_DEGs <- subset(SQ_Mono_DEGs, p_val_adj < 0.05)
write.table(SQ_Mono_DEGs, file="SQ_Mono_DEGs", sep=",")

SQ_Neu_DEGs <- FindMarkers(SQ_Neu, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 =  "SQ_MUO")
SQ_Neu_DEGs <- subset(SQ_Neu_DEGs, p_val_adj < 0.05)
write.table(SQ_Neu_DEGs, file="SQ_Neu_DEGs", sep=",")



## DEG analysis
SQ_LAM <- subset(SQ_Myeloid, idents = c("LAM"))
DimPlot(SQ_LAM, label = T)
SQ_LAM_DEGs <- FindMarkers(SQ_LAM, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_LAM_DEGs <- subset(SQ_LAM_DEGs, p_val_adj < 0.05)
write.table(SQ_LAM_DEGs, file="SQ_LAM_DEGs", sep=",")



## DEG analysis
SQ_PVM <- subset(SQ_Myeloid, idents = c("PVM"))
SQ_PVM_DEGs <- FindMarkers(SQ_PVM, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_PVM_DEGs <- subset(SQ_PVM_DEGs, p_val_adj < 0.05)
write.table(SQ_PVM_DEGs, file="SQ_PVM_DEGs", sep=",")



## DEG analysis
SQ_DCs <- subset(SQ_Myeloid, idents = c("cDC2A", "cDC2B", "cDC1"))
SQ_DCs_DEGs <- FindMarkers(SQ_DCs, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_DCs_DEGs <- subset(SQ_DCs_DEGs, p_val_adj < 0.05)
write.table(SQ_DCs_DEGs, file="SQ_DCs_DEGs", sep=",")


## DEG analysis
SQ_cDC1 <- subset(SQ_Myeloid, idents = c("cDC1"))
SQ_cDC1_DEGs <- FindMarkers(SQ_cDC1, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_cDC1_DEGs <- subset(SQ_cDC1_DEGs, p_val_adj < 0.05)
write.table(SQ_cDC1_DEGs, file="SQ_cDC1_DEGs", sep=",")


## DEG analysis
SQ_cDC2A <- subset(SQ_Myeloid, idents = c("cDC2A"))
SQ_cDC2A_DEGs <- FindMarkers(SQ_cDC2A, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_cDC2A_DEGs <- subset(SQ_cDC2A_DEGs, p_val_adj < 0.05)
write.table(SQ_cDC2A_DEGs, file="SQ_cDC2A_DEGs", sep=",")

## DEG analysis
SQ_cDC2B <- subset(SQ_Myeloid, idents = c("cDC2B"))
SQ_cDC2B_DEGs <- FindMarkers(SQ_cDC2B, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_cDC2B_DEGs <- subset(SQ_cDC2B_DEGs, p_val_adj < 0.05)
write.table(SQ_cDC2B_DEGs, file="SQ_cDC2B_DEGs", sep=",")

