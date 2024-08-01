##SQ Smooth muscle cells
SQ_SMC <- ScaleData(SQ_SMC, features = VariableFeatures(SQ_SMC))
SQ_SMC <- RunPCA(SQ_SMC, features = VariableFeatures(SQ_SMC), npcs = 100)
SQ_SMC <- RunUMAP(SQ_SMC, reduction = "pca", dims = 1:50)
ElbowPlot(SQ_SMC, ndims=50)

SQ_SMC <- FindNeighbors(SQ_SMC, dims = 1:10)
SQ_SMC <- FindClusters(SQ_SMC, resolution = 0.2)
SQ_SMC <- RunUMAP(SQ_SMC, dims = 1:10)
DimPlot(SQ_SMC)

SQ_SMC_markers <- FindAllMarkers(SQ_SMC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SQ_SMC_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

SQ_SMC_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> SQ_SMC_top10
DoHeatmap(SQ_SMC, features = SQ_SMC_top10$gene) + NoLegend()

SQ_SMC_DEGs <- FindMarkers(SQ_SMC, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_SMC_DEGs <- subset(SQ_SMC_DEGs, p_val_adj < 0.05)
write.table(SQ_SMC_DEGs, file="SQ_SMC_DEGs.csv", sep=",")

vln_plot <- VlnPlot(SQ_SMC, features = c("MTRNR2L1", "C1R", "ACTG2", "TXNIP", "TPM2", "RHOB", "MT1M", "C1S", "SOCS3", "NMB"), stack = TRUE, flip = TRUE)
vln_plot <- vln_plot + theme(legend.position = "none")
print(vln_plot)

##Om Smooth muscle cells
OM_SMC <- ScaleData(OM_SMC, features = VariableFeatures(OM_SMC))
OM_SMC <- RunPCA(OM_SMC, features = VariableFeatures(OM_SMC), npcs = 100)
OM_SMC <- RunUMAP(OM_SMC, reduction = "pca", dims = 1:50)
ElbowPlot(OM_SMC, ndims=50)

OM_SMC <- FindNeighbors(OM_SMC, dims = 1:8)
OM_SMC <- FindClusters(OM_SMC, resolution = 0.2)
OM_SMC <- RunUMAP(OM_SMC, dims = 1:8)
DimPlot(OM_SMC)

OM_SMC_markers <- FindAllMarkers(OM_SMC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
OM_SMC_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

OM_SMC_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> OM_SMC_top10
DoHeatmap(OM_SMC, features = OM_SMC_top10$gene) + NoLegend()

OM_SMC_DEGs <- FindMarkers(OM_SMC, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_SMC_DEGs <- subset(OM_SMC_DEGs, p_val_adj < 0.05)
write.table(OM_SMC_DEGs, file="OM_SMC_DEGs.csv", sep=",")

vln_plot <- VlnPlot(OM_SMC, features = c("JUN", "JUND", "ITLN", "S100A4", "SOD3", "IRI27", "RPS4Y1", "ECRG4", "MT1E"), stack = TRUE, flip = TRUE)
vln_plot <- vln_plot + theme(legend.position = "none")
print(vln_plot)


##SQ Endothelial cells
SQ_ECs <- ScaleData(SQ_ECs, features = VariableFeatures(SQ_ECs))
SQ_ECs <- RunPCA(SQ_ECs, features = VariableFeatures(SQ_ECs), npcs = 100)
SQ_ECs <- RunUMAP(SQ_ECs, reduction = "pca", dims = 1:50)
ElbowPlot(SQ_ECs, ndims=50)

SQ_ECs <- FindNeighbors(SQ_ECs, dims = 1:5)
SQ_ECs <- FindClusters(SQ_ECs, resolution = 0.2)
SQ_ECs <- RunUMAP(SQ_ECs, dims = 1:5)
DimPlot(SQ_ECs)

SQ_ECs_markers <- FindAllMarkers(SQ_ECs, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SQ_ECs_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

SQ_ECs_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> SQ_ECs_top10
DoHeatmap(SQ_ECs, features = SQ_ECs_top10$gene) + NoLegend()

SQ_ECs_DEGs <- FindMarkers(SQ_ECs, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_ECs_DEGs <- subset(SQ_ECs_DEGs, p_val_adj < 0.05)
write.table(SQ_ECs_DEGs, file="SQ_ECs_DEGs.csv", sep=",")

vln_plot <- VlnPlot(SQ_ECs, features = c("HLA-DRB5", "CYGB", "IGFBP4", "MYC", "FABP5", "MTRNR2L8", "TXNIP", "VWA1", "CD34"), stack = TRUE, flip = TRUE, split.by = 'Sample')
vln_plot <- vln_plot + theme(legend.position = "none")
print(vln_plot)

##OM Endothelial cells
OM_EC <- ScaleData(OM_EC, features = VariableFeatures(OM_EC))
OM_EC <- RunPCA(OM_EC, features = VariableFeatures(OM_EC), npcs = 100)
OM_EC <- RunUMAP(OM_EC, reduction = "pca", dims = 1:50)
ElbowPlot(OM_EC, ndims=50)

OM_EC <- FindNeighbors(OM_EC, dims = 1:6)
OM_EC <- FindClusters(OM_EC, resolution = 0.2)
OM_EC <- RunUMAP(OM_EC, dims = 1:6)
DimPlot(OM_EC)

OM_EC_markers <- FindAllMarkers(OM_EC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
OM_EC_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

OM_EC_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> OM_EC_top5
DoHeatmap(OM_EC, features = OM_EC_top5$gene) + NoLegend()

OM_EC_DEGs <- FindMarkers(OM_EC, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_EC_DEGs <- subset(OM_EC_DEGs, p_val_adj < 0.05)
write.table(OM_EC_DEGs, file="OM_EC_DEGs.csv", sep=",")

vln_plot <- VlnPlot(OM_EC, features = c("MT1X", "MT1E", "MT1M", "MT2A", "RPL41", "RPS8", "RPL14", "RPLP0", "RPL13"), stack = TRUE, flip = TRUE, split.by = 'Sample')
vln_plot <- vln_plot + theme(legend.position = "none")
print(vln_plot)
