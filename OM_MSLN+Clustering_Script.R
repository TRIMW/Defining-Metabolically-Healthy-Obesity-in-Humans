OM_MSCs <- subset(x = OM_merge, idents = c("MSC1s", "MSC2s", "MSC3s"))
OM_MSCs <- ScaleData(OM_MSCs, features = VariableFeatures(OM_MSCs))
OM_MSCs <- RunPCA(OM_MSCs, features = VariableFeatures(OM_MSCs), npcs = 100)
OM_MSCs <- RunUMAP(OM_MSCs, reduction = "pca", dims = 1:50)
OM_MSCs <- FindNeighbors(OM_MSCs, dims = 1:10)
OM_MSCs <- FindClusters(OM_MSCs, resolution = 0.5)
OM_MSCs <- RunUMAP(OM_MSCs, dims = 1:10)
DimPlot(OM_MSCs)

DimPlot(OM_MSCs)
OM_MSCs <- subset(x = OM_MSCs, idents = c("0", "1", "2", "3", "4", "5"))

DimPlot(OM_MSCs, split.by = "Sample")

OM_MSCs_markers <- FindAllMarkers(OM_MSCs, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
OM_MSCs_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

OM_MSCs_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> OM_MSCs_top10
DoHeatmap(OM_MSCs, features = OM_MSCs_top10$gene) + NoLegend()


VlnPlot(object = OM_MCs, features = 'HALLMARK_GLYCOLYSIS_features1', split.by = 'Sample', pt.size = 0)
ggsave("OM_MC_GLYCOLYSIS.tiff", path = '/Users/willtrim/Library/CloudStorage/OneDrive-HarvardUniversity/Manuscripts & Conferences/Manuscripts/Metabolic Obesity/Data/scRNAseq/R Image Outputs/', dpi=500)
VlnPlot(object = OM_MCs, features = 'HALLMARK_OXPHOS_features1', split.by = 'Sample', pt.size = 0)
ggsave("OM_MC_OXPHOS.tiff", path = '/Users/willtrim/Library/CloudStorage/OneDrive-HarvardUniversity/Manuscripts & Conferences/Manuscripts/Metabolic Obesity/Data/scRNAseq/R Image Outputs/', dpi=500)
VlnPlot(object = OM_MCs, features = 'HALLMARK_FA_METAB_features1', split.by = 'Sample', pt.size = 0)
ggsave("OM_MC_FA_METAB.tiff", path = '/Users/willtrim/Library/CloudStorage/OneDrive-HarvardUniversity/Manuscripts & Conferences/Manuscripts/Metabolic Obesity/Data/scRNAseq/R Image Outputs/', dpi=500)