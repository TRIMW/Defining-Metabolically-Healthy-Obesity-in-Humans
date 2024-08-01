##SUbQ PDGFRA and PI16 population clustering
SQ_PI16 <- ScaleData(SQ_PI16, features = VariableFeatures(SQ_PI16))
SQ_PI16 <- RunPCA(SQ_PI16, features = VariableFeatures(SQ_PI16), npcs = 100)
SQ_PI16 <- RunUMAP(SQ_PI16, reduction = "pca", dims = 1:50)
ElbowPlot(SQ_PI16, ndims=50)

SQ_PI16 <- FindNeighbors(SQ_PI16, dims = 1:10)
SQ_PI16 <- FindClusters(SQ_PI16, resolution = 0.2)
SQ_PI16 <- RunUMAP(SQ_PI16, dims = 1:10)
DimPlot(SQ_PI16)

SQ_PI16_markers <- FindAllMarkers(SQ_PI16, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SQ_PI16_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

SQ_PI16_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> SQ_PI16_top10
DoHeatmap(SQ_PI16, features = SQ_PI16_top10$gene) + NoLegend()

OM_PI16 <- ScaleData(OM_PI16, features = VariableFeatures(OM_PI16))
OM_PI16 <- RunPCA(OM_PI16, features = VariableFeatures(OM_PI16), npcs = 100)
OM_PI16 <- RunUMAP(OM_PI16, reduction = "pca", dims = 1:50)
ElbowPlot(OM_PI16, ndims=50)

OM_PI16 <- FindNeighbors(OM_PI16, dims = 1:10)
OM_PI16 <- FindClusters(OM_PI16, resolution = 0.2)
OM_PI16 <- RunUMAP(OM_PI16, dims = 1:10)
DimPlot(OM_PI16)

OM_PI16_markers <- FindAllMarkers(OM_PI16, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
OM_PI16_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

OM_PI16_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> OM_PI16_top5
DoHeatmap(OM_PI16, features = OM_PI16_top5$gene) + NoLegend()

##Subset PI16 from PDGFRA
OM_PI16_pos <- subset(x = OM_PI16, idents = c('2'))
SQ_PI16_pos <- subset(x = SQ_PI16, idents = c('3'))

OM_PDGFRA <- subset(x = OM_PI16, idents = c('0', '1', '3', '4'))
SQ_PDGFRA <- subset(x = SQ_PI16, idents = c('0', '1', '2', '4'))

##Re-cluster each population
OM_PI16_pos <- ScaleData(OM_PI16_pos, features = VariableFeatures(OM_PI16_pos))
OM_PI16_pos <- RunPCA(OM_PI16_pos, features = VariableFeatures(OM_PI16_pos), npcs = 100)
OM_PI16_pos <- RunUMAP(OM_PI16_pos, reduction = "pca", dims = 1:50)
ElbowPlot(OM_PI16_pos, ndims=50)

OM_PI16_pos <- FindNeighbors(OM_PI16_pos, dims = 1:4)
OM_PI16_pos <- FindClusters(OM_PI16_pos, resolution = 0.2)
OM_PI16_pos <- RunUMAP(OM_PI16_pos, dims = 1:4)
DimPlot(OM_PI16_pos)

OM_PI16_pos_markers <- FindAllMarkers(OM_PI16_pos, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
OM_PI16_pos_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

OM_PI16_pos_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> OM_PI16_pos_top10
DoHeatmap(OM_PI16_pos, features = OM_PI16_pos_top10$gene) + NoLegend()

write.table(OM_PI16_pos_markers, file="OM_PI16_pos_markers.csv", sep=",")

SQ_PI16_pos <- ScaleData(SQ_PI16_pos, features = VariableFeatures(SQ_PI16_pos))
SQ_PI16_pos <- RunPCA(SQ_PI16_pos, features = VariableFeatures(SQ_PI16_pos), npcs = 100)
SQ_PI16_pos <- RunUMAP(SQ_PI16_pos, reduction = "pca", dims = 1:50)
ElbowPlot(SQ_PI16_pos, ndims=50)

SQ_PI16_pos <- FindNeighbors(SQ_PI16_pos, dims = 1:4)
SQ_PI16_pos <- FindClusters(SQ_PI16_pos, resolution = 0.2)
SQ_PI16_pos <- RunUMAP(SQ_PI16_pos, dims = 1:4)
DimPlot(SQ_PI16_pos)

write.table(SQ_PI16_pos_markers, file="SQ_PI16_pos_markers.csv", sep=",")

OM_PDGFRA <- ScaleData(OM_PDGFRA, features = VariableFeatures(OM_PDGFRA))
OM_PDGFRA <- RunPCA(OM_PDGFRA, features = VariableFeatures(OM_PDGFRA), npcs = 100)
OM_PDGFRA <- RunUMAP(OM_PDGFRA, reduction = "pca", dims = 1:50)
ElbowPlot(OM_PDGFRA, ndims=50)

OM_PDGFRA <- FindNeighbors(OM_PDGFRA, dims = 1:5)
OM_PDGFRA <- FindClusters(OM_PDGFRA, resolution = 0.2)
OM_PDGFRA <- RunUMAP(OM_PDGFRA, dims = 1:5)
DimPlot(OM_PDGFRA)

OM_PDGFRA_markers <- FindAllMarkers(OM_PDGFRA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
OM_PDGFRA_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

OM_PDGFRA_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> OM_PDGFRA_top10
DoHeatmap(OM_PDGFRA, features = OM_PDGFRA_top10$gene) + NoLegend()

write.table(OM_PDGFRA_markers, file="OM_PDGFRA_markers.csv", sep=",")

SQ_PDGFRA <- ScaleData(SQ_PDGFRA, features = VariableFeatures(SQ_PDGFRA))
SQ_PDGFRA <- RunPCA(SQ_PDGFRA, features = VariableFeatures(SQ_PDGFRA), npcs = 100)
SQ_PDGFRA <- RunUMAP(SQ_PDGFRA, reduction = "pca", dims = 1:50)
ElbowPlot(SQ_PDGFRA, ndims=50)

SQ_PDGFRA <- FindNeighbors(SQ_PDGFRA, dims = 1:5)
SQ_PDGFRA <- FindClusters(SQ_PDGFRA, resolution = 0.2)
SQ_PDGFRA <- RunUMAP(SQ_PDGFRA, dims = 1:5)
DimPlot(SQ_PDGFRA)

SQ_PDGFRA_markers <- FindAllMarkers(SQ_PDGFRA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SQ_PDGFRA_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

SQ_PDGFRA_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> SQ_PDGFRA_top10
DoHeatmap(SQ_PDGFRA, features = SQ_PDGFRA_top10$gene) + NoLegend()

write.table(SQ_PDGFRA_markers, file="SQ_PDGFRA_markers.csv", sep=",")

##DEGs
SQ_PI16_pos_DEGs <- FindMarkers(SQ_PI16_pos, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_PI16_pos_DEGs <- subset(SQ_PI16_pos_DEGs, p_val_adj < 0.05)
write.table(SQ_PI16_pos_DEGs, file="SQ_PI16_pos_DEGs.csv", sep=",")

OM_PI16_pos_DEGs <- FindMarkers(OM_PI16_pos, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_PI16_pos_DEGs <- subset(OM_PI16_pos_DEGs, p_val_adj < 0.05)
write.table(OM_PI16_pos_DEGs, file="OM_PI16_pos_DEGs.csv", sep=",")

SQ_PDGFRA_DEGs <- FindMarkers(SQ_PDGFRA, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_PDGFRA_DEGs <- subset(SQ_PDGFRA_DEGs, p_val_adj < 0.05)
write.table(SQ_PDGFRA_DEGs, file="SQ_PDGFRA_DEGs.csv", sep=",")

OM_PDGFRA_DEGs <- FindMarkers(OM_PDGFRA, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_PDGFRA_DEGs <- subset(OM_PDGFRA_DEGs, p_val_adj < 0.05)
write.table(OM_PDGFRA_DEGs, file="OM_PDGFRA_DEGs.csv", sep=",")

##Metabolic module scoring from HALLMARK Pathways
OM_PI16 <- AddModuleScore(
  object = OM_PI16,
  features = HALLMARK_FA_METAB,
  ctrl = 5,
  name = 'HALLMARK_FA_METAB_features'
)

OM_PI16 <- AddModuleScore(
  object = OM_PI16,
  features = HALLMARK_GLYCOLYSIS,
  ctrl = 5,
  name = 'HALLMARK_GLYCOLYSIS_features'
)

OM_PI16 <- AddModuleScore(
  object = OM_PI16,
  features = HALLMARK_OXPHOS,
  ctrl = 5,
  name = 'HALLMARK_OXPHOS_features'
)

SQ_PI16 <- AddModuleScore(
  object = SQ_PI16,
  features = HALLMARK_FA_METAB,
  ctrl = 5,
  name = 'HALLMARK_FA_METAB_features'
)

SQ_PI16 <- AddModuleScore(
  object = SQ_PI16,
  features = HALLMARK_GLYCOLYSIS,
  ctrl = 5,
  name = 'HALLMARK_GLYCOLYSIS_features'
)

SQ_PI16 <- AddModuleScore(
  object = SQ_PI16,
  features = HALLMARK_OXPHOS,
  ctrl = 5,
  name = 'HALLMARK_OXPHOS_features'
)

OM_PDGFRA <- AddModuleScore(
  object = OM_PDGFRA,
  features = HALLMARK_FA_METAB,
  ctrl = 5,
  name = 'HALLMARK_FA_METAB_features'
)

OM_PDGFRA <- AddModuleScore(
  object = OM_PDGFRA,
  features = HALLMARK_GLYCOLYSIS,
  ctrl = 5,
  name = 'HALLMARK_GLYCOLYSIS_features'
)

OM_PDGFRA <- AddModuleScore(
  object = OM_PDGFRA,
  features = HALLMARK_OXPHOS,
  ctrl = 5,
  name = 'HALLMARK_OXPHOS_features'
)

SQ_PDGFRA <- AddModuleScore(
  object = SQ_PDGFRA,
  features = HALLMARK_FA_METAB,
  ctrl = 5,
  name = 'HALLMARK_FA_METAB_features'
)

SQ_PDGFRA <- AddModuleScore(
  object = SQ_PDGFRA,
  features = HALLMARK_GLYCOLYSIS,
  ctrl = 5,
  name = 'HALLMARK_GLYCOLYSIS_features'
)

SQ_PDGFRA <- AddModuleScore(
  object = SQ_PDGFRA,
  features = HALLMARK_OXPHOS,
  ctrl = 5,
  name = 'HALLMARK_OXPHOS_features'
)

OM_MCs <- AddModuleScore(
  object = OM_MCs,
  features = HALLMARK_FA_METAB,
  ctrl = 5,
  name = 'HALLMARK_FA_METAB_features'
)

OM_MCs <- AddModuleScore(
  object = OM_MCs,
  features = HALLMARK_GLYCOLYSIS,
  ctrl = 5,
  name = 'HALLMARK_GLYCOLYSIS_features'
)

OM_MCs <- AddModuleScore(
  object = OM_MCs,
  features = HALLMARK_OXPHOS,
  ctrl = 5,
  name = 'HALLMARK_OXPHOS_features'
)

OM_PI16 <- AddModuleScore(
  object = OM_PI16,
  features = HALLMARK_ADIPOGENESIS,
  ctrl = 5,
  name = 'HALLMARK_ADIPOGENESIS_features'
)

SQ_PI16 <- AddModuleScore(
  object = SQ_PI16,
  features = HALLMARK_ADIPOGENESIS,
  ctrl = 5,
  name = 'HALLMARK_ADIPOGENESIS_features'
)

OM_PDGFRA <- AddModuleScore(
  object = OM_PDGFRA,
  features = HALLMARK_ADIPOGENESIS,
  ctrl = 5,
  name = 'HALLMARK_ADIPOGENESIS_features'
)

SQ_PDGFRA <- AddModuleScore(
  object = SQ_PDGFRA,
  features = HALLMARK_ADIPOGENESIS,
  ctrl = 5,
  name = 'HALLMARK_ADIPOGENESIS_features'
)

VlnPlot(object = OM_PI16, features = 'HALLMARK_GLYCOLYSIS_features1', split.by = 'Sample', pt.size = 0)
ggsave("OM_PI16_GLYCOLYSIS.tiff", path = '/Users/willtrim/Library/CloudStorage/OneDrive-HarvardUniversity/Manuscripts & Conferences/Manuscripts/Metabolic Obesity/Data/scRNAseq/R Image Outputs/', dpi=500)
VlnPlot(object = OM_PI16, features = 'HALLMARK_OXPHOS_features1', split.by = 'Sample', pt.size = 0)
ggsave("OM_PI16_OXPHOS.tiff", path = '/Users/willtrim/Library/CloudStorage/OneDrive-HarvardUniversity/Manuscripts & Conferences/Manuscripts/Metabolic Obesity/Data/scRNAseq/R Image Outputs/', dpi=500)
VlnPlot(object = OM_PI16, features = 'HALLMARK_FA_METAB_features1', split.by = 'Sample', pt.size = 0)
ggsave("OM_PI16_FA_METAB.tiff", path = '/Users/willtrim/Library/CloudStorage/OneDrive-HarvardUniversity/Manuscripts & Conferences/Manuscripts/Metabolic Obesity/Data/scRNAseq/R Image Outputs/', dpi=500)

VlnPlot(object = OM_PDGFRA, features = 'HALLMARK_GLYCOLYSIS_features1', split.by = 'Sample', pt.size = 0)
ggsave("OM_PDGFRA_GLYCOLYSIS.tiff", path = '/Users/willtrim/Library/CloudStorage/OneDrive-HarvardUniversity/Manuscripts & Conferences/Manuscripts/Metabolic Obesity/Data/scRNAseq/R Image Outputs/', dpi=500)
VlnPlot(object = OM_PDGFRA, features = 'HALLMARK_OXPHOS_features1', split.by = 'Sample', pt.size = 0)
ggsave("OM_PDGFRA_OXPHOS.tiff", path = '/Users/willtrim/Library/CloudStorage/OneDrive-HarvardUniversity/Manuscripts & Conferences/Manuscripts/Metabolic Obesity/Data/scRNAseq/R Image Outputs/', dpi=500)
VlnPlot(object = OM_PDGFRA, features = 'HALLMARK_FA_METAB_features1', split.by = 'Sample', pt.size = 0)
ggsave("OM_PDGFRA_FA_METAB.tiff", path = '/Users/willtrim/Library/CloudStorage/OneDrive-HarvardUniversity/Manuscripts & Conferences/Manuscripts/Metabolic Obesity/Data/scRNAseq/R Image Outputs/', dpi=500)

VlnPlot(object = SQ_PI16, features = 'HALLMARK_GLYCOLYSIS_features1', split.by = 'Sample', pt.size = 0)
ggsave("SQ_PI16_GLYCOLYSIS.tiff", path = '/Users/willtrim/Library/CloudStorage/OneDrive-HarvardUniversity/Manuscripts & Conferences/Manuscripts/Metabolic Obesity/Data/scRNAseq/R Image Outputs/', dpi=500)
VlnPlot(object = SQ_PI16, features = 'HALLMARK_OXPHOS_features1', split.by = 'Sample', pt.size = 0)
ggsave("SQ_PI16_OXPHOS.tiff", path = '/Users/willtrim/Library/CloudStorage/OneDrive-HarvardUniversity/Manuscripts & Conferences/Manuscripts/Metabolic Obesity/Data/scRNAseq/R Image Outputs/', dpi=500)
VlnPlot(object = SQ_PI16, features = 'HALLMARK_FA_METAB_features1', split.by = 'Sample', pt.size = 0)
ggsave("SQ_PI16_FA_METAB.tiff", path = '/Users/willtrim/Library/CloudStorage/OneDrive-HarvardUniversity/Manuscripts & Conferences/Manuscripts/Metabolic Obesity/Data/scRNAseq/R Image Outputs/', dpi=500)

VlnPlot(object = SQ_PDGFRA, features = 'HALLMARK_GLYCOLYSIS_features1', split.by = 'Sample', pt.size = 0)
ggsave("SQ_PDGFRA_GLYCOLYSIS.tiff", path = '/Users/willtrim/Library/CloudStorage/OneDrive-HarvardUniversity/Manuscripts & Conferences/Manuscripts/Metabolic Obesity/Data/scRNAseq/R Image Outputs/', dpi=500)
VlnPlot(object = SQ_PDGFRA, features = 'HALLMARK_OXPHOS_features1', split.by = 'Sample', pt.size = 0)
ggsave("SQ_PDGFRA_OXPHOS.tiff", path = '/Users/willtrim/Library/CloudStorage/OneDrive-HarvardUniversity/Manuscripts & Conferences/Manuscripts/Metabolic Obesity/Data/scRNAseq/R Image Outputs/', dpi=500)
VlnPlot(object = SQ_PDGFRA, features = 'HALLMARK_FA_METAB_features1', split.by = 'Sample', pt.size = 0)
ggsave("SQ_PDGFRA_FA_METAB.tiff", path = '/Users/willtrim/Library/CloudStorage/OneDrive-HarvardUniversity/Manuscripts & Conferences/Manuscripts/Metabolic Obesity/Data/scRNAseq/R Image Outputs/', dpi=500)

VlnPlot(object = OM_PI16, features = 'HALLMARK_ADIPOGENESIS_features1', split.by = 'Sample', pt.size = 0)
ggsave("OM_PI16_ADIPOGENESIS.tiff", path = '/Users/willtrim/Library/CloudStorage/OneDrive-HarvardUniversity/Manuscripts & Conferences/Manuscripts/Metabolic Obesity/Data/scRNAseq/R Image Outputs/', dpi=500)
VlnPlot(object = SQ_PI16, features = 'HALLMARK_ADIPOGENESIS_features1', split.by = 'Sample', pt.size = 0)
ggsave("SQ_PI16_ADIPOGENESIS.tiff", path = '/Users/willtrim/Library/CloudStorage/OneDrive-HarvardUniversity/Manuscripts & Conferences/Manuscripts/Metabolic Obesity/Data/scRNAseq/R Image Outputs/', dpi=500)

VlnPlot(object = OM_PDGFRA, features = 'HALLMARK_ADIPOGENESIS_features1', split.by = 'Sample', pt.size = 0)
ggsave("OM_PDGFRA_ADIPOGENESIS.tiff", path = '/Users/willtrim/Library/CloudStorage/OneDrive-HarvardUniversity/Manuscripts & Conferences/Manuscripts/Metabolic Obesity/Data/scRNAseq/R Image Outputs/', dpi=500)
VlnPlot(object = SQ_PDGFRA, features = 'HALLMARK_ADIPOGENESIS_features1', split.by = 'Sample', pt.size = 0)
ggsave("SQ_PDGFRA_ADIPOGENESIS.tiff", path = '/Users/willtrim/Library/CloudStorage/OneDrive-HarvardUniversity/Manuscripts & Conferences/Manuscripts/Metabolic Obesity/Data/scRNAseq/R Image Outputs/', dpi=500)

VlnPlot(object = OM_MCs, features = 'HALLMARK_GLYCOLYSIS_features1', split.by = 'Sample', pt.size = 0)
ggsave("OM_MC_GLYCOLYSIS.tiff", path = '/Users/willtrim/Library/CloudStorage/OneDrive-HarvardUniversity/Manuscripts & Conferences/Manuscripts/Metabolic Obesity/Data/scRNAseq/R Image Outputs/', dpi=500)
VlnPlot(object = OM_MCs, features = 'HALLMARK_OXPHOS_features1', split.by = 'Sample', pt.size = 0)
ggsave("OM_MC_OXPHOS.tiff", path = '/Users/willtrim/Library/CloudStorage/OneDrive-HarvardUniversity/Manuscripts & Conferences/Manuscripts/Metabolic Obesity/Data/scRNAseq/R Image Outputs/', dpi=500)
VlnPlot(object = OM_MCs, features = 'HALLMARK_FA_METAB_features1', split.by = 'Sample', pt.size = 0)
ggsave("OM_MC_FA_METAB.tiff", path = '/Users/willtrim/Library/CloudStorage/OneDrive-HarvardUniversity/Manuscripts & Conferences/Manuscripts/Metabolic Obesity/Data/scRNAseq/R Image Outputs/', dpi=500)