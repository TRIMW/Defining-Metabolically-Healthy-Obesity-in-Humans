

BRI1456.data <- Read10X(data.dir = "/Users/willtrim/Library/CloudStorage/OneDrive-TCDUD.onmicrosoft.com/Post doc Trinity Harvard/Data/Two Depots Study/Data/scRNAseq - BTM/Data/GRCh38/BRI-1456/outs/raw_feature_bc_matrix")
BRI1456 <- CreateSeuratObject(counts = BRI1456.data, project = "SQfat10x", min.cells = 3, min.features = 200)
BRI1456$tissue <-"SQ"
BRI1456$disease <-"Healthy"
BRI1456$run <-"1456"
BRI1456[["percent.mt"]] <- PercentageFeatureSet(BRI1456, pattern = "^MT-")
VlnPlot(BRI1456, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(BRI1456, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(BRI1456, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
BRI1456 <- subset(BRI1456, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 20 & percent.mt > 1)
counts_BRI1456 <- GetAssayData(BRI1456, assay = "RNA", slot = "counts")

BRI1457.data <- Read10X(data.dir = "/Users/willtrim/Library/CloudStorage/OneDrive-TCDUD.onmicrosoft.com/Post doc Trinity Harvard/Data/Two Depots Study/Data/scRNAseq - BTM/Data/GRCh38/BRI-1457/outs/raw_feature_bc_matrix")
BRI1457 <- CreateSeuratObject(counts = BRI1457.data, project = "SQfat10x", min.cells = 3, min.features = 200)
BRI1457
BRI1457$tissue <-"SQ"
BRI1457$disease <-"Unhealthy"
BRI1457$run <-"1457"
BRI1457[["percent.mt"]] <- PercentageFeatureSet(BRI1457, pattern = "^MT-")
VlnPlot(BRI1457, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(BRI1457, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(BRI1457, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
BRI1457 <- subset(BRI1457, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 16.5 & percent.mt > 3)
counts_BRI1457 <- GetAssayData(BRI1457, assay = "RNA", slot = "counts")

BRI1458.data <- Read10X(data.dir = "/Users/willtrim/Library/CloudStorage/OneDrive-TCDUD.onmicrosoft.com/Post doc Trinity Harvard/Data/Two Depots Study/Data/scRNAseq - BTM/Data/GRCh38/BRI-1458/outs/raw_feature_bc_matrix")
BRI1458 <- CreateSeuratObject(counts = BRI1458.data, project = "OMfat10x", min.cells = 3, min.features = 200)
BRI1458
BRI1458$tissue <-"Omentum"
BRI1458$disease <-"Healthy"
BRI1458$run <-"1458"
BRI1458[["percent.mt"]] <- PercentageFeatureSet(BRI1458, pattern = "^MT-")
VlnPlot(BRI1458, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(BRI1458, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(BRI1458, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
BRI1458 <- subset(BRI1458, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 19 & percent.mt > 1)
counts_BRI1458 <- GetAssayData(BRI1458, assay = "RNA", slot = "counts")

BRI1459.data <- Read10X(data.dir = "/Users/willtrim/Library/CloudStorage/OneDrive-TCDUD.onmicrosoft.com/Post doc Trinity Harvard/Data/Two Depots Study/Data/scRNAseq - BTM/Data/GRCh38/BRI-1458/outs/raw_feature_bc_matrix")
BRI1459 <- CreateSeuratObject(counts = BRI1459.data, project = "SQfat10x", min.cells = 3, min.features = 200)
BRI1459
BRI1459$tissue <-"Omentum"
BRI1459$disease <-"Unhealthy"
BRI1459$run <-"1459"
BRI1459[["percent.mt"]] <- PercentageFeatureSet(BRI1459, pattern = "^MT-")
VlnPlot(BRI1459, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(BRI1459, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(BRI1459, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
BRI1459 <- subset(BRI1459, subset = nFeature_RNA > 600 & nFeature_RNA < 5750 & percent.mt < 18 & percent.mt > 2)
counts_BRI1459 <- GetAssayData(BRI1459, assay = "RNA", slot = "counts")

BRI1456counts <- CreateSeuratObject(counts = counts_BRI1456)
BRI1457counts <- CreateSeuratObject(counts = counts_BRI1457)
BRI1458counts <- CreateSeuratObject(counts = counts_BRI1458)
BRI1459counts <- CreateSeuratObject(counts = counts_BRI1459)
SQ_merge <- merge(BRI1456counts, y = c(BRI1457counts), merge.data = FALSE)
OM_merge <- merge(BRI1458counts, y = c(BRI1459counts), merge.data = FALSE)

SQ_merge@meta.data$disease <- c(rep("MHO", ncol(BRI1456counts)), rep("MUO", ncol(BRI1457counts)))
OM_merge@meta.data$disease <- c(rep("MHO", ncol(BRI1458counts)), rep("MUO", ncol(BRI1459counts)))
SQ_merge@meta.data$Sample <- c(rep("SQ_MHO", ncol(BRI1456)), rep("SQ_MUO", ncol(BRI1457)))
OM_merge@meta.data$Sample <- c(rep("OM_MHO", ncol(BRI1458)), rep("OM_MUO", ncol(BRI1459)))

SQ_merge <- NormalizeData(SQ_merge, normalization.method = "LogNormalize", scale.factor = 10000)
OM_merge <- NormalizeData(OM_merge, normalization.method = "LogNormalize", scale.factor = 10000)

SQ_merge <- FindVariableFeatures(SQ_merge, selection.method = "vst", nfeatures = 3000)
OM_merge <- FindVariableFeatures(OM_merge, selection.method = "vst", nfeatures = 3000)

SQ_merge[["percent.mt"]] <- PercentageFeatureSet(SQ_merge, pattern = "^MT-")
OM_merge[["percent.mt"]] <- PercentageFeatureSet(OM_merge, pattern = "^MT-")
SQ_top10 <- head(VariableFeatures(SQ_merge), 10)
OM_top10 <- head(VariableFeatures(OM_merge), 10)

SQ_merge <- ScaleData(SQ_merge, features = VariableFeatures(SQ_merge))
SQ_merge <- RunPCA(SQ_merge, features = VariableFeatures(SQ_merge), npcs = 100)
DimPlot(SQ_merge, reduction = "pca", split.by= "Sample")

OM_merge <- ScaleData(OM_merge, features = VariableFeatures(OM_merge))
OM_merge <- RunPCA(OM_merge, features = VariableFeatures(OM_merge), npcs = 100)
DimPlot(OM_merge, reduction = "pca", split.by= "Sample")

ElbowPlot(SQ_merge, ndims=50)
ElbowPlot(OM_merge, ndims=50)
SQ_merge <- RunUMAP(SQ_merge, reduction = "pca", dims = 1:50)
OM_merge <- RunUMAP(OM_merge, reduction = "pca", dims = 1:50)

DimPlot(SQ_merge, reduction = "umap", group.by = "Sample")
DimPlot(OM_merge, reduction = "umap", group.by = "Sample")

SQ_merge <- FindNeighbors(SQ_merge, dims = 1:12)
SQ_merge <- FindClusters(SQ_merge, resolution = 0.5)
DimPlot(SQ_merge, reduction = "umap")

OM_merge <- FindNeighbors(OM_merge, dims = 1:15)
OM_merge <- FindClusters(OM_merge, resolution = 0.5)
DimPlot(OM_merge, reduction = "umap")

SQ_merge <- RunUMAP(SQ_merge, dims = 1:12)
OM_merge <- RunUMAP(OM_merge, dims = 1:15)

SQ.markers <- FindAllMarkers(SQ_merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SQ.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

OM.markers <- FindAllMarkers(OM_merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
OM.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

SQ.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> SQ_top10
DoHeatmap(SQ_merge, features = SQ_top10$gene) + NoLegend()

OM.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> OM_top10
DoHeatmap(OM_merge, features = OM_top10$gene) + NoLegend()

SQ_IDs <- c('PDGFRA+1', 'PDGFRA+2', 'PDGFRA+3', 'T cells 1' , 'Myeloid cells 1' , 'NK cells', 'ECs' , 'PI16+' , 'SMCs', 'T cells 2', 'Myeloid cells 2', 'B cells')
names(SQ_IDs) <- levels(SQ_merge)
SQ_merge <- RenameIdents(SQ_merge, SQ_IDs)
DimPlot(SQ_merge, reduction = "umap", pt.size = 0.5) + NoLegend()

SQ_IDs <- c('PDGFRA+', 'PDGFRA+', 'PDGFRA+', 'T cells' , 'Myeloid cells' , 'NK cells', 'ECs' , 'PI16+' , 'SMCs', 'T cells', 'Myeloid cells', 'B cells')
names(SQ_IDs) <- levels(SQ_merge)
SQ_merge <- RenameIdents(SQ_merge, SQ_IDs)
DimPlot(SQ_merge, reduction = "umap", pt.size = 0.5)

OM_IDs <- c('T cells', 'PDGFRA+ 1', 'MSLN+ 1', 'PI16+' , 'PDGFRA+ 2' , 'MSLN+ 2', 'Myeloid cells' , 'NK cells' , 'PDGFRA+ 3', 'ECs', 'SMCs', 'B cells', 'PDGFRA+ 4', 'MSLN+ 3', 'Plasma cells')
names(OM_IDs) <- levels(OM_merge)
OM_merge <- RenameIdents(OM_merge, OM_IDs)
DimPlot(OM_merge, reduction = "umap", pt.size = 0.5) + NoLegend()

OM_IDs <- c('T cells', 'PDGFRA+', 'MSLN+', 'PI16+' , 'PDGFRA+' , 'MSLN+', 'Myeloid cells' , 'NK cells' , 'PDGFRA+', 'ECs', 'SMCs', 'B cells', 'PDGFRA+', 'MSLN+', 'Plasma cells')
names(OM_IDs) <- levels(OM_merge)
OM_merge <- RenameIdents(OM_merge, OM_IDs)
DimPlot(OM_merge, reduction = "umap", pt.size = 0.5)

write_xlsx(SQ.markers, "/Users/willtrim/Desktop/OM.markers.xlsx")
write_xlsx(SQ.markers, "/Users/willtrim/Desktop/SQ.markers.xlsx")


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

SQ_Lympho <- subset(x = SQ_merge, idents = c('T cells', 'B cells', 'NK cells'))

SQ_Lympho <- ScaleData(SQ_Lympho, features = VariableFeatures(SQ_Lympho))
SQ_Lympho <- RunPCA(SQ_Lympho, features = VariableFeatures(SQ_Lympho), npcs = 100)
SQ_Lympho <- RunUMAP(SQ_Lympho, reduction = "pca", dims = 1:50)
ElbowPlot(SQ_Lympho, ndims=50)

SQ_Lympho <- FindNeighbors(SQ_Lympho, dims = 1:8)
SQ_Lympho <- FindClusters(SQ_Lympho, resolution = 0.3)
SQ_Lympho <- RunUMAP(SQ_Lympho, dims = 1:8)

DimPlot(SQ_Lympho)

SQ_Lympho_markers <- FindAllMarkers(SQ_Lympho, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SQ_Lympho_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

SQ_Lympho_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> SQ_Lympho_top10
DoHeatmap(SQ_Lympho, features = SQ_Lympho_top10$gene) + NoLegend()



##Extract OM Immune cells
OM_Myeloid <- subset(x = OM_merge, idents = c('Myeloid cells'))


OM_Myeloid <- subset(x = OM_Myeloid, idents = c("0", "1", "2", "3", "4", "5", "6"))

OM_Myeloid <- subset(x = OM_Myeloid, idents = c("0", "2", "3", "4", "5", "7", "8"))
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


OM_Lympho <- subset(x = OM_merge, idents = c('T cells', 'B cells', 'Plasma cells', 'NK cells'))
OM_Lympho <- ScaleData(OM_Lympho, features = VariableFeatures(OM_Lympho))
OM_Lympho <- RunPCA(OM_Lympho, features = VariableFeatures(OM_Lympho), npcs = 100)
OM_Lympho <- RunUMAP(OM_Lympho, reduction = "pca", dims = 1:50)
OM_Lympho <- FindNeighbors(OM_Lympho, dims = 1:10)
OM_Lympho <- FindClusters(OM_Lympho, resolution = 0.5)
OM_Lympho <- RunUMAP(OM_Lympho, dims = 1:10)

DimPlot(OM_Lympho, split.by = "Sample") + theme(aspect.ratio = 1)

OM_Lympho_markers <- FindAllMarkers(OM_Lympho, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
OM_Lympho_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

OM_Lympho_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> OM_Lympho_top10
DoHeatmap(OM_Lympho, features = OM_Lympho_top10$gene) + NoLegend()

##DEG lists
OM_T <- subset(OM_merge, idents = c("T cells"))
DimPlot(OM_T, label = T)
OM_T_DEGs <- FindMarkers(OM_T, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_T_DEGs <- subset(OM_T_DEGs, p_val_adj < 0.05)
write.table(OM_T_DEGs, file="OM_T_DEGs.csv", sep=",")

OM_NK <- subset(OM_merge, idents = c("NK cells"))
DimPlot(OM_NK, label = T)
OM_NK_DEGs <- FindMarkers(OM_T, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_NK_DEGs <- subset(OM_NK_DEGs, p_val_adj < 0.05)
write.table(OM_NK_DEGs, file="OM_NK_DEGs.csv", sep=",")

OM_Myeloid <- subset(OM_merge, idents = c("Myeloid cells"))
DimPlot(OM_Myeloid, label = T)
OM_Myeloid_DEGs <- FindMarkers(OM_Myeloid, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_Myeloid_DEGs <- subset(OM_Myeloid_DEGs, p_val_adj < 0.05)
write.table(OM_Myeloid_DEGs, file="OM_Myeloid_DEGs.csv", sep=",")

OM_B <- subset(OM_merge, idents = c("B cells"))
DimPlot(OM_B, label = T)
OM_B_DEGs <- FindMarkers(OM_B, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_B_DEGs <- subset(OM_B_DEGs, p_val_adj < 0.05)
write.table(OM_B_DEGs, file="OM_B_DEGs.csv", sep=",")

OM_PC <- subset(OM_merge, idents = c("?"))
DimPlot(OM_PC, label = T)
OM_PC_DEGs <- FindMarkers(OM_PC, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_PC_DEGs <- subset(OM_PC_DEGs, p_val_adj < 0.05)
write.table(OM_PC_DEGs, file="OM_PlasmaCells_DEGs.csv", sep=",")

OM_PI16 <- subset(OM_merge, idents = c("PI16+"))
DimPlot(OM_PI16, label = T)
OM_PI16_DEGs <- FindMarkers(OM_PI16, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_PI16_DEGs <- subset(OM_PI16_DEGs, p_val_adj < 0.05)
write.table(OM_PI16_DEGs, file="OM_PI16_DEGs.csv", sep=",")

OM_PDGFRA <- subset(OM_merge, idents = c("PDGFRA+"))
DimPlot(OM_PDGFRA, label = T)
OM_PDGFRA_DEGs <- FindMarkers(OM_PDGFRA, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_PDGFRA_DEGs <- subset(OM_PDGFRA_DEGs, p_val_adj < 0.05)
write.table(OM_PDGFRA_DEGs, file="OM_PDGFRA_DEGs.csv", sep=",")

OM_MSLN <- subset(OM_merge, idents = c("MSLN+"))
DimPlot(OM_MSLN, label = T)
OM_MSLN_DEGs <- FindMarkers(OM_MSLN, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_MSLN_DEGs <- subset(OM_MSLN_DEGs, p_val_adj < 0.05)
write.table(OM_MSLN_DEGs, file="OM_MSLN_DEGs.csv", sep=",")

OM_EC <- subset(OM_merge, idents = c("ECs"))
DimPlot(OM_EC, label = T)
OM_EC_DEGs <- FindMarkers(OM_EC, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_EC_DEGs <- subset(OM_EC_DEGs, p_val_adj < 0.05)
write.table(OM_EC_DEGs, file="OM_EC_DEGs.csv", sep=",")

OM_SMC <- subset(OM_merge, idents = c("SMCs"))
DimPlot(OM_SMC, label = T)
OM_SMC_DEGs <- FindMarkers(OM_SMC, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_SMC_DEGs <- subset(OM_SMC_DEGs, p_val_adj < 0.05)
write.table(OM_SMC_DEGs, file="OM_SMC_DEGs.csv", sep=",")

SQ_SMC <- subset(SQ_merge, idents = c("SMCs"))
DimPlot(SQ_SMC, label = T)
SQ_SMC_DEGs <- FindMarkers(SQ_SMC, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_SMC_DEGs <- subset(SQ_SMC_DEGs, p_val_adj < 0.05)
write.table(SQ_SMC_DEGs, file="SQ_SMC_DEGs.csv", sep=",")

SQ_ECs <- subset(SQ_merge, idents = c("ECs"))
DimPlot(SQ_ECs, label = T)
SQ_ECs_DEGs <- FindMarkers(SQ_ECs, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_ECs_DEGs <- subset(SQ_ECs_DEGs, p_val_adj < 0.05)
write.table(SQ_ECs_DEGs, file="SQ_ECs_DEGs.csv", sep=",")

SQ_PDGFRA <- subset(SQ_merge, idents = c("PDGFRA+"))
DimPlot(SQ_PDGFRA, label = T)
SQ_PDGFRA_DEGs <- FindMarkers(SQ_PDGFRA, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_PDGFRA_DEGs <- subset(SQ_PDGFRA_DEGs, p_val_adj < 0.05)
write.table(SQ_PDGFRA_DEGs, file="SQ_PDGFRA_DEGs.csv", sep=",")

SQ_PI16 <- subset(SQ_merge, idents = c("PI16+"))
DimPlot(SQ_PI16, label = T)
SQ_PI16_DEGs <- FindMarkers(SQ_PI16, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_PI16_DEGs <- subset(SQ_PI16_DEGs, p_val_adj < 0.05)
write.table(SQ_PI16_DEGs, file="SQ_PI16_DEGs.csv", sep=",")

SQ_T <- subset(SQ_merge, idents = c("T cells"))
DimPlot(SQ_T, label = T)
SQ_T_DEGs <- FindMarkers(SQ_T, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_T_DEGs <- subset(SQ_T_DEGs, p_val_adj < 0.05)
write.table(SQ_T_DEGs, file="SQ_T_DEGs.csv", sep=",")

SQ_Myeloid <- subset(SQ_merge, idents = c("Myeloid cells"))
DimPlot(SQ_Myeloid, label = T)
SQ_Myeloid_DEGs <- FindMarkers(SQ_Myeloid, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_Myeloid_DEGs <- subset(SQ_Myeloid_DEGs, p_val_adj < 0.05)
write.table(SQ_Myeloid_DEGs, file="SQ_Myeloid_DEGs.csv", sep=",")

SQ_NK <- subset(SQ_merge, idents = c("NK cells"))
DimPlot(SQ_NK, label = T)
SQ_NK_DEGs <- FindMarkers(SQ_NK, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_NK_DEGs <- subset(SQ_NK_DEGs, p_val_adj < 0.05)
write.table(SQ_NK_DEGs, file="SQ_NK_DEGs.csv", sep=",")

SQ_B <- subset(SQ_merge, idents = c("B cells"))
DimPlot(SQ_B, label = T)
SQ_B_DEGs <- FindMarkers(SQ_B, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_B_DEGs <- subset(SQ_B_DEGs, p_val_adj < 0.05)
write.table(SQ_B_DEGs, file="SQ_B_DEGs.csv", sep=",")

##OM Lymphocyte Processing
OM_Lympho <- subset(x = OM_Immune, idents = c("0", "1", "2", "5", "6", "7"))
OM_Lympho <- ScaleData(OM_Lympho, features = VariableFeatures(OM_Lympho))
OM_Lympho <- RunPCA(OM_Lympho, features = VariableFeatures(OM_Lympho), npcs = 100)
OM_Lympho <- RunUMAP(OM_Lympho, reduction = "pca", dims = 1:50)
OM_Lympho <- FindNeighbors(OM_Lympho, dims = 1:10)
OM_Lympho <- FindClusters(OM_Lympho, resolution = 0.5)
OM_Lympho <- RunUMAP(OM_Lympho, dims = 1:10)

DimPlot(OM_Lympho, split.by = "Sample") + theme(aspect.ratio = 1)

OM_Lympho_IDs <- c("CD4+ T", "CD8+ T", "mNK", "NK-like", "B", "NA CD4+ T", "GD T", "ILC-like", "Plasma cells")
names(OM_Lympho_IDs) <- levels(OM_Lympho)
OM_Lympho <- RenameIdents(OM_Lympho, OM_Lympho_IDs)
DimPlot(OM_Lympho, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

OM_Lympho_markers <- FindAllMarkers(OM_Lympho, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
OM_Lympho_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.table(OM_Lympho_markers, file="OM_Lympho_markers.csv", sep=",")


##SQ Lymphocyte Processing

SQ_Lympho_IDs <- c("CD4+ T", "NK-like", "CD8+ T", "mNK", "Treg", "B")
names(SQ_Lympho_IDs) <- levels(SQ_Lympho)
SQ_Lympho <- RenameIdents(SQ_Lympho, SQ_Lympho_IDs)
DimPlot(SQ_Lympho, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

SQ_Lympho_markers <- FindAllMarkers(SQ_Lympho, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SQ_Lympho_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.table(SQ_Lympho_markers, file="SQ_Lympho_markers.csv", sep=",")

##Lymphocyte DEGs
SQ_B_cells <- subset(SQ_Lympho, idents = c("B"))
DimPlot(SQ_B_cells, label = T)
SQ_B_cells_DEGs <- FindMarkers(SQ_B_cells, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_B_cells_DEGs <- subset(SQ_B_cells_DEGs, p_val_adj < 0.05)
write.table(SQ_B_cells_DEGs, file="SQ_B_cells_DEGs.csv", sep=",")

SQ_CD4T_cells <- subset(SQ_Lympho, idents = c("CD4+ T"))
DimPlot(SQ_CD4T_cells, label = T)
SQ_CD4T_cells_DEGs <- FindMarkers(SQ_CD4T_cells, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_CD4T_cells_DEGs <- subset(SQ_CD4T_cells_DEGs, p_val_adj < 0.05)
write.table(SQ_CD4T_cells_DEGs, file="SQ_CD4T_cells_DEGs.csv", sep=",")

SQ_CD8T_cells <- subset(SQ_Lympho, idents = c("CD8+ T"))
DimPlot(SQ_CD8T_cells, label = T)
SQ_CD8T_cells_DEGs <- FindMarkers(SQ_CD8T_cells, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_CD8T_cells_DEGs <- subset(SQ_CD8T_cells_DEGs, p_val_adj < 0.05)
write.table(SQ_CD8T_cells_DEGs, file="SQ_CD8T_cells_DEGs.csv", sep=",")

SQ_Treg <- subset(SQ_Lympho, idents = c("Treg"))
DimPlot(SQ_Treg, label = T)
SQ_Treg_DEGs <- FindMarkers(SQ_Treg, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_Treg_DEGs <- subset(SQ_Treg_DEGs, p_val_adj < 0.05)
write.table(SQ_Treg_DEGs, file="SQ_Treg_DEGs.csv", sep=",")

SQ_mNK <- subset(SQ_Lympho, idents = c("mNK"))
DimPlot(SQ_mNK, label = T)
SQ_mNK_DEGs <- FindMarkers(SQ_mNK, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_mNK_DEGs <- subset(SQ_mNK_DEGs, p_val_adj < 0.05)
write.table(SQ_mNK_DEGs, file="SQ_mNK_DEGs.csv", sep=",")

SQ_NK_like <- subset(SQ_Lympho, idents = c("NK-like"))
DimPlot(SQ_NK_like, label = T)
SQ_NK_like_DEGs <- FindMarkers(SQ_NK_like, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_NK_like_DEGs <- subset(SQ_NK_like_DEGs, p_val_adj < 0.05)
write.table(SQ_NK_like_DEGs, file="SQ_NK_like_DEGs.csv", sep=",")

OM_NK_like <- subset(OM_Lympho, idents = c("NK-like"))
DimPlot(OM_NK_like, label = T)
OM_NK_like_DEGs <- FindMarkers(OM_NK_like, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_NK_like_DEGs <- subset(OM_NK_like_DEGs, p_val_adj < 0.05)
write.table(OM_NK_like_DEGs, file="OM_NK_like_DEGs.csv", sep=",")

OM_CD4 <- subset(OM_Lympho, idents = c("CD4+ T"))
DimPlot(OM_CD4, label = T)
OM_CD4_DEGs <- FindMarkers(OM_CD4, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_CD4_DEGs <- subset(OM_CD4_DEGs, p_val_adj < 0.05)
write.table(OM_CD4_DEGs, file="OM_CD4_DEGs.csv", sep=",")

OM_CD8 <- subset(OM_Lympho, idents = c("CD8+ T"))
DimPlot(OM_CD8, label = T)
OM_CD8_DEGs <- FindMarkers(OM_CD8, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_CD8_DEGs <- subset(OM_CD8_DEGs, p_val_adj < 0.05)
write.table(OM_CD8_DEGs, file="OM_CD8_DEGs.csv", sep=",")

OM_mNK <- subset(OM_Lympho, idents = c("mNK"))
DimPlot(OM_mNK, label = T)
OM_mNK_DEGs <- FindMarkers(OM_mNK, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_mNK_DEGs <- subset(OM_mNK_DEGs, p_val_adj < 0.05)
write.table(OM_mNK_DEGs, file="OM_mNK_DEGs.csv", sep=",")

OM_B_cells <- subset(OM_Lympho, idents = c("B"))
DimPlot(OM_B_cells, label = T)
OM_B_cells_DEGs <- FindMarkers(OM_B_cells, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_B_cells_DEGs <- subset(OM_B_cells_DEGs, p_val_adj < 0.05)
write.table(OM_B_cells_DEGs, file="OM_B_cells_DEGs.csv", sep=",")

OM_NACD4 <- subset(OM_Lympho, idents = c("NA CD4+ T"))
DimPlot(OM_NACD4, label = T)
OM_NACD4_DEGs <- FindMarkers(OM_NACD4, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_NACD4_DEGs <- subset(OM_NACD4_DEGs, p_val_adj < 0.05)
write.table(OM_NACD4_DEGs, file="OM_NACD4_DEGs.csv", sep=",")

OM_GD <- subset(OM_Lympho, idents = c("GD T"))
DimPlot(OM_GD, label = T)
OM_GD_DEGs <- FindMarkers(OM_GD, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_GD_DEGs <- subset(OM_GD_DEGs, p_val_adj < 0.05)
write.table(OM_GD_DEGs, file="OM_GD_DEGs.csv", sep=",")

OM_ILC <- subset(OM_Lympho, idents = c("ILC-like"))
DimPlot(OM_ILC, label = T)
OM_ILC_DEGs <- FindMarkers(OM_ILC, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_ILC_DEGs <- subset(OM_ILC_DEGs, p_val_adj < 0.05)
write.table(OM_ILC_DEGs, file="OM_ILC_DEGs.csv", sep=",")

OM_Plasma_cells <- subset(OM_Lympho, idents = c("Plasma cells"))
DimPlot(OM_Plasma_cells, label = T)
OM_Plasma_cells_DEGs <- FindMarkers(OM_Plasma_cells, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_Plasma_cells_DEGs <- subset(OM_Plasma_cells_DEGs, p_val_adj < 0.05)
write.table(OM_Plasma_cells_DEGs, file="OM_Plasma_cells_DEGs.csv", sep=",")


##Myeloid clustering
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

OM_Myeloid_markers <- FindAllMarkers(OM_Myeloid, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
OM_Myeloid_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

OM_Myeloid_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> OM_Myeloid_top10
DoHeatmap(OM_Myeloid, features = OM_Myeloid_top10$gene) + NoLegend()

write.table(OM_Myeloid_markers, file="OM_Myeloid_markers.csv", sep=",")

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

##Myeloid DEGs
OM_THBS1 <- subset(OM_Myeloid, idents = c("THBS1+"))
OM_THBS1_DEGs <- FindMarkers(OM_THBS1, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_THBS1_DEGs <- subset(OM_THBS1_DEGs, p_val_adj < 0.05)
write.table(OM_THBS1_DEGs, file="OM_THBS1_DEGs.csv", sep=",")

OM_LAMs <- subset(OM_Myeloid, idents = c("LAMs"))
OM_LAMs_DEGs <- FindMarkers(OM_LAMs, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_LAMs_DEGs <- subset(OM_LAMs_DEGs, p_val_adj < 0.05)
write.table(OM_LAMs_DEGs, file="OM_LAMs_DEGs.csv", sep=",")

OM_cDC2B <- subset(OM_Myeloid, idents = c("cDC2Bs"))
OM_cDC2B_DEGs <- FindMarkers(OM_cDC2B, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_cDC2B_DEGs <- subset(OM_cDC2B_DEGs, p_val_adj < 0.05)
write.table(OM_cDC2B_DEGs, file="OM_cDC2B_DEGs.csv", sep=",")

OM_Monocytes <- subset(OM_Myeloid, idents = c("Monocytes"))
OM_Monocytes_DEGs <- FindMarkers(OM_Monocytes, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_Monocytes_DEGs <- subset(OM_Monocytes_DEGs, p_val_adj < 0.05)
write.table(OM_Monocytes_DEGs, file="OM_Monocytes_DEGs.csv", sep=",")

OM_PVMs <- subset(OM_Myeloid, idents = c("PVMs"))
OM_PVMs_DEGs <- FindMarkers(OM_PVMs, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_PVMs_DEGs <- subset(OM_PVMs_DEGs, p_val_adj < 0.05)
write.table(OM_PVMs_DEGs, file="OM_PVMs_DEGs.csv", sep=",")

OM_FBLN1 <- subset(OM_Myeloid, idents = c("FBLN1+"))
OM_FBLN1_DEGs <- FindMarkers(OM_FBLN1, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_FBLN1_DEGs <- subset(OM_FBLN1_DEGs, p_val_adj < 0.05)
write.table(OM_FBLN1_DEGs, file="OM_FBLN1_DEGs.csv", sep=",")

OM_cDC2A <- subset(OM_Myeloid, idents = c("cDC2As"))
OM_cDC2A_DEGs <- FindMarkers(OM_cDC2A, group.by = "Sample", ident.1 = "OM_MHO", ident.2 = "OM_MUO")
OM_cDC2A_DEGs <- subset(OM_cDC2A_DEGs, p_val_adj < 0.05)
write.table(OM_cDC2A_DEGs, file="OM_cDC2A_DEGs.csv", sep=",")

SQ_LAMs <- subset(SQ_Myeloid, idents = c("LAMs"))
SQ_LAMs_DEGs <- FindMarkers(SQ_LAMs, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_LAMs_DEGs <- subset(SQ_LAMs_DEGs, p_val_adj < 0.05)
write.table(SQ_LAMs_DEGs, file="SQ_LAMs_DEGs.csv", sep=",")

SQ_cDC2Bs <- subset(SQ_Myeloid, idents = c("cDC2Bs"))
SQ_cDC2Bs_DEGs <- FindMarkers(SQ_cDC2Bs, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_cDC2Bs_DEGs <- subset(SQ_cDC2Bs_DEGs, p_val_adj < 0.05)
write.table(SQ_cDC2Bs_DEGs, file="SQ_cDC2Bs_DEGs.csv", sep=",")

SQ_Monocytes <- subset(SQ_Myeloid, idents = c("Monocytes"))
SQ_Monocytes_DEGs <- FindMarkers(SQ_Monocytes, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_Monocytes_DEGs <- subset(SQ_Monocytes_DEGs, p_val_adj < 0.05)
write.table(SQ_Monocytes_DEGs, file="SQ_Monocytes_DEGs.csv", sep=",")

SQ_PVMs <- subset(SQ_Myeloid, idents = c("PVMs"))
SQ_PVMs_DEGs <- FindMarkers(SQ_PVMs, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_PVMs_DEGs <- subset(SQ_PVMs_DEGs, p_val_adj < 0.05)
write.table(SQ_PVMs_DEGs, file="SQ_PVMs_DEGs.csv", sep=",")

SQ_cDC1 <- subset(SQ_Myeloid, idents = c("cDC1s"))
SQ_cDC1_DEGs <- FindMarkers(SQ_cDC1, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_cDC1_DEGs <- subset(SQ_cDC1_DEGs, p_val_adj < 0.05)
write.table(SQ_cDC1_DEGs, file="SQ_cDC1_DEGs.csv", sep=",")

SQ_cDC2A <- subset(SQ_Myeloid, idents = c("cDC2As"))
SQ_cDC2A_DEGs <- FindMarkers(SQ_cDC2A, group.by = "Sample", ident.1 = "SQ_MHO", ident.2 = "SQ_MUO")
SQ_cDC2A_DEGs <- subset(SQ_cDC2A_DEGs, p_val_adj < 0.05)
write.table(SQ_cDC2A_DEGs, file="SQ_cDC2A_DEGs.csv", sep=",")