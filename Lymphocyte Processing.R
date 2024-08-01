##Rename clusters and re-order
OM_Lympho_IDs <- c('CD4+ T', 'CD8+ T', 'mNKs', 'NK-like' , 'B cells' , 'NA CD4+ T', 'GD T' , 'ILC-like' , 'Plasma cells')
names(OM_Lympho_IDs) <- levels(OM_Lymph)
OM_Lymph <- RenameIdents(OM_Lymph, OM_Lympho_IDs)
DimPlot(OM_Lymph, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

OM_Lympho_levels <- c('CD4+ T', 'NA CD4+ T', 'CD8+ T', 'GD T' , 'mNKs', 'NK-like' , 'B cells' ,'Plasma cells', 'ILC-like')
Idents(OM_Lymph) <- factor(Idents(OM_Lymph), levels= OM_Lympho_levels)

##Set colour scheme
OM_Lympho_colors <- c('CD4+ T' = '#005f73', 'CD8+ T' = 'deeppink', 'mNKs' = 'peru', 'NK-like' = 'chartreuse3', 'B cells' = 'skyblue2', 'NA CD4+ T' = 'gold', 'GD T' = 'purple', 'ILC-like' = 'dodgerblue', 'Plasma cells' = '#e9d8a6')

SCpubr::do_DimPlot(sample = OM_Lymph,
                   colors.use = OM_Lympho_colors, plot_cell_borders = TRUE, border.size = 2, border.color = 'black', pt.size = 0.6, label = TRUE, repel = TRUE)
vln_plot <- VlnPlot(OM_Lymph, 
                    features = c("EML4", 'JCHAIN', "CD79A", "FCER1G", "FCGR3A", "TRGC2", "CD8A",  "LEF1", "IL7R", "CD3E"),
                    stack = TRUE, 
                    flip = TRUE, 
                    cols = c("EML4" = 'dodgerblue', "IL7R" = '#005f73', "LEF1" = 'gold',"CD8A" = 'deeppink',"TRGC2" = 'purple', "FCGR3A" = 'peru', "FCER1G" = 'chartreuse3', "CD79A" = 'skyblue2', 'JCHAIN' = '#e9d8a6', "CD3E" = 'ivory' ))

vln_plot <- vln_plot + theme(legend.position = "none")
print(vln_plot)

SCpubr::do_DimPlot(sample = OM_Lymph,
                   colors.use = c('OM_MHO' = 'grey', 'OM_MUO' = 'gold'), group.by = 'Sample',
                   plot_cell_borders = TRUE, border.size = 2, border.color = 'black', pt.size = 0.6, label = FALSE, repel = TRUE)

##Rename clusters and re-order
SQ_Lympho_IDs <- c('CD4+ T', 'NK-like', 'CD8+ T', 'mNKs', 'Tregs', 'B cells')
names(SQ_Lympho_IDs) <- levels(SQ_Lymph)
SQ_Lymph <- RenameIdents(SQ_Lymph, SQ_Lympho_IDs)
DimPlot(SQ_Lymph, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

SQ_Lympho_levels <- c('CD4+ T', 'CD8+ T', 'Tregs', 'mNKs', 'NK-like', 'B cells')
Idents(SQ_Lymph) <- factor(Idents(SQ_Lymph), levels= SQ_Lympho_levels)


##Set colour scheme
SQ_Lympho_colors <- c('CD4+ T' = '#005f73', 'CD8+ T' = 'deeppink', 'NK-like' = 'chartreuse3', 'B cells' = 'skyblue2', 'mNKs' = 'peru', 'Tregs' = 'coral')

SCpubr::do_DimPlot(sample = SQ_Lymph,
                   colors.use = SQ_Lympho_colors, plot_cell_borders = TRUE, border.size = 2, border.color = 'black', pt.size = 0.6, label = TRUE, repel = TRUE)

SCpubr::do_DimPlot(sample = SQ_Lymph,
                   colors.use = c('SQ_MHO' = 'grey', 'SQ_MUO' = 'gold'), group.by = 'Sample',
                   plot_cell_borders = TRUE, border.size = 2, border.color = 'black', pt.size = 0.9, label = FALSE, repel = TRUE)

vln_plot <- VlnPlot(SQ_Lymph, 
                    features = c("JCHAIN", 'FCER1G', "KLRC1", "FOXP3", "CTLA4", "CD8A", "IL7R", "CD3E"),
                    stack = TRUE,  
                    flip = TRUE, 
                    cols = c('IL7R' = '#005f73', 'CD8A' = 'deeppink', 'FCER1G' = 'chartreuse3', 'JCHAIN' = 'skyblue2', 'KLRC1' = 'peru', 'CTLA4' = 'coral', 'FOXP3' = 'coral', 'CD3E' = 'ivory'))

vln_plot <- vln_plot + theme(legend.position = "none")
print(vln_plot)


##Split Data sets
OM_Lymph_Copy <- OM_Lymph
Idents(OM_Lymph_Copy) <- "Sample"
DimPlot(OM_Lymph_Copy)

##Isolate Healthy Obese Group
OM_Lymph_Copy_MHO <- subset(OM_Lymph_Copy, idents = c("OM_MHO"))
Idents(OM_Lymph_Copy_MHO) <- "seurat_clusters"
DimPlot(OM_Lymph_Copy_MHO)

##Isolate unhealthy obese group
OM_Lymph_MUO_Copy <- subset(OM_Lymph_Copy, idents = c("OM_MUO"))
Idents(OM_Lymph_MUO_Copy) <- "seurat_clusters"
DimPlot(OM_Lymph_MUO_Copy, label = T)

##Split Data sets
SQ_Lymph_copy <- SQ_Lymph
Idents(SQ_Lymph_copy) <- "Sample"
DimPlot(SQ_Lymph_copy)

##Isolate Healthy Obese Group
SQ_Lymph_Copy_MHO <- subset(SQ_Lymph_copy, idents = c("SQ_MHO"))
Idents(SQ_Lymph_Copy_MHO) <- "seurat_clusters"
DimPlot(SQ_Lymph_Copy_MHO)

##Isolate unhealthy obese group
SQ_Lymph_MUO_Copy <- subset(SQ_Lymph_copy, idents = c("SQ_MUO"))
Idents(SQ_Lymph_MUO_Copy) <- "seurat_clusters"
DimPlot(SQ_Lymph_MUO_Copy, label = T)


OM_Lympho_Split_IDs <- c('CD4+ T', 'CD8+ T', 'CD4+ T', 'CD8+ T', 'mNKs', 'NK-like' , 'B cells' , 'NA CD4+ T', 'GD T' , 'ILC-like' , 'Plasma cells')
names(OM_Lympho_Split_IDs) <- levels(OM_Lymph_Copy_MHO)
OM_Lymph_Copy_MHO <- RenameIdents(OM_Lymph_Copy_MHO, OM_Lympho_Split_IDs)
DimPlot(OM_Lymph_Copy_MHO, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

OM_Lympho_Split_IDs <- c('CD4+ T', 'CD8+ T', 'CD4+ T', 'CD8+ T', 'mNKs', 'NK-like' , 'B cells' , 'NA CD4+ T', 'GD T' , 'ILC-like' , 'Plasma cells')
names(OM_Lympho_Split_IDs) <- levels(OM_Lymph_MUO_Copy)
OM_Lymph_MUO_Copy <- RenameIdents(OM_Lymph_MUO_Copy, OM_Lympho_Split_IDs)
DimPlot(OM_Lymph_MUO_Copy, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

OM_Lympho_levels <- c('CD4+ T', 'NA CD4+ T', 'CD8+ T', 'GD T' , 'mNKs', 'NK-like' , 'B cells' ,'Plasma cells', 'ILC-like')
Idents(OM_Lymph_Copy_MHO) <- factor(Idents(OM_Lymph_Copy_MHO), levels= OM_Lympho_levels)

OM_Lympho_levels <- c('CD4+ T', 'NA CD4+ T', 'CD8+ T', 'GD T' , 'mNKs', 'NK-like' , 'B cells' ,'Plasma cells', 'ILC-like')
Idents(OM_Lymph_MUO_Copy) <- factor(Idents(OM_Lymph_MUO_Copy), levels= OM_Lympho_levels)


SQ_Lympho_Split_IDs <- c('CD4+ T', 'NKs', 'CD8+ T', 'GD T', 'Tregs', 'B cells', 'B cells')
names(SQ_Lympho_Split_IDs) <- levels(SQ_Lymph_Copy_MHO)
SQ_Lymph_Copy_MHO <- RenameIdents(SQ_Lymph_Copy_MHO, SQ_Lympho_Split_IDs)
DimPlot(SQ_Lymph_Copy_MHO, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

SQ_Lympho_Split_IDs <- c('CD4+ T', 'NKs', 'CD8+ T', 'GD T', 'Tregs', 'B cells', 'B cells')
names(SQ_Lympho_Split_IDs) <- levels(SQ_Lymph_MUO_Copy)
SQ_Lymph_MUO_Copy <- RenameIdents(SQ_Lymph_MUO_Copy, SQ_Lympho_Split_IDs)
DimPlot(SQ_Lymph_MUO_Copy, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


SQ_Lympho_levels <- c('CD4+ T', 'CD8+ T', 'Tregs', 'GD T' , 'NKs', 'B cells')
Idents(SQ_Lymph_Copy_MHO) <- factor(Idents(SQ_Lymph_Copy_MHO), levels= SQ_Lympho_levels)

SQ_Lympho_levels <- c('CD4+ T', 'CD8+ T', 'Tregs', 'GD T' , 'NKs', 'B cells')
Idents(SQ_Lymph_MUO_Copy) <- factor(Idents(SQ_Lymph_MUO_Copy), levels= SQ_Lympho_levels)


##Get cell counts
OM_Lymph_MHO_Cells <- OM_Lymph_Copy_MHO@active.ident %>% as.data.table
View(OM_Lymph_MHO_Cells)
write.xlsx2(OM_Lymph_MHO_Cells, 'OM_Lymph_MHO_Counts.xlsx', sheetName = "Sheet1",
            col.names = TRUE, row.names = TRUE, append = FALSE)

OM_Lymph_MUO_Cells <- OM_Lymph_MUO_Copy@active.ident %>% as.data.table
View(OM_Lymph_MUO_Cells)
write.xlsx2(OM_Lymph_MUO_Cells, 'OM_Lymph_MUO_Counts.xlsx', sheetName = "Sheet1",
            col.names = TRUE, row.names = TRUE, append = FALSE)


##Exhaustion Markers Module Score
Exhaustion_Genes <- list(c('TIGIT', 'HAVCR2', 'PDCD1', 'TOX', 'TOX2', 'CD27', 'TCF7', 'LAG3', 'KLRG1', 'CTLA4', 'CD200R1', 'TNFSF4', 'CD86', '2B4', 'CD160', 'BTLA4', 'TIM3'))
OM_Lymph <- AddModuleScore(object = OM_Lymph, features = Exhaustion_Genes, name = 'Exhaustion_Genes')
SQ_Lymph <- AddModuleScore(object = SQ_Lymph, features = Exhaustion_Genes, name = 'Exhaustion_Genes')

OM_Lymph_Copy_MHO <- AddModuleScore(object = OM_Lymph_Copy_MHO, features = Exhaustion_Genes, name = 'Exhaustion_Genes')
OM_Lymph_MUO_Copy <- AddModuleScore(object = OM_Lymph_MUO_Copy, features = Exhaustion_Genes, name = 'Exhaustion_Genes')


SQ_Lymph_Copy_MHO <- AddModuleScore(object = SQ_Lymph_Copy_MHO, features = Exhaustion_Genes, name = 'Exhaustion_Genes')
SQ_Lymph_MUO_Copy <- AddModuleScore(object = SQ_Lymph_MUO_Copy, features = Exhaustion_Genes, name = 'Exhaustion_Genes')

Metal_Ion_Genes <- list(c('MT1A', 'MT1B', 'MT1CP', 'MT1DP', 'MT1E', 'MT1F', 'MT1G', 'MT1H', 'MT1HL1', 'MT1IP', 'MT1JP', 'MT1L', 'MT1M', 'MT1P1', 'MT1P3', 'MT1X', 'MT2A', 'MT3', 'MT4', 'CP', 'DMT1', 'FPN1'))
OM_Lymph <- AddModuleScore(object = OM_Lymph, features = Metal_Ion_Genes, name = 'Metal_Ion_Genes')
SQ_Lymph <- AddModuleScore(object = SQ_Lymph, features = Metal_Ion_Genes, name = 'Metal_Ion_Genes')

OM_Lymph_Copy_MHO <- AddModuleScore(object = OM_Lymph_Copy_MHO, features = Metal_Ion_Genes, name = 'Metal_Ion_Genes')
OM_Lymph_MUO_Copy <- AddModuleScore(object = OM_Lymph_MUO_Copy, features = Metal_Ion_Genes, name = 'Metal_Ion_Genes')

SQ_Lymph_Copy_MHO <- AddModuleScore(object = SQ_Lymph_Copy_MHO, features = Metal_Ion_Genes, name = 'Metal_Ion_Genes')
SQ_Lymph_MUO_Copy <- AddModuleScore(object = SQ_Lymph_MUO_Copy, features = Metal_Ion_Genes, name = 'Metal_Ion_Genes')


###HALLMARK Processing
OM_B_MHO <- seurat_extract(OM_B,
                           meta1 = "Sample", value_meta1 = "OM_MHO")
OM_B_MUO <- seurat_extract(OM_B,
                           meta1 = "Sample", value_meta1 = "OM_MUO")
OM_B_Comparison <- compare_pathways(samples = list(OM_B_MHO, OM_B_MUO), pathways = pathways)

OM_CD4_MHO <- seurat_extract(OM_CD4, meta1 = "Sample", value_meta1 = "OM_MHO")
OM_CD4_MUO <- seurat_extract(OM_CD4, meta1 = "Sample", value_meta1 = "OM_MUO")
OM_CD4_Comparison <- compare_pathways(samples = list(OM_CD4_MHO, OM_CD4_MUO), pathways = pathways)

OM_CD8_MHO <- seurat_extract(OM_CD8, meta1 = "Sample", value_meta1 = "OM_MHO")
OM_CD8_MUO <- seurat_extract(OM_CD8, meta1 = "Sample", value_meta1 = "OM_MUO")
OM_CD8_Comparison <- compare_pathways(samples = list(OM_CD8_MHO, OM_CD8_MUO), pathways = pathways)

OM_ILCs_MHO <- seurat_extract(OM_ILCs, meta1 = "Sample", value_meta1 = "OM_MHO")
OM_ILCs_MUO <- seurat_extract(OM_ILCs, meta1 = "Sample", value_meta1 = "OM_MUO")
OM_ILCs_Comparison <- compare_pathways(samples = list(OM_ILCs_MHO, OM_ILCs_MUO), pathways = pathways)

OM_Lymph_GDT_MHO <- seurat_extract(OM_Lymph_GDT, meta1 = "Sample", value_meta1 = "OM_MHO")
OM_Lymph_GDT_MUO <- seurat_extract(OM_Lymph_GDT, meta1 = "Sample", value_meta1 = "OM_MUO")
OM_GDT_Comparison <- compare_pathways(samples = list(OM_Lymph_GDT_MHO, OM_Lymph_GDT_MUO), pathways = pathways)

OM_Lymph_NACD4T_MHO <- seurat_extract(OM_Lymph_NACD4T, meta1 = "Sample", value_meta1 = "OM_MHO")
OM_Lymph_NACD4T_MUO <- seurat_extract(OM_Lymph_NACD4T, meta1 = "Sample", value_meta1 = "OM_MUO")
OM_NACD4T_Comparison <- compare_pathways(samples = list(OM_Lymph_NACD4T_MHO, OM_Lymph_NACD4T_MUO), pathways = pathways)

OM_Lymph_PlasmaCells_MHO <- seurat_extract(OM_Lymph_PlasmaCells, meta1 = "Sample", value_meta1 = "OM_MHO")
OM_Lymph_PlasmaCells_MUO <- seurat_extract(OM_Lymph_PlasmaCells, meta1 = "Sample", value_meta1 = "OM_MUO")
OM_PlasmaCells_Comparison <- compare_pathways(samples = list(OM_Lymph_PlasmaCells_MHO, OM_Lymph_PlasmaCells_MUO), pathways = pathways)

OM_mNK_MHO <- seurat_extract(OM_mNK, meta1 = "Sample", value_meta1 = "OM_MHO")
OM_mNK_MUO <- seurat_extract(OM_mNK, meta1 = "Sample", value_meta1 = "OM_MUO")
OM_mNK_Comparison <- compare_pathways(samples = list(OM_mNK_MHO, OM_mNK_MUO), pathways = pathways)

OM_NK_like_MHO <- seurat_extract(OM_NK_like, meta1 = "Sample", value_meta1 = "OM_MHO")
OM_NK_like_MUO <- seurat_extract(OM_NK_like, meta1 = "Sample", value_meta1 = "OM_MUO")
OM_NK_like_Comparison <- compare_pathways(samples = list(OM_NK_like_MHO, OM_NK_like_MUO), pathways = pathways)

SQ_B_MHO <- seurat_extract(SQ_B,
                           meta1 = "Sample", value_meta1 = "SQ_MHO")
SQ_B_MUO <- seurat_extract(SQ_B,
                           meta1 = "Sample", value_meta1 = "SQ_MUO")
SQ_B_Comparison <- compare_pathways(samples = list(SQ_B_MHO, SQ_B_MUO), pathways = pathways)

SQ_CD4_MHO <- seurat_extract(SQ_CD4, meta1 = "Sample", value_meta1 = "SQ_MHO")
SQ_CD4_MUO <- seurat_extract(SQ_CD4, meta1 = "Sample", value_meta1 = "SQ_MUO")
SQ_CD4_Comparison <- compare_pathways(samples = list(SQ_CD4_MHO, SQ_CD4_MUO), pathways = pathways)

SQ_CD8_MHO <- seurat_extract(SQ_CD8, meta1 = "Sample", value_meta1 = "SQ_MHO")
SQ_CD8_MUO <- seurat_extract(SQ_CD8, meta1 = "Sample", value_meta1 = "SQ_MUO")
SQ_CD8_Comparison <- compare_pathways(samples = list(SQ_CD8_MHO, SQ_CD8_MUO), pathways = pathways)

SQ_GDT_MHO <- seurat_extract(SQ_GDT, meta1 = "Sample", value_meta1 = "SQ_MHO")
SQ_GDT_MUO <- seurat_extract(SQ_GDT, meta1 = "Sample", value_meta1 = "SQ_MUO")
SQ_GDT_Comparison <- compare_pathways(samples = list(SQ_GDT_MHO, SQ_GDT_MUO), pathways = pathways)

SQ_NK_MHO <- seurat_extract(SQ_NK, meta1 = "Sample", value_meta1 = "SQ_MHO")
SQ_NK_MUO <- seurat_extract(SQ_NK, meta1 = "Sample", value_meta1 = "SQ_MUO")
SQ_NK_Comparison <- compare_pathways(samples = list(SQ_NK_MHO, SQ_NK_MUO), pathways = pathways)

SQ_Tregs_MHO <- seurat_extract(SQ_Tregs, meta1 = "Sample", value_meta1 = "SQ_MHO")
SQ_Tregs_MUO <- seurat_extract(SQ_Tregs, meta1 = "Sample", value_meta1 = "SQ_MUO")
SQ_Tregs_Comparison <- compare_pathways(samples = list(SQ_Tregs_MHO, SQ_Tregs_MUO), pathways = pathways)

write.xlsx2(OM_B_Comparison, 'OM_B_Comparison.xlsx', sheetName = "Sheet1",
            col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx2(OM_CD4_Comparison, 'OM_CD4_Comparison.xlsx', sheetName = "Sheet1",
            col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx2(OM_CD8_Comparison, 'OM_CD8_Comparison.xlsx', sheetName = "Sheet1",
            col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx2(OM_ILCs_Comparison, 'OM_ILCs_Comparison.xlsx', sheetName = "Sheet1",
            col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx2(OM_GDT_Comparison, 'OM_GDT_Comparison.xlsx', sheetName = "Sheet1",
            col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx2(OM_NACD4T_Comparison, 'OM_NACD4T_Comparison.xlsx', sheetName = "Sheet1",
            col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx2(OM_PlasmaCells_Comparison, 'OM_PlasmaCells_Comparison.xlsx', sheetName = "Sheet1",
            col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx2(OM_mNK_Comparison, 'OM_mNK_Comparison.xlsx', sheetName = "Sheet1",
            col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx2(OM_NK_like_Comparison, 'OM_NK_like_Comparison.xlsx', sheetName = "Sheet1",
            col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx2(SQ_B_Comparison, 'SQ_B_Comparison.xlsx', sheetName = "Sheet1",
            col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx2(SQ_CD4_Comparison, 'SQ_CD4_Comparison.xlsx', sheetName = "Sheet1",
            col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx2(SQ_CD8_Comparison, 'SQ_CD8_Comparison.xlsx', sheetName = "Sheet1",
            col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx2(SQ_GDT_Comparison, 'SQ_GDT_Comparison.xlsx', sheetName = "Sheet1",
            col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx2(SQ_NK_Comparison, 'SQ_NK_Comparison.xlsx', sheetName = "Sheet1",
            col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx2(SQ_Tregs_Comparison, 'SQ_Tregs_Comparison.xlsx', sheetName = "Sheet1",
            col.names = TRUE, row.names = TRUE, append = FALSE)
