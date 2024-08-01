##Splitting a file by group
SQ_merge_copy <- SQ_merge
Idents(SQ_merge_copy) <- "Sample"
DimPlot(SQ_merge_copy)
SQ_MHO_merge_copy <- subset(SQ_merge_copy, idents = c("SQ_MHO"))
Idents(SQ_MHO_merge_copy) <- "seurat_clusters"
DimPlot(SQ_MHO_merge_copy)
SQ_MUO_merge_copy <- subset(SQ_merge_copy, idents = c("SQ_MUO"))
Idents(SQ_MUO_merge_copy) <- "seurat_clusters"
DimPlot(SQ_MUO_merge_copy)
