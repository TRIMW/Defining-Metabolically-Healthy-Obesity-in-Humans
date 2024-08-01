# SplitDotPlotGG has been replaced with the `split.by` parameter for DotPlot
DotPlot(myeloid.cells, features = “THBS1”, “FCGR3A, split.by = "groups") + RotatedAxis()

##how to partition samples to run DE between MHO and mUO for each population
# step 1: make copy of seurat object to manipulate. 
myeloid.cells.copy <- myeloid.cells
#Step 2: Split identity by tissue. 
Idents(myeloid.cells.copy) <- "Tissue"
DimPlot(myeloid.cells.copy)
#Step 3: subset out SQ (or OM), make seurat clusters now for new tissue specific file. 
SQ_myeloid.cells.copy <- subset(myeloid.cells.copy, idents = c("SQ"))
Idents(SQ_myeloid.cells.copy) <- "seurat_clusters"
DimPlot(SQ_myeloid.cells.copy)
# Step 4: Partition out each respective cell cluster into new seurat object for isolated DEG identification
myeloid_Cluster0 <- subset(SQ_myeloid.cells.copy, idents = "0")
#Step 5: find DE markers between conditions
Myeloid_Markers_Cluster_0_SQ_by_disease <- FindMarkers(myeloid_Cluster0, group.by = "Disease", ident.1 = "MHO", ident.2 = "MUO")
#Step 6: Isolate only P<0.05 DEGs, overwriting original file to exclude genes >0.05
Myeloid_Markers_Cluster_0_SQ_by_disease <- subset(Myeloid_Markers_Cluster_0_SQ_by_disease, p_val_adj < 0.05)
#Now iterate through each population (Steps 4–6) to get DEGs for each population