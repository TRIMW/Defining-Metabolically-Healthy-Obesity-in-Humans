## https://msraredon.github.io/Connectome/articles/01%20Connectome%20Workflow.html

library(Seurat)
library(SeuratData)
library(Connectome)
library(ggplot2)
library(cowplot)

## Split by Disease status
OM_merge_copy <- OM_merge
Idents(OM_merge_copy) <- "Sample"
OM_MHO_merge_copy <- subset(OM_merge_copy, idents = c("OM_MHO"))
Idents(OM_MHO_merge_copy) <- "seurat_clusters"
OM_MUO_merge_copy <- subset(OM_merge_copy, idents = c("OM_MUO"))
Idents(OM_MUO_merge_copy) <- "seurat_clusters"
DimPlot(OM_MHO_merge_copy)
DimPlot(OM_MUO_merge_copy)

## Rename Clusters
OM_MHO_MUO_IDs <- c("T cells", "PDGFRA+", "T cells", "MSLN+", "PI16+", "PDGFRA+", "MSLN+", "Myeloid cells", "Myeloid cells", "NK cells", "PDGFRA+", "ECs", "Myofibroblasts", "B cells", "PDGFRA+", "MSLN+", "Plasma cells")
names(OM_MHO_MUO_IDs) <- levels(OM_MHO_merge_copy)
OM_MHO_merge_copy <- RenameIdents(OM_MHO_merge_copy, OM_MHO_MUO_IDs)
DimPlot(OM_MHO_merge_copy, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(OM_MHO_merge_copy, reduction = "umap", pt.size = 0.5, label = T) + NoLegend()

OM_MHO_MUO_IDs <- c("T cells", "PDGFRA+", "T cells", "MSLN+", "PI16+", "PDGFRA+", "MSLN+", "Myeloid cells", "Myeloid cells", "NK cells", "PDGFRA+", "ECs", "Myofibroblasts", "B cells", "PDGFRA+", "MSLN+", "Plasma cells")
names(OM_MHO_MUO_IDs) <- levels(OM_MUO_merge_copy)
OM_MUO_merge_copy <- RenameIdents(OM_MUO_merge_copy, OM_MHO_MUO_IDs)
DimPlot(OM_MUO_merge_copy, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(OM_MUO_merge_copy, reduction = "umap", pt.size = 0.5) + NoLegend()

## MHO CONNECTOME ANALYSIS PERFORMED FIRST - MUO REQUIRES COPY AND CHANGE MHO TO MUO THROUGHOUT ###

## Scale and make connectome
OM_MHO_All <- NormalizeData(OM_MHO_merge_copy)
connectome.genes <- union(Connectome::ncomms8866_human$Ligand.ApprovedSymbol,Connectome::ncomms8866_human$Receptor.ApprovedSymbol)
genes <- connectome.genes[connectome.genes %in% rownames(OM_MHO_All)]
OM_MHO_All <- ScaleData(OM_MHO_All,features = genes)
OM_MHO_All.con <- CreateConnectome(OM_MHO_All,species = 'human',min.cells.per.ident = 75,p.values = F,calculate.DOR = F)

## Filter edges
p1 <- ggplot(OM_MHO_All.con, aes(x=ligand.scale)) + geom_density() + ggtitle('Ligand.scale')
p2 <- ggplot(OM_MHO_All.con, aes(x=recept.scale)) + geom_density() + ggtitle('Recept.scale')
p3 <- ggplot(OM_MHO_All.con, aes(x=percent.target)) + geom_density() + ggtitle('Percent.target')
p4 <- ggplot(OM_MHO_All.con, aes(x=percent.source)) + geom_density() + ggtitle('Percent.source')
plot_grid(p1,p2,p3,p4)

OM_MHO_All.con2 <- FilterConnectome(OM_MHO_All.con,min.pct = 0.1,min.z = 0.25,remove.na = T)

## If you want to view a specific Receptor/Ligand, use this command
p1 <- NetworkPlot(OM_MHO_All.con2,features = 'VEGFA',min.pct = 0.1,weight.attribute = 'weight_sc',include.all.nodes = T)
p2 <- NetworkPlot(OM_MHO_All.con2,features = 'VEGFA',min.pct = 0.75,weight.attribute = 'weight_sc',include.all.nodes = T)
plot_grid(p1,p2,nrow=1)

## Centrality Analysis
Centrality(OM_MHO_All.con2,
           modes.include = NULL,
           min.z = NULL,
           weight.attribute = 'weight_sc',
           group.by = 'mode', cols.use = OM_colors)

## Centrality analysis but for specific receptor/ligand pairs
Centrality(OM_MHO_All.con2,
           modes.include = c('VEGF'),
           weight.attribute = 'weight_sc',
           min.z = 0,
           group.by = 'mechanism')

## CellCell Scatter to identify specific cell-to-cell interactions
p1 <- CellCellScatter(OM_MHO_All.con2,sources.include = 'endothelial',targets.include = 'activated_stellate',
                      label.threshold = 3,
                      weight.attribute = 'weight_sc',min.pct = 0.25,min.z = 2)
p1 <- p1 + xlim(0,NA) + ylim(0,NA)
p1

## SignalScatter clusters like CellCellScatter but links cells based on specific receptor-ligand pairs of interest
p2 <- SignalScatter(OM_MHO_All.con2, features = c('JAG1','JAG2','DLL4'),label.threshold = 1,weight.attribute = 'weight_sc',min.z = 2)
p2 <- p2 + xlim(2,NA) + ylim(2,NA)
p2

## CircosPlot step 1 = select top x (5) vectors for each cell-cell vector
test_OM_MHO <- OM_MHO_All.con2
test_OM_MHO <- data.frame(test_OM_MHO %>% group_by(vector) %>% top_n(5,weight_sc))

## then highlight cells of interest
cells.of.interest <- c('MSLN+','PI16+','PDGFRA+','Myofibroblasts', 'ECs', 'Myeloid cells', 'B cells', 'Plasma cells', 'T cells', 'NK cells')

## Then build four charts for these relationships
# Using edgeweights from normalized slot:
CircosPlot(test_OM_MHO,weight.attribute = 'weight_norm',sources.include = cells.of.interest,targets.include = cells.of.interest,lab.cex = 0.6,title = 'Edgeweights from normalized slot')

# Using edgeweights from scaled slot:
CircosPlot(test_OM_MHO,weight.attribute = 'weight_sc',sources.include = cells.of.interest,targets.include = cells.of.interest,lab.cex = 0.6,title = 'Edgeweights from scaled slot')

# Using separate ligand and receptor expression (from normalized slot)
CircosPlot(test_OM_MHO,weight.attribute = 'weight_norm',sources.include = cells.of.interest,targets.include = cells.of.interest,balanced.edges = F,lab.cex = 0.6,title = 'Ligand vs. receptor expression (from normalized slot)')

# Using separate ligand and receptor expression (from scaled slot)
CircosPlot(test_OM_MHO,weight.attribute = 'weight_sc',sources.include = cells.of.interest,targets.include = cells.of.interest,balanced.edges = F,lab.cex = 0.6,title = 'Ligand vs. receptor expression (from scaled slot)')

## Niche-wise investigations are to show all communication pathways converging on a given cell (i.e., focus on a specific cluster and how do all other cells in the tissue communicate specifically with that)
CircosPlot(test_OM_MHO,targets.include = 'MSLN+',lab.cex = 0.6)

## Or, use the same analysis but in reverse with Source-Wise, where it shows how your specific cell type of interest interacts with all other cell types
CircosPlot(test_OM_MHO,sources.include = 'MSLN+',lab.cex = 0.6)

## Large scale visualisation shows all connections between all receptor ligand pairs for all cells in the tissue
CircosPlot(OM_MHO_All.con2,min.z=1,lab.cex = 0.4,gap.degree = 0.1)



####START MUO ANALYSIS HERE

## Scale and make connectome
OM_MUO_All <- NormalizeData(OM_MUO_merge_copy)
connectome.genes <- union(Connectome::ncomms8866_human$Ligand.ApprovedSymbol,Connectome::ncomms8866_human$Receptor.ApprovedSymbol)
genes <- connectome.genes[connectome.genes %in% rownames(OM_MUO_All)]
OM_MUO_All <- ScaleData(OM_MUO_All,features = genes)
OM_MUO_All.con <- CreateConnectome(OM_MUO_All,species = 'human',min.cells.per.ident = 75,p.values = F,calculate.DOR = F)

## Filter edges
p1 <- ggplot(OM_MUO_All.con, aes(x=ligand.scale)) + geom_density() + ggtitle('Ligand.scale')
p2 <- ggplot(OM_MUO_All.con, aes(x=recept.scale)) + geom_density() + ggtitle('Recept.scale')
p3 <- ggplot(OM_MUO_All.con, aes(x=percent.target)) + geom_density() + ggtitle('Percent.target')
p4 <- ggplot(OM_MUO_All.con, aes(x=percent.source)) + geom_density() + ggtitle('Percent.source')
plot_grid(p1,p2,p3,p4)

OM_MUO_All.con2 <- FilterConnectome(OM_MUO_All.con,min.pct = 0.1,min.z = 0.25,remove.na = T)

## If you want to view a specific Receptor/Ligand, use this command
p1 <- NetworkPlot(OM_MUO_All.con2,features = 'VEGFA',min.pct = 0.1,weight.attribute = 'weight_sc',include.all.nodes = T)
p2 <- NetworkPlot(OM_MUO_All.con2,features = 'VEGFA',min.pct = 0.75,weight.attribute = 'weight_sc',include.all.nodes = T)
plot_grid(p1,p2,nrow=1)

## Centrality Analysis
Centrality(OM_MUO_All.con2,
           modes.include = NULL,
           min.z = NULL,
           weight.attribute = 'weight_sc',
           group.by = 'mode')

## Centrality analysis but for specific receptor/ligand pairs
Centrality(OM_MUO_All.con2,
           modes.include = c('VEGF'),
           weight.attribute = 'weight_sc',
           min.z = 0,
           group.by = 'mechanism')

## CellCell Scatter to identify specific cell-to-cell interactions
p1 <- CellCellScatter(OM_MUO_All.con2,sources.include = 'MSLN+',targets.include = 'T cells',
                      label.threshold = 3,
                      weight.attribute = 'weight_sc',min.pct = 0.25,min.z = 2)
p1 <- p1 + xlim(0,NA) + ylim(0,NA)
p1

## SignalScatter clusters like CellCellScatter but links cells based on specific receptor-ligand pairs of interest
p2 <- SignalScatter(OM_MUO_All.con2, features = c('JAG1','JAG2','DLL4'),label.threshold = 1,weight.attribute = 'weight_sc',min.z = 2)
p2 <- p2 + xlim(2,NA) + ylim(2,NA)
p2

## CircosPlot step 1 = select top x (5) vectors for each cell-cell vector
test_MUO_OM <- OM_MUO_All.con2
test_MUO_OM <- data.frame(test_MUO_OM %>% group_by(vector) %>% top_n(5,weight_sc))

## then highlight cells of interest
cells.of.interest <- c('MSLN+','PI16+','PDGFRA+','Myofibroblasts', 'ECs', 'Myeloid cells', 'B cells', 'Plasma cells', 'T cells', 'NK cells')

## Then build four charts for these relationships
# Using edgeweights from normalized slot:
CircosPlot(test_MUO_OM,weight.attribute = 'weight_norm',sources.include = cells.of.interest,targets.include = cells.of.interest,lab.cex = 0.6,title = 'Edgeweights from normalized slot')
# Using edgeweights from scaled slot:
CircosPlot(test_MUO_OM,weight.attribute = 'weight_sc',sources.include = cells.of.interest,targets.include = cells.of.interest,lab.cex = 0.6,title = 'Edgeweights from scaled slot')
# Using separate ligand and receptor expression (from normalized slot)
CircosPlot(test_MUO_OM,weight.attribute = 'weight_norm',sources.include = cells.of.interest,targets.include = cells.of.interest,balanced.edges = F,lab.cex = 0.6,title = 'Ligand vs. receptor expression (from normalized slot)')
# Using separate ligand and receptor expression (from scaled slot)
CircosPlot(test_MUO_OM,weight.attribute = 'weight_sc',sources.include = cells.of.interest,targets.include = cells.of.interest,balanced.edges = F,lab.cex = 0.6,title = 'Ligand vs. receptor expression (from scaled slot)')

## Niche-wise investigations are to show all communication pathways converging on a given cell (i.e., focus on a specific cluster and how do all other cells in the tissue communicate specifically with that)
CircosPlot(test_MUO_OM,targets.include = 'MSLN+',lab.cex = 0.6)

## Or, use the same analysis but in reverse with Source-Wise, where it shows how your specific cell type of interest interacts with all other cell types
CircosPlot(test_MUO_OM,sources.include = 'MSLN+',lab.cex = 0.6)

## Large scale visualisation shows all connections between all receptor ligand pairs for all cells in the tissue
CircosPlot(OM_MUO_All.con2,min.z=1,lab.cex = 0.4,gap.degree = 0.1)




################### SUBCUT CONNECTOME ANALYSIS

SQ_merge_copy <- SQ_merge
Idents(SQ_merge_copy) <- "Sample"
SQ_MHO_merge_copy <- subset(SQ_merge_copy, idents = c("SQ_MHO"))
Idents(SQ_MHO_merge_copy) <- "seurat_clusters"
SQ_MUO_merge_copy <- subset(SQ_merge_copy, idents = c("SQ_MUO"))
Idents(SQ_MUO_merge_copy) <- "seurat_clusters"
DimPlot(SQ_MHO_merge_copy)
DimPlot(SQ_MUO_merge_copy)

## Rename Clusters
SQ_MHO_MUO_IDs <- c("PDGFRA+", "PDGFRA+", "PDGFRA+", "T Cell", "Myeloid", "NK", "EC", "APC", "SMC", "T Cell", "Myeloid", "?")
names(SQ_MHO_MUO_IDs) <- levels(SQ_MHO_merge_copy)
SQ_MHO_merge_copy <- RenameIdents(SQ_MHO_merge_copy, SQ_MHO_MUO_IDs)
DimPlot(SQ_MHO_merge_copy, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(SQ_MHO_merge_copy, reduction = "umap", pt.size = 0.5) + NoLegend()


## MHO CONNECTOME ANALYSIS PERFORMED FIRST - MUO REQUIRES COPY AND CHANGE MHO TO MUO THROUGHOUT ###

## Scale and make connectome
SQ_MHO_All <- NormalizeData(SQ_MHO_merge_copy)
connectome.genes <- union(Connectome::ncomms8866_human$Ligand.ApprovedSymbol,Connectome::ncomms8866_human$Receptor.ApprovedSymbol)
genes <- connectome.genes[connectome.genes %in% rownames(SQ_MHO_All)]
SQ_MHO_All <- ScaleData(SQ_MHO_All,features = genes)
SQ_MHO_All.con <- CreateConnectome(SQ_MHO_All,species = 'human',min.cells.per.ident = 75,p.values = F,calculate.DOR = F)

## Filter edges
p1 <- ggplot(SQ_MHO_All.con, aes(x=ligand.scale)) + geom_density() + ggtitle('Ligand.scale')
p2 <- ggplot(SQ_MHO_All.con, aes(x=recept.scale)) + geom_density() + ggtitle('Recept.scale')
p3 <- ggplot(SQ_MHO_All.con, aes(x=percent.target)) + geom_density() + ggtitle('Percent.target')
p4 <- ggplot(SQ_MHO_All.con, aes(x=percent.source)) + geom_density() + ggtitle('Percent.source')
plot_grid(p1,p2,p3,p4)

SQ_MHO_All.con2 <- FilterConnectome(SQ_MHO_All.con,min.pct = 0.1,min.z = 0.25,remove.na = T)

## If you want to view a specific Receptor/Ligand, use this command
p1 <- NetworkPlot(SQ_MHO_All.con2,features = 'C3',min.pct = 0.1,weight.attribute = 'weight_sc',include.all.nodes = T)
p2 <- NetworkPlot(SQ_MHO_All.con2,features = 'C3',min.pct = 0.75,weight.attribute = 'weight_sc',include.all.nodes = T)
plot_grid(p1,p2,nrow=1)

## Centrality Analysis
Centrality(SQ_MHO_All.con2,
           modes.include = NULL,
           min.z = NULL,
           weight.attribute = 'weight_sc',
           group.by = 'mode')

## Centrality analysis but for specific receptor/ligand pairs
Centrality(SQ_MHO_All.con2,
           modes.include = c('VEGF'),
           weight.attribute = 'weight_sc',
           min.z = 0,
           group.by = 'mechanism')

## CellCell Scatter to identify specific cell-to-cell interactions
p1 <- CellCellScatter(SQ_MHO_All.con2,sources.include = 'Myeloid',targets.include = 'T Cell',
                      label.threshold = 3,
                      weight.attribute = 'weight_sc',min.pct = 0.25,min.z = 2)
p1 <- p1 + xlim(0,NA) + ylim(0,NA)
p1

## SignalScatter clusters like CellCellScatter but links cells based on specific receptor-ligand pairs of interest
p2 <- SignalScatter(OM_MHO_All.con2, features = c('JAG1','JAG2','DLL4'),label.threshold = 1,weight.attribute = 'weight_sc',min.z = 2)
p2 <- p2 + xlim(2,NA) + ylim(2,NA)
p2

## CircosPlot step 1 = select top x (5) vectors for each cell-cell vector
test_MHO_All <- SQ_MHO_All.con2
test_MHO_All <- data.frame(test_MHO_All %>% group_by(vector) %>% top_n(5,weight_sc))

## then highlight cells of interest
cells.of.interest_SQ <- c('T Cell','APC','PDGFRA+','SMC', 'EC', 'Myeloid', '?', 'NK')

## Then build four charts for these relationships
# Using edgeweights from normalized slot:
CircosPlot(test_MHO_All,weight.attribute = 'weight_norm',sources.include = cells.of.interest_SQ,targets.include = cells.of.interest_SQ,lab.cex = 0.6,title = 'Edgeweights from normalized slot')

# Using edgeweights from scaled slot:
CircosPlot(test_MHO_All,weight.attribute = 'weight_sc',sources.include = cells.of.interest_SQ,targets.include = cells.of.interest_SQ,lab.cex = 0.6,title = 'Edgeweights from scaled slot')

# Using separate ligand and receptor expression (from normalized slot)
CircosPlot(test_MHO_All,weight.attribute = 'weight_norm',sources.include = cells.of.interest_SQ,targets.include = cells.of.interest_SQ,balanced.edges = F,lab.cex = 0.6,title = 'Ligand vs. receptor expression (from normalized slot)')

# Using separate ligand and receptor expression (from scaled slot)
CircosPlot(test_MHO_All,weight.attribute = 'weight_sc',sources.include = cells.of.interest_SQ,targets.include = cells.of.interest_SQ,balanced.edges = F,lab.cex = 0.6,title = 'Ligand vs. receptor expression (from scaled slot)')

## Niche-wise investigations are to show all communication pathways converging on a given cell (i.e., focus on a specific cluster and how do all other cells in the tissue communicate specifically with that)
CircosPlot(test_MHO_All,targets.include = 'endothelial',lab.cex = 0.6)

## Or, use the same analysis but in reverse with Source-Wise, where it shows how your specific cell type of interest interacts with all other cell types
CircosPlot(test_MHO_All,sources.include = 'endothelial',lab.cex = 0.6)

## Large scale visualisation shows all connections between all receptor ligand pairs for all cells in the tissue
CircosPlot(SQ_MHO_All.con2,min.z=1,lab.cex = 0.4,gap.degree = 0.1)



####START MUO ANALYSIS HERE

names(SQ_MHO_MUO_IDs) <- levels(SQ_MUO_merge_copy)
SQ_MUO_merge_copy <- RenameIdents(SQ_MUO_merge_copy, SQ_MHO_MUO_IDs)
DimPlot(SQ_MUO_merge_copy, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(SQ_MUO_merge_copy, reduction = "umap", pt.size = 0.5) + NoLegend()


## Scale and make connectome
SQ_MUO_All <- NormalizeData(SQ_MUO_merge_copy)
connectome.genes <- union(Connectome::ncomms8866_human$Ligand.ApprovedSymbol,Connectome::ncomms8866_human$Receptor.ApprovedSymbol)
genes <- connectome.genes[connectome.genes %in% rownames(SQ_MUO_All)]
SQ_MUO_All <- ScaleData(SQ_MUO_All,features = genes)
SQ_MUO_All.con <- CreateConnectome(SQ_MUO_All,species = 'human',min.cells.per.ident = 75,p.values = F,calculate.DOR = F)

## Filter edges
p1 <- ggplot(SQ_MUO_All.con, aes(x=ligand.scale)) + geom_density() + ggtitle('Ligand.scale')
p2 <- ggplot(SQ_MUO_All.con, aes(x=recept.scale)) + geom_density() + ggtitle('Recept.scale')
p3 <- ggplot(SQ_MUO_All.con, aes(x=percent.target)) + geom_density() + ggtitle('Percent.target')
p4 <- ggplot(SQ_MUO_All.con, aes(x=percent.source)) + geom_density() + ggtitle('Percent.source')
plot_grid(p1,p2,p3,p4)

SQ_MUO_All.con2 <- FilterConnectome(SQ_MUO_All.con,min.pct = 0.1,min.z = 0.25,remove.na = T)

## If you want to view a specific Receptor/Ligand, use this command
p1 <- NetworkPlot(SQ_MUO_All.con2,features = 'VEGFA',min.pct = 0.1,weight.attribute = 'weight_sc',include.all.nodes = T)
p2 <- NetworkPlot(SQ_MUO_All.con2,features = 'VEGFA',min.pct = 0.75,weight.attribute = 'weight_sc',include.all.nodes = T)
plot_grid(p1,p2,nrow=1)

## Centrality Analysis
Centrality(SQ_MUO_All.con2,
           modes.include = NULL,
           min.z = NULL,
           weight.attribute = 'weight_sc',
           group.by = 'mode')

## Centrality analysis but for specific receptor/ligand pairs
Centrality(OM_MUO_All.con2,
           modes.include = c('VEGF'),
           weight.attribute = 'weight_sc',
           min.z = 0,
           group.by = 'mechanism')

## CellCell Scatter to identify specific cell-to-cell interactions
p1 <- CellCellScatter(SQ_MUO_All.con2,sources.include = 'endothelial',targets.include = 'activated_stellate',
                      label.threshold = 3,
                      weight.attribute = 'weight_sc',min.pct = 0.25,min.z = 2)
p1 <- p1 + xlim(0,NA) + ylim(0,NA)
p1

## SignalScatter clusters like CellCellScatter but links cells based on specific receptor-ligand pairs of interest
p2 <- SignalScatter(SQ_MUO_All.con2, features = c('JAG1','JAG2','DLL4'),label.threshold = 1,weight.attribute = 'weight_sc',min.z = 2)
p2 <- p2 + xlim(2,NA) + ylim(2,NA)
p2

## CircosPlot step 1 = select top x (5) vectors for each cell-cell vector
test_MUO_SQ <- SQ_MUO_All.con2
test_MUO_SQ <- data.frame(test_MUO_SQ %>% group_by(vector) %>% top_n(5,weight_sc))

## then highlight cells of interest
cells.of.interest_MUO_SQ <- c('T Cell','APC','PDGFRA+','SMC', 'EC', 'Myeloid', '?', 'NK')

## Then build four charts for these relationships
# Using edgeweights from normalized slot:
CircosPlot(test_MUO_SQ,weight.attribute = 'weight_norm',sources.include = cells.of.interest_MUO_SQ,targets.include = cells.of.interest_MUO_SQ,lab.cex = 0.6,title = 'Edgeweights from normalized slot')

# Using edgeweights from scaled slot:
CircosPlot(test_MUO_SQ,weight.attribute = 'weight_sc',sources.include = cells.of.interest_MUO_SQ,targets.include = cells.of.interest_MUO_SQ,lab.cex = 0.6,title = 'Edgeweights from scaled slot')

# Using separate ligand and receptor expression (from normalized slot)
CircosPlot(test_MUO_SQ,weight.attribute = 'weight_norm',sources.include = cells.of.interest_MUO_SQ,targets.include = cells.of.interest_MUO_SQ,balanced.edges = F,lab.cex = 0.6,title = 'Ligand vs. receptor expression (from normalized slot)')

# Using separate ligand and receptor expression (from scaled slot)
CircosPlot(test_MUO_SQ,weight.attribute = 'weight_sc',sources.include = cells.of.interest_MUO_SQ,targets.include = cells.of.interest_MUO_SQ,balanced.edges = F,lab.cex = 0.6,title = 'Ligand vs. receptor expression (from scaled slot)')

## Niche-wise investigations are to show all communication pathways converging on a given cell (i.e., focus on a specific cluster and how do all other cells in the tissue communicate specifically with that)
CircosPlot(test_MUO_SQ,targets.include = 'MSLN+',lab.cex = 0.6)

## Or, use the same analysis but in reverse with Source-Wise, where it shows how your specific cell type of interest interacts with all other cell types
CircosPlot(test_MUO_SQ,sources.include = 'MSLN+',lab.cex = 0.6)

## Large scale visualisation shows all connections between all receptor ligand pairs for all cells in the tissue
CircosPlot(test_MUO_SQ,min.z=1,lab.cex = 0.4)
