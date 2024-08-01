##Omentum processing:

##Identifying and labelling clusters - Whole Sample
OM_levels <- c(OM_levels <- c('PreAd1', 'PreAd2', 'PreAd3', 'PreAd4', 'APCs', 'MSC1s', 'MSC2s', 'MSC3s', 'SMCs', 'ECs', 'T cells', 'B cells', '?', 'NK cells', 'Myeloid cells'))
Idents(OM) <- factor(Idents(OM), levels= OM_levels)

OM_new.cluster.ids <- c('PDGFRA+', 'PDGFRA+', 'PDGFRA+', 'PDGFRA+', 'PI16+', 'MSLN+', 'MSLN+', 'MSLN+', 'SMCs', 'ECs', 'T cells', 'B cells', 'Plasma cells', 'NK cells', 'Myeloid cells')
names(OM_new.cluster.ids) <- levels(OM)
OM <- RenameIdents(OM, OM_new.cluster.ids)

##Colour scheme
OM_colors <- c("PDGFRA+" = "coral",
               "PI16+" = "#94d2bd",
               "MSLN+" = "violet",
               "SMCs" = "dodgerblue",
               "ECs" = "#e9d8a6",
               "T cells" = "#ee9b00",
               "B cells" = "#ca6702", 
               "Plasma cells" = "purple2",
               "NK cells" = "#005f73",
               "Myeloid cells" = "#ae2012")

##Plot UMAP:
SCpubr::do_DimPlot(sample = OM,
                   colors.use = OM_colors, plot_cell_borders = TRUE, border.size = 2, border.color = 'black', pt.size = 0.6, label = TRUE, repel = TRUE)

##Make violin plot of cell marker genes
library(cowplot)

vln_plot <- VlnPlot(OM, 
                    features = c("CD14", "FCGR3A", "JCHAIN", "CD79A", "CD3E", "PTPRC", "VWF", "RGS5", "MSLN", "PI16", "PDGFRA"),
                    stack = TRUE, 
                    flip = TRUE, 
                    cols = c("PDGFRA" = "coral", "PI16" = "#94D2BD", "MSLN" = "violet", "RGS5" = "dodgerblue", 
                             "VWF" = "#e9d8a6", "PTPRC" = "yellow", "CD3E" = "#ee9b00", 
                             "CD79A" = "#ca6702", "JCHAIN" = "purple2", "FCGR3A" = "#005f73", "CD14" = "#ae2012"))

OM_colors <- c("PDGFRA+" = "coral",
               "PI16+" = "#94d2bd",
               "MSLN+" = "violet",
               "SMCs" = "dodgerblue",
               "ECs" = "#e9d8a6",
               "T cells" = "#ee9b00",
               "B cells" = "#ca6702", 
               "Plasma cells" = "purple2",
               "NK cells" = "#005f73",
               "Myeloid cells" = "#ae2012")


# Remove the legend
vln_plot <- vln_plot + theme(legend.position = "none")
print(vln_plot)


##Splitting OM data file into groups and re-clustering
##Splitting a file by group
OM_Copy <- OM
Idents(OM_Copy) <- "Sample"
DimPlot(OM_Copy)

##Isolate Healthy Obese Group
OM_MHO_Copy <- subset(OM_Copy, idents = c("OM_MHO"))
Idents(OM_MHO_Copy) <- "seurat_clusters"
DimPlot(OM_MHO_Copy)

##Isolate unhealthy obese group
OM_MUO_Copy <- subset(OM_Copy, idents = c("OM_MUO"))
Idents(OM_MUO_Copy) <- "seurat_clusters"
DimPlot(OM_MUO_Copy)

##Relabel Clusters MHO then MUO (Same code)
OM_new.cluster.ids <- c("T cells", "PDGFRA+", "T cells", "MSLN+", "PI16+", "PDGFRA+", "MSLN+", "Myeloid cells", "Myeloid cells", "NK cells", "PDGFRA+", "ECs", "SMCs", "B cells", "PDGFRA+", "MSLN+", "Plasma cells")

names(OM_new.cluster.ids) <- levels(OM_MHO_Copy)

OM_MHO_Copy <- RenameIdents(OM_MHO_Copy, OM_new.cluster.ids)

names(OM_new.cluster.ids) <- levels(OM_MUO_Copy)

OM_MUO_Copy <- RenameIdents(OM_MUO_Copy, OM_new.cluster.ids)

##Re-order clusters
OM_levels <- c('PDGFRA+', 'PI16+', 'MSLN+', 'ECs', 'SMCs', 'NK cells', 'T cells', 'B cells', 'Plasma cells', 'Myeloid cells')

Idents(OM_MHO_Copy) <- factor(Idents(OM_MHO_Copy), levels= OM_levels)

SCpubr::do_DimPlot(sample = OM_MHO_Copy,
                   colors.use = OM_colors, plot_cell_borders = TRUE, border.size = 2, border.color = 'black', pt.size = 0.6, label = TRUE, repel = TRUE)

Idents(OM_MUO_Copy) <- factor(Idents(OM_MUO_Copy), levels= OM_levels)

SCpubr::do_DimPlot(sample = OM_MUO_Copy,
                   colors.use = OM_colors, plot_cell_borders = TRUE, border.size = 2, border.color = 'black', pt.size = 0.6, label = TRUE, repel = TRUE)

##Extract cell numbers per cluster
library(dplyr)
library(magrittr)
library(xlxs)


OM_MHO_Cells <- OM_MHO_Copy@active.ident %>% as.data.table
View(OM_MHO_Cells)
write.xlsx2(OM_MHO_Cells, 'OM_MHO_Counts.xlsx', sheetName = "Sheet1",
            col.names = TRUE, row.names = TRUE, append = FALSE)

OM_MUO_Cells <- OM_MUO_Copy@active.ident %>% as.data.table
View(OM_MUO_Cells)
write.xlsx2(OM_MUO_Cells, 'OM_MUO_Counts.xlsx', sheetName = "Sheet1",
            col.names = TRUE, row.names = TRUE, append = FALSE)

