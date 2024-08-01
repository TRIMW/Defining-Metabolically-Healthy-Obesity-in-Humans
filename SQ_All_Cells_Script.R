##Splitting OM data file into groups and re-clustering
##Splitting a file by group
SQ_Copy <- SQ
Idents(SQ_Copy) <- "Sample"
DimPlot(SQ_Copy)

##Isolate Healthy Obese Group
SQ_MHO_Copy <- subset(SQ_Copy, idents = c("SQ_MHO"))
Idents(SQ_MHO_Copy) <- "seurat_clusters"
DimPlot(SQ_MHO_Copy)

##Isolate unhealthy obese group
SQ_MUO_Copy <- subset(SQ_Copy, idents = c("SQ_MUO"))
Idents(SQ_MUO_Copy) <- "seurat_clusters"
DimPlot(SQ_MUO_Copy)

##Relabel Clusters MHO then MUO (Same code)
SQ_new.cluster.ids <- c("PDGFRA+", "PDGFRA+", "PDGFRA+", "T cells", "Myeloid cells", "NK cells", "ECs", "PI16+", "SMCs", "T cells", "Myeloid cells", "B cells")

names(SQ_new.cluster.ids) <- levels(SQ_MHO_Copy)

SQ_MHO_Copy <- RenameIdents(SQ_MHO_Copy, SQ_new.cluster.ids)

names(SQ_new.cluster.ids) <- levels(SQ_MUO_Copy)

SQ_MUO_Copy <- RenameIdents(SQ_MUO_Copy, SQ_new.cluster.ids)

names(SQ_new.cluster.ids) <- levels(SQ)

SQ <- RenameIdents(SQ, SQ_new.cluster.ids)

##Colour Scheme
SQ_colors <- c("PDGFRA+" = "coral",
               "PI16+" = "#94d2bd",
               "SMCs" = "dodgerblue",
               "ECs" = "#e9d8a6",
               "T cells" = "#ee9b00",
               "B cells" = "#ca6702", 
               "NK cells" = "#005f73",
               "Myeloid cells" = "#ae2012")

##Re-order clusters
SQ_levels <- c('PDGFRA+', 'PI16+', 'SMCs', 'ECs', 'T cells', 'NK cells', 'B cells', 'Myeloid cells')

Idents(SQ) <- factor(Idents(SQ), levels= SQ_levels)

Idents(SQ_MHO_Copy) <- factor(Idents(SQ_MHO_Copy), levels= SQ_levels)

SCpubr::do_DimPlot(sample = SQ_MHO_Copy,
                   colors.use = SQ_colors, plot_cell_borders = TRUE, border.size = 2, border.color = 'black', pt.size = 0.7, label = TRUE, repel = TRUE)

Idents(SQ_MUO_Copy) <- factor(Idents(SQ_MUO_Copy), levels= SQ_levels)

SCpubr::do_DimPlot(sample = SQ_MUO_Copy,
                   colors.use = SQ_colors, plot_cell_borders = TRUE, border.size = 2, border.color = 'black', pt.size = 0.6, label = TRUE, repel = TRUE)


vln_plot <- VlnPlot(SQ, 
                    features = c("CD14", "JCHAIN", "FCGR3A", "CD3E", "PTPRC", "VWF", "RGS5", "PI16", "PDGFRA"),
                    stack = TRUE, 
                    flip = TRUE, 
                    cols = c("PDGFRA" = "coral", "PI16" = "#94D2BD", "RGS5" = "dodgerblue", 
                             "VWF" = "#e9d8a6", "PTPRC" = "yellow", "CD3E" = "#ee9b00", 
                             "JCHAIN" = "#ca6702", "FCGR3A" = "#005f73", "CD14" = "#ae2012"))



SQ_MHO_Cells <- SQ_MHO_Copy@active.ident %>% as.data.table
View(SQ_MHO_Cells)
write.xlsx2(SQ_MHO_Cells, 'SQ_MHO_Counts.xlsx', sheetName = "Sheet1",
            col.names = TRUE, row.names = TRUE, append = FALSE)

SQ_MUO_Cells <- SQ_MUO_Copy@active.ident %>% as.data.table
View(SQ_MUO_Cells)
write.xlsx2(SQ_MUO_Cells, 'SQ_MUO_Counts.xlsx', sheetName = "Sheet1",
            col.names = TRUE, row.names = TRUE, append = FALSE)
