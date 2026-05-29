library(data.table)
library(fgsea)
library(grid)
library(dplyr)

RES_DIR <- "/Users/willtrim/Documents/projs/Wills_analysis/obesity_scRNAseq/results/DE/subq/"

cell_type_col <- "cell_type2_"
gene_id_col <- "entrez_id"
n_top <- 15

out_dir <- paste0(RES_DIR, "fgsea/", cell_type_col, "/MAST")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

setwd(out_dir)

de_res <- read.csv(
  paste0(
    RES_DIR,
    "MAST_results_per_",
    cell_type_col,
    ".csv"),
)


de_res_by_ct <- split(de_res, de_res$cell_type)


plot_fgsea_top <- function(fgseaRes, pathways, ranks, n_top) {
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=n_top), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=n_top), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  plotGseaTable(pathways[topPathways], ranks, fgseaRes, gseaParam = 0.5)
}

###########################################################################################
###########################################################################################

msigdb_gmt_fp <- "/Users/willtrim/Documents/reference/msigdb/msigdb_v2025.1.Hs_GMTs/"


### Reactome


gmt.file <- paste(msigdb_gmt_fp, "c2.cp.reactome.v2025.1.Hs.entrez.gmt", sep="")
pathways <- gmtPathways(gmt.file)
str(head(pathways))

###  Hallmark 

hallmark.gmt.file <- paste(msigdb_gmt_fp, "h.all.v2025.1.Hs.entrez.gmt", sep="")
hallmark.pathways <- gmtPathways(hallmark.gmt.file)
str(head(hallmark.pathways))

### Wikipathways
wiki.gmt.file <- paste(msigdb_gmt_fp, "c2.cp.wikipathways.v2025.1.Hs.entrez.gmt", sep="")
wiki.pathways <- gmtPathways(wiki.gmt.file)
str(head(wiki.pathways))


### KEGG

library(clusterProfiler)

kegg_pathways_df <- download_KEGG(species = "hsa", keggType = "KEGG", keyType = "kegg")

kegg_pathways_list <- kegg_pathways_df["KEGGPATHID2EXTID"][[1]] %>%
  group_by(from) %>%
  group_split()

library(purrr)
kegg_pathways <- map(kegg_pathways_list, ~pull(.x, to))
kegg_pathway_ids <-  map(kegg_pathways_list, ~ .x[[1, 1]])
kegg_pathways_names_map <- kegg_pathways_df[["KEGGPATHID2NAME"]]$to
names(kegg_pathways_names_map) <- kegg_pathways_df[["KEGGPATHID2NAME"]]$from
kegg_pathways_names <- lapply(kegg_pathway_ids, function(x) {kegg_pathways_names_map[[x]]})
names(kegg_pathways) <- kegg_pathways_names

### GTRD

gtrd.gmt.file <- paste(msigdb_gmt_fp, "c3.tft.gtrd.v2025.1.Hs.entrez.gmt", sep="")
gtrd.pathways <- gmtPathways(gtrd.gmt.file)
str(head(gtrd.pathways))


### GO

#### BP
bp.gmt.file <- paste(msigdb_gmt_fp, "c5.go.bp.v2025.1.Hs.entrez.gmt", sep="")
bp.pathways <- gmtPathways(bp.gmt.file)
str(head(bp.pathways))


#### CC
cc.gmt.file <- paste(msigdb_gmt_fp, "c5.go.cc.v2025.1.Hs.entrez.gmt", sep="")
cc.pathways <- gmtPathways(cc.gmt.file)
str(head(cc.pathways))


#### MF

mf.gmt.file <- paste(msigdb_gmt_fp, "c5.go.mf.v2025.1.Hs.entrez.gmt", sep="")
mf.pathways <- gmtPathways(mf.gmt.file)
str(head(mf.pathways))

###########################################################################################
###########################################################################################

already_done <- c()

plot_top <- FALSE

# for(cell_type in names(de_res_by_ct)) {
for(cell_type in c("VECs")) {
  if(cell_type %in% already_done) {
    next
  }
  cell_type_de <- de_res_by_ct[[cell_type]]
  cell_type_de <- cell_type_de[!is.na(cell_type_de[,gene_id_col]), ]
  cell_type_de <- cell_type_de[!duplicated(cell_type_de[,gene_id_col]), ]
  
  # rownames(cell_type_de) <- cell_type_de[,gene_id_col]
  
  ranks <- cell_type_de[,"score"]
  names(ranks) <- cell_type_de[,gene_id_col]

  fgseaRes_react <- fgsea(pathways = pathways, stats = ranks, minSize  = 15, maxSize  = 500, eps=0)
  fwrite(fgseaRes_react, paste0(cell_type, ".reactome.csv"), row.names = FALSE)
  
  fgseaRes_hallmark <- fgsea(pathways = hallmark.pathways, stats = ranks, minSize  = 15, maxSize  = 500, eps=0)
  fwrite(fgseaRes_hallmark, paste0(cell_type, ".hallmark.csv"), row.names = FALSE)
  
  fgseaRes_wp <- fgsea(pathways = wiki.pathways, stats = ranks, minSize  = 15, maxSize  = 500, eps=0)
  fwrite(fgseaRes_wp, paste0(cell_type, ".wiki.csv"), row.names = FALSE)
  
  fgseaRes_kegg <- fgsea(pathways = kegg_pathways, stats = ranks, minSize  = 15, maxSize  = 500, eps=0)
  fwrite(fgseaRes_kegg, paste0(cell_type, ".kegg.csv"), row.names = FALSE)
  
  fgseaRes_bp <- fgsea(pathways = bp.pathways, stats = ranks, minSize  = 15, maxSize  = 500, eps=0)
  fwrite(fgseaRes_bp, paste0(cell_type, ".bp.csv"), row.names = FALSE)
  
  fgseaRes_cc <- fgsea(pathways = cc.pathways, stats = ranks, minSize  = 15, maxSize  = 500, eps=0)
  fwrite(fgseaRes_cc, paste0(cell_type, ".cc.csv"), row.names = FALSE)
  
  fgseaRes_mf <- fgsea(pathways = mf.pathways, stats = ranks, minSize  = 15, maxSize  = 500, eps=0)
  fwrite(fgseaRes_mf, paste0(cell_type, ".mf.csv"), row.names = FALSE)

  if(plot_top) {
    p1 <- plot_fgsea_top(fgseaRes_hallmark, hallmark.pathways, ranks, n_top)
    p2 <- plot_fgsea_top(fgseaRes_react, pathways, ranks, n_top)
    p3 <- plot_fgsea_top(fgseaRes_wp, wiki.pathways, ranks, n_top)
    p4 <- plot_fgsea_top(fgseaRes_kegg, kegg_pathways, ranks, n_top)
    p5 <- plot_fgsea_top(fgseaRes_bp, bp.pathways, ranks, n_top)
  
  
    plots <- list(p1, p2, p3, p4, p5)
  
    pdf(paste0(cell_type, ".pdf"), width = 16, height = 16)
    for (i in seq_along(plots)) {
        # if (i > 1) grid.newpage()
        grid.draw(plots[[i]])
    }
    dev.off()
  }
}

