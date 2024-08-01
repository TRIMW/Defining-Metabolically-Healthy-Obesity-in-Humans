##HALLMARK Analysis between clusters

library(msigdbr)
library(SCPA)
library(tidyverse)

pathways <- msigdbr("Homo sapiens", "H") %>% format_pathways()

OM_Lympho_MHO <- seurat_extract(OM_Lymph,
                             meta1 = "Sample", value_meta1 = "OM_MHO")
OM_Lympho_MUO <- seurat_extract(OM_Lymph,
                             meta1 = "Sample", value_meta1 = "OM_MUO")

OM_Lympho_Comparison <- compare_pathways(samples = list(OM_Lympho_MHO, OM_Lympho_MUO), pathways = pathways)

