##HALLMARK Analysis between clusters for MSLN+ cells

library(msigdbr)
library(SCPA)
library(tidyverse)

pathways <- msigdbr("Homo sapiens", "H") %>% format_pathways()

OM_MSC_MHO <- seurat_extract(OM_MSCs,
                             meta1 = "Sample", value_meta1 = "OM_MHO")

OM_MSC_MUO <- seurat_extract(OM_MSCs,
                             meta1 = "Sample", value_meta1 = "OM_MUO")



OM_MSC_Comparison <- compare_pathways(samples = list(OM_MSC_MHO, OM_MSC_MUO), pathways = pathways)