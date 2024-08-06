##Reference codes for manuscript colour schemes - Load in wesanderson palette first 
library(wesanderson)

# Discrete color
bp + scale_fill_manual(values = wes_palette("GrandBudapest1", n = 3))

# Gradient color
pal <- wes_palette("Zissou1", 100, type = "continuous")
ggplot(heatmap, aes(x = X2, y = X1, fill = value)) +
  geom_tile() + 
  scale_fill_gradientn(colours = pal) + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  coord_equal() 

###Colours Used:
OM_colors <- c("PDGFRA+" = "coral",
               "PI16+" = "#94d2bd",
               "MSLN+" = "violet",
               "Myofibroblasts" = "dodgerblue",
               "ECs" = "#e9d8a6",
               "T cells" = "#ee9b00",
               "B cells" = "#ca6702", 
               "Plasma cells" = "purple2",
               "NK cells" = "#005f73",
               "Myeloid cells" = "#ae2012")

SQ_colors <- c("PDGFRA+" = "coral",
               "PI16+" = "#94d2bd",
               "Myofibroblasts" = "dodgerblue",
               "ECs" = "#e9d8a6",
               "T cells" = "#ee9b00",
               "B cells" = "#ca6702", 
               "NK cells" = "#005f73",
               "Myeloid cells" = "#ae2012")

OM_Myeloid_colors <- c('FBLN1+' = '#005f73', 'LAM' = 'deeppink', 'cDC2As' = 'peru', 'cDC2Bs' = 'chartreuse3', 'PVMs' = 'skyblue2', 'THBS1+' = 'gold', 'Monocytes' = 'purple')
SQ_Myeloid_colors <- c('LAM' = 'deeppink', 'cDC2As' = 'peru', 'cDC2Bs' = 'chartreuse3', 'PVM' = 'skyblue2', 'Monocytes' = 'purple', 'cDC1s' = 'dodgerblue')

SQ_Lympho_colors <- c('CD4+ T' = '#005f73', 'CD8+ T' = 'deeppink', 'NKs' = 'chartreuse3', 'B cells' = 'skyblue2', 'GD T' = 'purple', 'Tregs' = 'coral')
OM_Lympho_colors <- c('CD4+ T' = '#005f73', 'CD8+ T' = 'deeppink', 'mNKs' = 'peru', 'NK-like' = 'chartreuse3', 'B cells' = 'skyblue2', 'NA CD4+ T' = 'gold', 'GD T' = 'purple', 'ILC-like' = 'dodgerblue', 'Plasma cells' = '#e9d8a6')

##General MHO MUO colours

## MHO = mediumseagreen
## MUO = royalblue1