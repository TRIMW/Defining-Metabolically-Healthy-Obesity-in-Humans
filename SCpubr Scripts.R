https://enblacar.github.io/SCpubr-book/functions/DimPlots.html

Using SCpubr

We can also highlight whole identities with idents.highlight parameter. For this, just provide the desired identities to be selected. It can also work in combination with cells.highlight.

SCpubr::do_DimPlot(sample = sample, plot_cell_borders = TRUE,
                   border.size = 1.25)

library(SCpubr)
library(devtools)
library(ggplot2)
library(dplyr)

##Group by variable
SCpubr::do_DimPlot(sample = OM_merge,
                   group.by = "Sample", plot_cell_borders = TRUE,
                   border.size = 1.4)

##replace 'group.by = ' with 'split.by = ' to get dual plots
SCpubr::do_DimPlot(sample = OM_merge,
                   split.by = "Sample", plot_cell_borders = TRUE,
                   border.size = 1.4)

##Set colours
# Define a set of colors.
colors <- c("0" = "#001219",
            "1" = "#005f73",
            "2" = "#0a9396",
            "3" = "#94d2bd",
            "4" = "#e9d8a6",
            "5" = "#ee9b00",
            "6" = "#ca6702",
            "7" = "#bb3e03",
            "8" = "#ae2012",
            "9" = "#9b2226")

# Label the clusters - text geom.
SCpubr::do_DimPlot(sample = sample,
                   colors.use = colors)

##Show only specific subsets of clusters:
SCpubr::do_DimPlot(sample = OM_merge,
                   idents.keep = c("0", "2"))

# Using cells.highlight.
p1 <- SCpubr::do_DimPlot(sample = OM_merge, 
                         cells.highlight = cells.use)

# Using idents.highlight.
p2 <- SCpubr::do_DimPlot(sample = OM_merge, 
                         idents.highlight = c("6"))

# Using both.
p3 <- SCpubr::do_DimPlot(sample = OM_merge, 
                         cells.highlight = cells.use, 
                         idents.highlight = c("6"), plot_cell_borders = TRUE,
                         border.size = 1.25)

p <- p1 | p2 | p3
p

## For different colours across a gradient, add in:
diverging.palette = 

  ### Followed by:
'Spectral', 'RdYlGn', 'RdGr', 'RdBu', 'PuOr', 'PRGn', 'PiYG' , 'PiYG', 'BrBG' 

##For example"

SCpubr::do_FeaturePlot(sample,
                       features = "PC_1", 
                       enforce_symmetry = TRUE,
                       diverging.palette = "BrBG")


##Axes on bar plots can be flipped:
(flip = TRUE)

##Custom titles can be added:
SCpubr::do_BarPlot(sample = sample,
                   group.by = "seurat_clusters",
                   split.by = "annotation",
                   position = "fill",
                   flip = TRUE,
                   plot.title = "This is a title",
                   plot.subtitle = "This is a subtitle",
                   plot.caption = "This is a caption",
                   xlab = "My X axis title",
                   ylab = "My Y axis title",
                   legend.title = "My custom title")

##Plot fonts can be customised
SCpubr::do_BarPlot(sample = sample,
                   group.by = "seurat_clusters",
                   split.by = "annotation",
                   position = "fill",
                   flip = TRUE,
                   plot.title = "This is a title",
                   plot.subtitle = "This is a subtitle",
                   plot.caption = "This is a caption",
                   xlab = "My X axis title",
                   ylab = "My Y axis title",
                   legend.title = "My custom title",
           SCpubr::package_report(startup = TRUE,
                       extended = TRUE)        plot.title.face = "italic",
                   plot.subtitle.face = "bold.italic",
                   plot.caption.face = "bold",
                   axis.title.face = "italic",
                   axis.text.face = "plain",
                   legend.title.face = "italic",
                   legend.text.face = "bold.italic",
                   font.type = "mono",
                   font.size = 15)
