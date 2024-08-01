library(dplyr)
library(ggplot2)
library(Seurat)
library(tidyverse)

df <- read.csv('/Users/johnscanlan/Documents/PhD/RNA-SEQ analsyis/MSCs_Inferred_TF_activity_A5.csv')
sc <- readRDS('/Users/johnscanlan/Documents/PhD/Will_BITFAM/OM_merge.rds')

# Subset the Seurat object
clusters_to_include <- c("MSC1s", "MSC2s", "MSC3s")
sc <- subset(x = sc, ident = clusters_to_include)

# Assuming df is your DataFrame
rownames(df) <- df$barcode
tf_mat <- df[, !(names(df) %in% c("clusters", "health", "barcode"))]
# Convert tf_mat to a numeric matrix
tf_mat <- as.data.frame(tf_mat)
sc@meta.data$barcode <- rownames(sc@meta.data)

tf_mat <- t(tf_mat)

tf_mat <- CreateSeuratObject(tf_mat, project = 'tfs')
tf_mat@meta.data <- sc@meta.data
tf_mat@reductions <- sc@reductions
tf_mat@graphs <- sc@graphs
tf_mat@neighbors <- sc@neighbors
tf_mat@images <- sc@images

DefaultAssay(object = MSC_BIT) <- 'RNA'
MSC_BIT <- NormalizeData(object = MSC_BIT, assay = 'RNA')
MSC_BIT <- ScaleData(MSC_BIT, assay = 'RNA')

DimPlot(tf_mat, label = T, group.by = 'seurat_clusters')
DimPlot(sc, label = T)

FeaturePlot(tf_mat, features = 'ESR1')

clus3 <- FindMarkers(tf_mat, ident.1 = '3', ident.2 = c('6', '15'), group.by = 'seurat_clusters')
clus6 <- FindMarkers(tf_mat, ident.1 = '6', ident.2 = c('3', '15'), group.by = 'seurat_clusters')
clus15 <- FindMarkers(tf_mat, ident.1 = '15', ident.2 = c('3', '6'), group.by = 'seurat_clusters')
h_vs_u <- FindMarkers(tf_mat, ident.1 = 'OM_MHO', ident.2 = c('OM_MUO'), group.by = 'Sample')

View(h_vs_u)

# Load the required libraries if not already loaded
# install.packages("ggplot2")
# library(ggplot2)

# Convert 'clusters' to a factor with labels '3', '6', and '15'
df$clusters <- factor(df$clusters, levels = c(3, 6, 15), labels = c('3', '6', '15'))

# Replace numeric values with strings in the 'clusters' column of 'df'
df <- df %>%
  mutate(clusters = case_when(
    clusters == 3 ~ "MSC2s",
    clusters == 6 ~ "MSC1s",
    clusters == 15 ~ "MSC3s",
    TRUE ~ as.character(clusters)  # Keep other values as is
  ))


# Create separate boxplots for each cell type
ggplot(df, aes(x = clusters, y = MAF, group = clusters)) +
  geom_boxplot() +
  labs(x = "Cell Types", y = "Expression of MAF") +
  ggtitle("Boxplot of MAF Expression by Cell Type")

DimPlot(sc, split.by  = 'seurat_clusters')
FeaturePlot(sc, features = 'MAF')

# Create a publication-quality boxplot with no gridlines
ggplot(df, aes(x = clusters, y = MAF, group = clusters)) +
  geom_boxplot(width = 0.7, fill = "lightblue", color = "darkblue", alpha = 0.7) +
  labs(x = "Cell Types", y = "Expression of MAF") +
  ggtitle("Boxplot of MAF Expression by Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(size = 14, face = "bold")) +
  theme(axis.title = element_text(size = 12, face = "bold")) +
  theme(axis.text = element_text(size = 10)) +
  theme(legend.position = "none") +
  theme(panel.grid = element_blank())  # Remove gridlines

# Load the ggplot2 library if not already loaded
library(ggplot2)

# Load the ggplot2 library if not already loaded
library(ggplot2)

# Create a publication-quality boxplot with x and y-axis labels, and different colors for gene groups
ggplot(df, aes(x = clusters, y = MAF, fill = gene_group)) +
  geom_boxplot(width = 0.7, alpha = 0.7) +
  labs(x = "Cell Types", y = "Expression of MAF") +
  xlab("Cell Types") +
  ylab("Expression of MAF") +
  ggtitle("Boxplot of MAF Expression by Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(size = 14, face = "bold")) +
  theme(axis.title = element_text(size = 12, face = "bold")) +
  theme(axis.text = element_text(size = 10)) +
  theme(legend.position = "top") +  # Move legend to the top
  scale_fill_manual(values = c("GroupA" = "lightblue", "GroupB" = "lightgreen", "GroupC" = "lightpink"))

DotPlot(tf_mat, features = c('MAF', 'AR', 'JUN', 'FOXC1', 'ESR1', 'MYB'), group.by = 'seurat_clusters', split.by = 'Sample') + RotatedAxis()

# Calculate column averages
averages <- colMeans(df)

# Create a new data frame with column names and averages
averages_df <- data.frame(
  ColumnName = names(averages),
  Average = averages
)

# Print the new data frame
print(averages_df)


VlnPlot(tf_mat, features = c('MAF', 'AR', 'JUN', 'FOXC1', 'ESR1', 'MYB'), stack = T, flip = T,adjust = T, group.by = 'Sample')

View(df)

#Most active TFs

# Calculate the average expression for each column
df <- df[, !(names(df) %in% c("clusters", "health", "barcode"))]
average_expression <- colMeans(df)

# Create a data frame with column names and their average expression
average_expression_df <- data.frame(TF = names(average_expression), Average_Expression = average_expression)

# Sort the data frame by average expression in descending order
sorted_df <- average_expression_df[order(average_expression_df$Average_Expression, decreasing = TRUE), ]

# Print the sorted data frame (top N columns with highest average expression)
top_n <- 10  # Replace with your desired number of top columns
top_columns <- sorted_df[1:top_n, ]

print("Top Columns with Highest Average Expression:")
print(top_columns)

View(clus3)
VlnPlot(tf_mat, features = c('FOXC1', 'ESR1', 'RUNX1', 'JUN', 'AR', 'NR3C1', 'MAF', 'MYB', 'FOS'), stack = T, flip = T, group.by = 'seurat_clusters')


# Save clus3 with row names
write.csv(clus3, '/Users/johnscanlan/Documents/PhD/Will_BITFAM/MSC2s_markers.csv', row.names = TRUE)

# Save clus6 with row names
write.csv(clus6, '/Users/johnscanlan/Documents/PhD/Will_BITFAM/MSC1s_markers.csv', row.names = TRUE)

# Save clus15 with row names
write.csv(clus15, '/Users/johnscanlan/Documents/PhD/Will_BITFAM/MSC3s_markers.csv', row.names = TRUE)

h <- FindMarkers(tf_mat, ident.1 = 'OM_MUO', ident.2 = 'OM_MHO', group.by = 'Sample')
View(h)

VlnPlot(tf_mat, features = 'PPARG', pt.size = 0.1, group.by = 'Sample')

saveRDS(tf_mat, '/Users/johnscanlan/Documents/PhD/Will_BITFAM/tf_mat.rds')
