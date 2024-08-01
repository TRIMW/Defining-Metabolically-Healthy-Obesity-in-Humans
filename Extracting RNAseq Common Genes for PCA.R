library(dlpyr)
library(writexl)

##Combining Bulk RNAseq datasets for PCA analysis
##Step 1 load in data
df1 <- read.csv("/Users/willtrim/Desktop/Liver.csv")
df2 <- read.csv("/Users/willtrim/Desktop/Omentum.csv")
df3 <- read.csv("/Users/willtrim/Desktop/Subq.csv")

##Extract row names (genes)
row_names1 <- df1$Gene
row_names2 <- df2$Gene
row_names3 <- df3$Gene

##Create data.frame of all data across tissues with only gene names consistent across all
common_data <- inner_join(df1, inner_join(df2, df3, by = "Gene"), by = "Gene")

##remove duplicate gene names
common_data <- common_data %>% distinct(Gene, .keep_all = TRUE)