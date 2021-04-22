#install.packages("BiocManager")
#BiocManager::install("DESeq2")
#BiocManager::install("GenomeInfoDb")
#BiocManager::install("latticeExtra")
#library(DESeq2)

getwd()
require(DESeq2)
setwd("F:/临时/class/")
data <- read.table("P7vsPAO1.txt", row.names = 1,header = T, na.strings = c("","NA"), skipNul=TRUE)
dim(data)
colnames(data)
rownames(data)
nrow(data)
colSums(data)
Sample <- c("PAO1", "PAO1", "P7", "P7")

samples <- data.frame(row.names=colnames(data), Group=as.factor(Sample))

DS_Table <- DESeqDataSetFromMatrix(countData = data, colData=samples, design=~Group)

rowSums(counts(DS_Table))

DS_Table_sort <- DS_Table[ rowSums(counts(DS_Table)) > 1, ]

dim(DS_Table)
dim(DS_Table_sort)

DS_Table_sort <- estimateSizeFactors(DS_Table_sort)

normalized_counts <- counts(DS_Table_sort, normalized=TRUE)
write.table(normalized_counts, file="PAO1_P7_2samples_normalized.csv")
DS <- DESeq(DS_Table_sort)

# testing two transformation and one non-transformed
rld <- rlogTransformation(DS, blind=TRUE) # transformation
vsd <- varianceStabilizingTransformation(DS, blind=TRUE) # transformation
nt <- normTransform(DS) # non transform


comparison11 <- results(DS, contrast=c("Group","P7","PAO1"))# ,alpha=0.05)
summary(comparison11)
comparison11 <- subset(comparison11, padj < 0.05)
comparison11 <- comparison11[abs(comparison11$log2FoldChange) >1,]
comparison11_df <- as.data.frame(comparison11)
head(comparison11_df)
dim(comparison11_df)
write.table(comparison11_df, file="P7 vs PAO1.csv", sep=",")

#### 
library(pheatmap)
pheatmap(assay(nt),
         kmeans_k = NA, breaks = NA, border_color = "white",
         cellwidth = NA, cellheight = NA, scale = "none", cluster_rows = TRUE,
         cluster_cols = F, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", clustering_method = "average",
         cutree_rows = NA, cutree_cols = NA,
         legend = TRUE, legend_breaks = NA,
         legend_labels = NA, annotation_row = NA,
         annotation = NA, annotation_colors = NA, annotation_legend = TRUE,
         annotation_names_row = T, annotation_names_col = TRUE,
         drop_levels = TRUE, show_rownames = F, show_colnames = T, main = NA,
         fontsize = 10, fontsize_row = 4.5, fontsize_col = 10, display_numbers = F,
         gaps_row = NULL, gaps_col = NULL, labels_row = NULL,
         labels_col = NULL, filename = NA, width = NA, height = NA,
         silent = FALSE, na_col = "#DDDDDD") #color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(256),
dev.off() ###sometimes PCA command gives errors because previous command might have overloaded the graphics. Here the heatmap gives errors as its too big. thats why run this command and then plot the PCA###

plotPCA((nt),intgroup = 'Group')

library(ggplot2)
library(ggrepel)
p = plotPCA((nt),intgroup = 'Group')
p <- p + theme(legend.position = 'none') + geom_point(size = 2) +geom_text_repel(aes_string(label = "name"), size = 5)
print(p)
dev.off()