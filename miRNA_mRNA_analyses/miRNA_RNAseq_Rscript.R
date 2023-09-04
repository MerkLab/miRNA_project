library(org.Hs.eg.db)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(grDevices)
library(tidyverse)
library(dplyr)
library(ggvenn)
library(DEGreport)
library(smplot2)
library(clusterProfiler)
library(enrichplot)
library(forcats)
library(ggpubr)
library(SuperExactTest)
library(VennDiagram)


#get TPM data and convert ENSG to gene symbol
getwd()
data_TPM <- read.table("data/gene_counts_TPM_mod.txt", header=T, row.names=1)
ENSG.vec = as.vector(row.names(data_TPM))
annots = select(org.Hs.eg.db, keys = ENSG.vec, columns = "SYMBOL", keytype = "ENSEMBL")
write.table(annots, file = "results/ENSG_gene_map.txt")
data_TPM = read.table("data/gene_counts_symbol_TPM_mod.txt", header=T)
data_TPM <- aggregate(data_TPM[-1], by = list(data_TPM$gene), FUN=sum)
row.names(data_TPM) <- data_TPM$Group.1
data_TPM <- data_TPM[, -1]
write.table(data_TPM, file = "results/counts_TPM.txt")

####get raw read counts data and meta data for all cells
meta = read.table("meta/meta.txt", header=T, row.names=1)
meta[] = lapply(meta, factor)
str(meta)
meta_sub = meta[, c("cell_line", "condition", "RNA_isolation")]
levels(meta_sub$RNA_isolation) = c("first", "second")

data_raw = read.table("data/raw_reads.txt", header=T, row.names=1)
data_raw <- aggregate(data_raw[-1], by = list(data_raw$gene), FUN=sum)
row.names(data_raw) <- data_raw$Group.1
data_raw = data_raw[,-1]

### Check that sample names match in both files
# Checks if names are equal
all(names(data_raw) %in% rownames(meta))

# Checks if names are same order
all(names(data_raw) == rownames(meta))

## Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = data_raw, colData = meta, design = ~cell_line + condition)
View(counts(dds))
dim(dds)

###Normalization
# inserts size factors (gene length and seq depth) into dds
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

#filter out low/non expressed genes (genes out where less than 3 samples have norm counts greater/equal to 10, before 56.638 genes, after  21.036 genes)
dim(dds)
idx <- rowSums(counts(dds, normalized=TRUE) >= 10) >= 3
dds <- dds[idx,]
dim(dds)

### Normalized data counts
normalized_counts <- counts(dds, normalized=TRUE)
setwd("/Users/lab/Documents/Daniel/ETMR_miRNA/RNAseq/results")
write.table(normalized_counts, file="normalized_counts.txt", sep="\t", quote=F, col.names=NA)

###Transformation prior to QC
rld <- rlog(dds, blind = TRUE)

### Hierarchical clustering
## Extract the rlog matrix from the object
# assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2
rld_mat <- assay(rld)
## Compute pairwise corrrelation values
rld_cor <- cor(rld_mat)
### Plot heatmap
ann_colors = list(cell_line = c(MCF7="blue3", HT29="darkgreen", PC9="darkorchid3", HL60="darkorange2"),
                  condition = c(empty="darkgrey", mir483="goldenrod2", mir663a="firebrick3"),
                  RNA_isolation=c(first = "cyan3", second = "chocolate2"))

setwd("/Users/lab/Documents/Daniel/ETMR_miRNA/RNAseq/results")
cairo_pdf(filename = "Correlation_plot.pdf", width = 9, height = 7)
pheatmap(rld_cor, annotation = meta_sub, annotation_colors = ann_colors)
dev.off()

###PCA
#make PCA plot for any PC
pca <- prcomp(t(rld_mat))
df <- cbind(meta, pca$x)

setwd("/Users/lab/Documents/Daniel/ETMR_miRNA/RNAseq/results")
cairo_pdf(filename = "PC1-2_line_condition.pdf",width = 8, height = 7)
ggplot(df, aes(x=PC1, y=PC2, color = condition, shape = cell_line, 
               size=cell_line))+
  geom_jitter(width=10, height=10) +
  scale_shape_manual(values=c(15, 16, 17,18))+
  scale_color_manual(values=c('darkblue','darkmagenta', 'darkorange2'))+
  scale_size_manual(values = c(5,5,5,6))+
  theme_classic()
dev.off()

cairo_pdf(filename = "PC1-3_line_condition.pdf",width = 8, height = 7)
ggplot(df, aes(x=PC1, y=PC3, color = condition, shape = cell_line, 
               size=cell_line))+
  geom_jitter(width=10, height=10) +
  scale_shape_manual(values=c(15, 16, 17,18))+
  scale_color_manual(values=c('darkblue','darkmagenta', 'darkorange2'))+
  scale_size_manual(values = c(5,5,5,6))+
  theme_classic()
dev.off()


#####determine DEGs per cell line
###MCF7 data
data_MCF7 = data_raw[,1:9]
dim(data_MCF7)
meta_MCF7 = meta[1:9,]
str(meta_MCF7)
### Check that sample names match in both files
# Checks if names are equal
all(names(data_MCF7) %in% rownames(meta_MCF7))
# Checks if names are same order
all(names(data_MCF7) == rownames(meta_MCF7))

## Create DESeq2Dataset object
dds_MCF7 <- DESeqDataSetFromMatrix(countData = data_MCF7, colData = meta_MCF7, design = ~ condition)

# inserts size factors (gene length and seq depth) into dds
dds_MCF7 <- estimateSizeFactors(dds_MCF7)
sizeFactors(dds_MCF7)

#filter out low/non expressed genes (genes out where less than 3 samples have norm counts greater/equal to 10, before 56.638 genes, after  16.234 genes)
dim(dds_MCF7)
idx <- rowSums(counts(dds_MCF7, normalized=TRUE) >= 10) >= 3
dds_MCF7 <- dds_MCF7[idx,]
dim(dds_MCF7)

### Normalized data counts
normalized_counts_MCF7 <- counts(dds_MCF7, normalized=TRUE)
View(normalized_counts_MCF7)
###Transformation prior to QC
rld_MCF7 <- rlog(dds_MCF7, blind = TRUE)
### Plot PCA
plotPCA(rld_MCF7, intgroup="condition")
rld_mat_MCF7 <- assay(rld_MCF7)
pca_MCF7 <- prcomp(t(rld_mat_MCF7))
df_MCF7 <- cbind(meta_MCF7, pca_MCF7$x)
setwd("/Users/lab/Documents/Daniel/ETMR_miRNA/RNAseq/results")
cairo_pdf(filename = "PC1-2_MCF7.pdf",width = 8, height = 7)
ggplot(df_MCF7)+
  geom_point(aes(x=PC1, y=PC2, color = condition, size=4))+
  scale_color_manual(values=c('darkblue','darkmagenta', 'darkorange2'))+
  theme_classic()
dev.off()

### Running DESeq2 for MCF7
##Fill further add slots to dds object
dds_MCF7 <- DESeq(dds_MCF7)
## Coefficients available for testing
resultsNames(dds_MCF7)
#get results
res_table_MCF7_483 <- results(dds_MCF7, name ="condition_mir483_vs_empty")
summary(res_table_MCF7_483,alpha=0.05)
res_table_MCF7_663a <- results(dds_MCF7, name ="condition_mir663a_vs_empty")
summary(res_table_MCF7_663a,alpha=0.05)
res_MCF7_483_df <- data.frame(res_table_MCF7_483)
res_MCF7_663a_df <- data.frame(res_table_MCF7_663a)
write.table(res_MCF7_483_df, file="DEG_MCF7_mir483.txt")
write.table(res_MCF7_663a_df, file="DEG_MCF7_mir663a.txt")

## Set thresholds for more stringency
padj.cutoff <- 0.05
lfc.cutoff <- 1

threshold_MCF7_483 <- res_table_MCF7_483$padj < padj.cutoff & abs(res_table_MCF7_483$log2FoldChange) > lfc.cutoff
length(which(threshold_MCF7_483))
res_table_MCF7_483$threshold <- threshold_MCF7_483

threshold_MCF7_663a <- res_table_MCF7_663a$padj < padj.cutoff & abs(res_table_MCF7_663a$log2FoldChange) > lfc.cutoff
length(which(threshold_MCF7_663a))
res_table_MCF7_663a$threshold <- threshold_MCF7_663a

#show changes in Volcano plot
LFCs_MCF7_483 <- data.frame(res_table_MCF7_483)
LFCs_MCF7_483$gene = row.names(LFCs_MCF7_483)
length(which(LFCs_MCF7_483$threshold == "TRUE" & LFCs_MCF7_483$log2FoldChange > 1))
length(which(LFCs_MCF7_483$threshold == "TRUE" & LFCs_MCF7_483$log2FoldChange < 1))

setwd("/Users/lab/Documents/Daniel/ETMR_miRNA/RNAseq/results")
cairo_pdf(filename = "Volcano_MCF7_483.pdf",width = 5, height = 5)
LFCs_MCF7_483 %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_MCF7_483 %>% filter(threshold == "TRUE" & log2FoldChange > 1), color = "red") +
  geom_point(data = LFCs_MCF7_483 %>% filter(threshold == "TRUE" & log2FoldChange < 1), color = "blue") +
  xlim(c(-3,3)) +
  ggtitle('MCF7_mir483_vs_empty') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()

LFCs_MCF7_663a <- data.frame(res_table_MCF7_663a)
LFCs_MCF7_663a$gene = row.names(LFCs_MCF7_663a)
length(which(LFCs_MCF7_663a$threshold == "TRUE" & LFCs_MCF7_663a$log2FoldChange > 1))
length(which(LFCs_MCF7_663a$threshold == "TRUE" & LFCs_MCF7_663a$log2FoldChange < 1))

setwd("/Users/lab/Documents/Daniel/ETMR_miRNA/RNAseq/results")
cairo_pdf(filename = "Volcano_MCF7_663a.pdf",width = 5, height = 5)
LFCs_MCF7_663a %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_MCF7_663a %>% filter(threshold == "TRUE" & log2FoldChange > 1), color = "red") +
  geom_point(data = LFCs_MCF7_663a %>% filter(threshold == "TRUE" & log2FoldChange < 1), color = "blue") +
  xlim(c(-3,3)) +
  ggtitle('MCF7_mir663a_vs_empty') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()



###HT29 data
data_HT29 = data_raw[,10:18]
dim(data_HT29)
meta_HT29 = meta[10:18,]
str(meta_HT29)
### Check that sample names match in both files
# Checks if names are equal
all(names(data_HT29) %in% rownames(meta_HT29))
# Checks if names are same order
all(names(data_HT29) == rownames(meta_HT29))

## Create DESeq2Dataset object
dds_HT29 <- DESeqDataSetFromMatrix(countData = data_HT29, colData = meta_HT29, design = ~ condition)

# inserts size factors (gene length and seq depth) into dds
dds_HT29 <- estimateSizeFactors(dds_HT29)
sizeFactors(dds_HT29)

#filter out low/non expressed genes (genes out where less than 3 samples have norm counts greater/equal to 10, before 56.638 genes, after  15.671 genes)
dim(dds_HT29)
idx <- rowSums(counts(dds_HT29, normalized=TRUE) >= 10) >= 3
dds_HT29 <- dds_HT29[idx,]
dim(dds_HT29)

### Normalized data counts
normalized_counts_HT29 <- counts(dds_HT29, normalized=TRUE)
View(normalized_counts_HT29)
###Transformation prior to QC
rld_HT29 <- rlog(dds_HT29, blind = TRUE)
### Plot PCA
plotPCA(rld_HT29, intgroup="condition")
rld_mat_HT29 <- assay(rld_HT29)
pca_HT29 <- prcomp(t(rld_mat_HT29))
df_HT29 <- cbind(meta_HT29, pca_HT29$x)
setwd("/Users/lab/Documents/Daniel/ETMR_miRNA/RNAseq/results")
cairo_pdf(filename = "PC1-2_HT29.pdf",width = 8, height = 7)
ggplot(df_HT29)+
  geom_point(aes(x=PC1, y=PC2, color = condition, size=4))+
  scale_color_manual(values=c('darkblue','darkmagenta', 'darkorange2'))+
  theme_classic()
dev.off()

### Running DESeq2 for HT29
##Fill further add slots to dds object
dds_HT29 <- DESeq(dds_HT29)
## Coefficients available for testing
resultsNames(dds_HT29)
#get results
res_table_HT29_483 <- results(dds_HT29, name ="condition_mir483_vs_empty")
summary(res_table_HT29_483,alpha=0.05)
res_table_HT29_663a <- results(dds_HT29, name ="condition_mir663a_vs_empty")
summary(res_table_HT29_663a,alpha=0.05)
res_HT29_483_df <- data.frame(res_table_HT29_483)
res_HT29_663a_df <- data.frame(res_table_HT29_663a)
write.table(res_HT29_483_df, file="DEG_HT29_mir483.txt")
write.table(res_HT29_663a_df, file="DEG_HT29_mir663a.txt")

threshold_HT29_483 <- res_table_HT29_483$padj < padj.cutoff & abs(res_table_HT29_483$log2FoldChange) > lfc.cutoff
length(which(threshold_HT29_483))
res_table_HT29_483$threshold <- threshold_HT29_483

threshold_HT29_663a <- res_table_HT29_663a$padj < padj.cutoff & abs(res_table_HT29_663a$log2FoldChange) > lfc.cutoff
length(which(threshold_HT29_663a))
res_table_HT29_663a$threshold <- threshold_HT29_663a

#show changes in Volcano plot
LFCs_HT29_483 <- data.frame(res_table_HT29_483)
LFCs_HT29_483$gene = row.names(LFCs_HT29_483)
length(which(LFCs_HT29_483$threshold == "TRUE" & LFCs_HT29_483$log2FoldChange > 1))
length(which(LFCs_HT29_483$threshold == "TRUE" & LFCs_HT29_483$log2FoldChange < 1))

setwd("/Users/lab/Documents/Daniel/ETMR_miRNA/RNAseq/results")
cairo_pdf(filename = "Volcano_HT29_483.pdf",width = 5, height = 5)
LFCs_HT29_483 %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_HT29_483 %>% filter(threshold == "TRUE" & log2FoldChange > 1), color = "red") +
  geom_point(data = LFCs_HT29_483 %>% filter(threshold == "TRUE" & log2FoldChange < 1), color = "blue") +
  xlim(c(-3,3)) +
  ggtitle('HT29_mir483_vs_empty') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()

LFCs_HT29_663a <- data.frame(res_table_HT29_663a)
LFCs_HT29_663a$gene = row.names(LFCs_HT29_663a)
length(which(LFCs_HT29_663a$threshold == "TRUE" & LFCs_HT29_663a$log2FoldChange > 1))
length(which(LFCs_HT29_663a$threshold == "TRUE" & LFCs_HT29_663a$log2FoldChange < 1))

setwd("/Users/lab/Documents/Daniel/ETMR_miRNA/RNAseq/results")
cairo_pdf(filename = "Volcano_HT29_663a.pdf",width = 5, height = 5)
LFCs_HT29_663a %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_HT29_663a %>% filter(threshold == "TRUE" & log2FoldChange > 1), color = "red") +
  geom_point(data = LFCs_HT29_663a %>% filter(threshold == "TRUE" & log2FoldChange < 1), color = "blue") +
  xlim(c(-3,3)) +
  ggtitle('HT29_mir663a_vs_empty') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()



###PC9 data
data_PC9 = data_raw[,19:27]
dim(data_PC9)
meta_PC9 = meta[19:27,]
str(meta_PC9)
### Check that sample names match in both files
# Checks if names are equal
all(names(data_PC9) %in% rownames(meta_PC9))
# Checks if names are same order
all(names(data_PC9) == rownames(meta_PC9))

## Create DESeq2Dataset object
dds_PC9 <- DESeqDataSetFromMatrix(countData = data_PC9, colData = meta_PC9, design = ~ condition)

# inserts size factors (gene length and seq depth) into dds
dds_PC9 <- estimateSizeFactors(dds_PC9)
sizeFactors(dds_PC9)

#filter out low/non expressed genes (genes out where less than 3 samples have norm counts greater/equal to 10, before 56.638 genes, after  16.252 genes)
dim(dds_PC9)
idx <- rowSums(counts(dds_PC9, normalized=TRUE) >= 10) >= 3
dds_PC9 <- dds_PC9[idx,]
dim(dds_PC9)

### Normalized data counts
normalized_counts_PC9 <- counts(dds_PC9, normalized=TRUE)
View(normalized_counts_PC9)
###Transformation prior to QC
rld_PC9 <- rlog(dds_PC9, blind = TRUE)
### Plot PCA
plotPCA(rld_PC9, intgroup="condition")
rld_mat_PC9 <- assay(rld_PC9)
pca_PC9 <- prcomp(t(rld_mat_PC9))
df_PC9 <- cbind(meta_PC9, pca_PC9$x)
setwd("/Users/lab/Documents/Daniel/ETMR_miRNA/RNAseq/results")
cairo_pdf(filename = "PC1-2_PC9.pdf",width = 8, height = 7)
ggplot(df_PC9)+
  geom_point(aes(x=PC1, y=PC2, color = condition, size=4))+
  scale_color_manual(values=c('darkblue','darkmagenta', 'darkorange2'))+
  theme_classic()
dev.off()

### Running DESeq2 for PC9
##Fill further add slots to dds object
dds_PC9 <- DESeq(dds_PC9)
## Coefficients available for testing
resultsNames(dds_PC9)
#get results
res_table_PC9_483 <- results(dds_PC9, name ="condition_mir483_vs_empty")
summary(res_table_PC9_483,alpha=0.05)
res_table_PC9_663a <- results(dds_PC9, name ="condition_mir663a_vs_empty")
summary(res_table_PC9_663a,alpha=0.05)
res_PC9_483_df <- data.frame(res_table_PC9_483)
res_PC9_663a_df <- data.frame(res_table_PC9_663a)
write.table(res_PC9_483_df, file="DEG_PC9_mir483.txt")
write.table(res_PC9_663a_df, file="DEG_PC9_mir663a.txt")

threshold_PC9_483 <- res_table_PC9_483$padj < padj.cutoff & abs(res_table_PC9_483$log2FoldChange) > lfc.cutoff
length(which(threshold_PC9_483))
res_table_PC9_483$threshold <- threshold_PC9_483

threshold_PC9_663a <- res_table_PC9_663a$padj < padj.cutoff & abs(res_table_PC9_663a$log2FoldChange) > lfc.cutoff
length(which(threshold_PC9_663a))
res_table_PC9_663a$threshold <- threshold_PC9_663a

#show changes in Volcano plot
LFCs_PC9_483 <- data.frame(res_table_PC9_483)
LFCs_PC9_483$gene = row.names(LFCs_PC9_483)
length(which(LFCs_PC9_483$threshold == "TRUE" & LFCs_PC9_483$log2FoldChange > 1))
length(which(LFCs_PC9_483$threshold == "TRUE" & LFCs_PC9_483$log2FoldChange < 1))

setwd("/Users/lab/Documents/Daniel/ETMR_miRNA/RNAseq/results")
cairo_pdf(filename = "Volcano_PC9_483.pdf",width = 5, height = 5)
LFCs_PC9_483 %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_PC9_483 %>% filter(threshold == "TRUE" & log2FoldChange > 1), color = "red") +
  geom_point(data = LFCs_PC9_483 %>% filter(threshold == "TRUE" & log2FoldChange < 1), color = "blue") +
  xlim(c(-3,3)) +
  ggtitle('PC9_mir483_vs_empty') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()

LFCs_PC9_663a <- data.frame(res_table_PC9_663a)
LFCs_PC9_663a$gene = row.names(LFCs_PC9_663a)
length(which(LFCs_PC9_663a$threshold == "TRUE" & LFCs_PC9_663a$log2FoldChange > 1))
length(which(LFCs_PC9_663a$threshold == "TRUE" & LFCs_PC9_663a$log2FoldChange < 1))

setwd("/Users/lab/Documents/Daniel/ETMR_miRNA/RNAseq/results")
cairo_pdf(filename = "Volcano_PC9_663a.pdf",width = 5, height = 5)
LFCs_PC9_663a %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_PC9_663a %>% filter(threshold == "TRUE" & log2FoldChange > 1), color = "red") +
  geom_point(data = LFCs_PC9_663a %>% filter(threshold == "TRUE" & log2FoldChange < 1), color = "blue") +
  xlim(c(-3,3)) +
  ggtitle('PC9_mir663a_vs_empty') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()




###HL60 data
data_HL60 = data_raw[,28:36]
dim(data_HL60)
meta_HL60 = meta[28:36,]
str(meta_HL60)
### Check that sample names match in both files
# Checks if names are equal
all(names(data_HL60) %in% rownames(meta_HL60))
# Checks if names are same order
all(names(data_HL60) == rownames(meta_HL60))

## Create DESeq2Dataset object
dds_HL60 <- DESeqDataSetFromMatrix(countData = data_HL60, colData = meta_HL60, design = ~ condition)

# inserts size factors (gene length and seq depth) into dds
dds_HL60 <- estimateSizeFactors(dds_HL60)
sizeFactors(dds_HL60)

#filter out low/non expressed genes (genes out where less than 3 samples have norm counts greater/equal to 10, before 56.638 genes, after  14.528 genes)
dim(dds_HL60)
idx <- rowSums(counts(dds_HL60, normalized=TRUE) >= 10) >= 3
dds_HL60 <- dds_HL60[idx,]
dim(dds_HL60)

### Normalized data counts
normalized_counts_HL60 <- counts(dds_HL60, normalized=TRUE)
View(normalized_counts_HL60)
###Transformation prior to QC
rld_HL60 <- rlog(dds_HL60, blind = TRUE)
### Plot PCA
plotPCA(rld_HL60, intgroup="condition")
rld_mat_HL60 <- assay(rld_HL60)
pca_HL60 <- prcomp(t(rld_mat_HL60))
df_HL60 <- cbind(meta_HL60, pca_HL60$x)
setwd("/Users/lab/Documents/Daniel/ETMR_miRNA/RNAseq/results")
cairo_pdf(filename = "PC1-2_HL60.pdf",width = 8, height = 7)
ggplot(df_HL60)+
  geom_point(aes(x=PC1, y=PC2, color = condition, size=4))+
  scale_color_manual(values=c('darkblue','darkmagenta', 'darkorange2'))+
  theme_classic()
dev.off()

### Running DESeq2 for HL60
##Fill further add slots to dds object
dds_HL60 <- DESeq(dds_HL60)
## Coefficients available for testing
resultsNames(dds_HL60)
#get results
res_table_HL60_483 <- results(dds_HL60, name ="condition_mir483_vs_empty")
summary(res_table_HL60_483,alpha=0.05)
res_table_HL60_663a <- results(dds_HL60, name ="condition_mir663a_vs_empty")
summary(res_table_HL60_663a,alpha=0.05)
res_HL60_483_df <- data.frame(res_table_HL60_483)
res_HL60_663a_df <- data.frame(res_table_HL60_663a)
write.table(res_HL60_483_df, file="DEG_HL60_mir483.txt")
write.table(res_HL60_663a_df, file="DEG_HL60_mir663a.txt")

threshold_HL60_483 <- res_table_HL60_483$padj < padj.cutoff & abs(res_table_HL60_483$log2FoldChange) > lfc.cutoff
length(which(threshold_HL60_483))
res_table_HL60_483$threshold <- threshold_HL60_483

threshold_HL60_663a <- res_table_HL60_663a$padj < padj.cutoff & abs(res_table_HL60_663a$log2FoldChange) > lfc.cutoff
length(which(threshold_HL60_663a))
res_table_HL60_663a$threshold <- threshold_HL60_663a

#show changes in Volcano plot
LFCs_HL60_483 <- data.frame(res_table_HL60_483)
LFCs_HL60_483$gene = row.names(LFCs_HL60_483)
length(which(LFCs_HL60_483$threshold == "TRUE" & LFCs_HL60_483$log2FoldChange > 1))
length(which(LFCs_HL60_483$threshold == "TRUE" & LFCs_HL60_483$log2FoldChange < 1))

setwd("/Users/lab/Documents/Daniel/ETMR_miRNA/RNAseq/results")
cairo_pdf(filename = "Volcano_HL60_483.pdf",width = 5, height = 5)
LFCs_HL60_483 %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_HL60_483 %>% filter(threshold == "TRUE" & log2FoldChange > 1), color = "red") +
  geom_point(data = LFCs_HL60_483 %>% filter(threshold == "TRUE" & log2FoldChange < 1), color = "blue") +
  xlim(c(-3,3)) +
  ggtitle('HL60_mir483_vs_empty') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()

LFCs_HL60_663a <- data.frame(res_table_HL60_663a)
LFCs_HL60_663a$gene = row.names(LFCs_HL60_663a)
length(which(LFCs_HL60_663a$threshold == "TRUE" & LFCs_HL60_663a$log2FoldChange > 1))
length(which(LFCs_HL60_663a$threshold == "TRUE" & LFCs_HL60_663a$log2FoldChange < 1))

cairo_pdf(filename = "Volcano_HL60_663a.pdf",width = 5, height = 5)
LFCs_HL60_663a %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_HL60_663a %>% filter(threshold == "TRUE" & log2FoldChange > 1), color = "red") +
  geom_point(data = LFCs_HL60_663a %>% filter(threshold == "TRUE" & log2FoldChange < 1), color = "blue") +
  xlim(c(-3,3)) +
  ggtitle('HL60_mir663a_vs_empty') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()


####check overlap of genes
##check overlap of upregulated genes in mir483 cells
#make vector of upregulated genes for each cell line

MCF7_483_up = LFCs_MCF7_483 %>% filter(threshold == "TRUE" & log2FoldChange > 1)
MCF7_483_up = MCF7_483_up$gene
HT29_483_up = LFCs_HT29_483 %>% filter(threshold == "TRUE" & log2FoldChange > 1)
HT29_483_up = HT29_483_up$gene
PC9_483_up = LFCs_PC9_483 %>% filter(threshold == "TRUE" & log2FoldChange > 1)
PC9_483_up = PC9_483_up$gene
HL60_483_up = LFCs_HL60_483 %>% filter(threshold == "TRUE" & log2FoldChange > 1)
HL60_483_up = HL60_483_up$gene
ggvenn(list(MCF7 = MCF7_483_up, HT29 = HT29_483_up, PC9 = PC9_483_up, HL60 = HL60_483_up), show_percentage = FALSE)

#make vector of downregulated genes for each cell line
MCF7_483_down = LFCs_MCF7_483 %>% filter(threshold == "TRUE" & log2FoldChange < 1)
MCF7_483_down = MCF7_483_down$gene
HT29_483_down = LFCs_HT29_483 %>% filter(threshold == "TRUE" & log2FoldChange < 1)
HT29_483_down = HT29_483_down$gene
PC9_483_down = LFCs_PC9_483 %>% filter(threshold == "TRUE" & log2FoldChange < 1)
PC9_483_down = PC9_483_down$gene
HL60_483_down = LFCs_HL60_483 %>% filter(threshold == "TRUE" & log2FoldChange < 1)
HL60_483_down = HL60_483_down$gene
ggvenn(list(MCF7 = MCF7_483_down, HT29 = HT29_483_down, PC9 = PC9_483_down, HL60 = HL60_483_down), show_percentage = FALSE)

##check overlap of upregulated genes in mir663a cells
#make vector of upregulated genes for each cell line

MCF7_663a_up = LFCs_MCF7_663a %>% filter(threshold == "TRUE" & log2FoldChange > 1)
MCF7_663a_up = MCF7_663a_up$gene
HT29_663a_up = LFCs_HT29_663a %>% filter(threshold == "TRUE" & log2FoldChange > 1)
HT29_663a_up = HT29_663a_up$gene
PC9_663a_up = LFCs_PC9_663a %>% filter(threshold == "TRUE" & log2FoldChange > 1)
PC9_663a_up = PC9_663a_up$gene
HL60_663a_up = LFCs_HL60_663a %>% filter(threshold == "TRUE" & log2FoldChange > 1)
HL60_663a_up = HL60_663a_up$gene
ggvenn(list(MCF7 = MCF7_663a_up, HT29 = HT29_663a_up, PC9 = PC9_663a_up, HL60 = HL60_663a_up), show_percentage = FALSE)

#make vector of downregulated genes for each cell line
MCF7_663a_down = LFCs_MCF7_663a %>% filter(threshold == "TRUE" & log2FoldChange < 1)
MCF7_663a_down = MCF7_663a_down$gene
HT29_663a_down = LFCs_HT29_663a %>% filter(threshold == "TRUE" & log2FoldChange < 1)
HT29_663a_down = HT29_663a_down$gene
PC9_663a_down = LFCs_PC9_663a %>% filter(threshold == "TRUE" & log2FoldChange < 1)
PC9_663a_down = PC9_663a_down$gene
HL60_663a_down = LFCs_HL60_663a %>% filter(threshold == "TRUE" & log2FoldChange < 1)
HL60_663a_down = HL60_663a_down$gene
ggvenn(list(MCF7 = MCF7_663a_down, HT29 = HT29_663a_down, PC9 = PC9_663a_down, HL60 = HL60_663a_down), show_percentage = FALSE)

#####use DESeq2's LRT to check for genes which change in any direction across all factor levels
##perform likelihood ratio test

dds_lrt = DESeq(dds, test = "LRT", reduced = ~ cell_line)
resultsNames(dds_lrt)
res_lrt = results(dds_lrt)
res_lrt
#create tibble for LRT results
res_lrt_tb = res_lrt %>% data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()
write.table(res_lrt_tb, file="gene_overview_LRT_analysis.txt")
#subset to return genes with padj < 0.01
padj.cutoff.stringent = 0.01
sigLRT_genes = res_lrt_tb %>%
  dplyr::filter(padj < padj.cutoff.stringent)
dim(sigLRT_genes)

#get number of sig genes
nrow(sigLRT_genes)

#get rlog values for sigLRT_genes
cluster_rlog = rld_mat[sigLRT_genes$gene,]
dim(cluster_rlog)

#identify clusters of genes with the degPatterns function
clusters = degPatterns(cluster_rlog, metadata = meta, time = "condition", col = "cell_line")
clusters_line_merge = degPatterns(cluster_rlog, metadata = meta, time = "condition", col = NULL)
#get clustered LRT genes, 4 groups based on 3 conditions (empty, 483, and 663a)
head(clusters_line_merge$df)
setwd("/Users/lab/Documents/Daniel/ETMR_miRNA/RNAseq/results")
LRT_G1 = clusters_line_merge$df %>% dplyr::filter(cluster == 1)
write.table(LRT_G1, file = "LRT_cluster_group1.txt")
LRT_G2 = clusters_line_merge$df %>% dplyr::filter(cluster == 2)
write.table(LRT_G2, file = "LRT_cluster_group2.txt")
LRT_G3 = clusters_line_merge$df %>% dplyr::filter(cluster == 3)
write.table(LRT_G3, file = "LRT_cluster_group3.txt")
LRT_G4 = clusters_line_merge$df %>% dplyr::filter(cluster == 4)
write.table(LRT_G4, file = "LRT_cluster_group4.txt")

#make own boxplot
cluster_df =clusters_line_merge$normalized
cluster_df_G1 = cluster_df %>% dplyr::filter(cluster == 1 & genes %in% LRT_G1$genes)
cluster_df_G2 = cluster_df %>% dplyr::filter(cluster == 2 & genes %in% LRT_G2$genes)
cluster_df_G3 = cluster_df %>% dplyr::filter(cluster == 3 & genes %in% LRT_G3$genes)
cluster_df_G4 = cluster_df %>% dplyr::filter(cluster == 4 & genes %in% LRT_G4$genes)
dim(cluster_df_G1)
dim(cluster_df_G2)
dim(cluster_df_G3)
dim(cluster_df_G4)

setwd("/Users/lab/Documents/Daniel/ETMR_miRNA/RNAseq/results")
cairo_pdf(filename = "LRT_degPattern_G1.pdf", width = 3, height = 7)
ggplot(cluster_df_G1, 
  aes(condition, value,fill = condition))+
  scale_fill_manual(values=c('darkblue','darkmagenta', 'darkorange2'))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(), alpha=0.05, aes(color=condition))+
  scale_color_manual(values=c('darkblue','darkmagenta', 'darkorange2'))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

setwd("/Users/lab/Documents/Daniel/ETMR_miRNA/RNAseq/results")
cairo_pdf(filename = "LRT_degPattern_G2.pdf", width = 3, height = 7)
ggplot(cluster_df_G2, 
       aes(condition, value,fill = condition))+
  scale_fill_manual(values=c('darkblue','darkmagenta', 'darkorange2'))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(), alpha=0.05, aes(color=condition))+
  scale_color_manual(values=c('darkblue','darkmagenta', 'darkorange2'))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

setwd("/Users/lab/Documents/Daniel/ETMR_miRNA/RNAseq/results")
cairo_pdf(filename = "LRT_degPattern_G3.pdf", width = 3, height = 7)
ggplot(cluster_df_G3, 
       aes(condition, value,fill = condition))+
  scale_fill_manual(values=c('darkblue','darkmagenta', 'darkorange2'))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(), alpha=0.1, aes(color=condition))+
  scale_color_manual(values=c('darkblue','darkmagenta', 'darkorange2'))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

cairo_pdf(filename = "LRT_degPattern_G4.pdf", width = 3, height = 7)
ggplot(cluster_df_G4, 
       aes(condition, value,fill = condition))+
  scale_fill_manual(values=c('darkblue','darkmagenta', 'darkorange2'))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(), alpha=0.1, aes(color=condition))+
  scale_color_manual(values=c('darkblue','darkmagenta', 'darkorange2'))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()


#project LRT groups on volcano plots
#G1
cairo_pdf(filename = "Volcano_HL60_483_LRT_G1.pdf", width = 4, height = 7)
LFCs_HL60_483 %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_HL60_483 %>% filter(gene %in% LRT_G1$genes), color = "darkmagenta") +
  xlim(c(-3,3)) +
  ggtitle('HL60_483_LRT_G1') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()

cairo_pdf(filename = "Volcano_MCF7_483_LRT_G1.pdf", width = 4, height = 7)
LFCs_MCF7_483 %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_MCF7_483 %>% filter(gene %in% LRT_G1$genes), color = "darkmagenta") +
  xlim(c(-3,3)) +
  ggtitle('MCF7_483_LRT_G1') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()

cairo_pdf(filename = "Volcano_HT29_483_LRT_G1.pdf", width = 4, height = 7)
LFCs_HT29_483 %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_HT29_483 %>% filter(gene %in% LRT_G1$genes), color = "darkmagenta") +
  xlim(c(-3,3)) +
  ggtitle('HT29_483_LRT_G1') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()

cairo_pdf(filename = "Volcano_PC9_483_LRT_G1.pdf", width = 4, height = 7)
LFCs_PC9_483 %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_PC9_483 %>% filter(gene %in% LRT_G1$genes), color = "darkmagenta") +
  xlim(c(-3,3)) +
  ggtitle('PC9_483_LRT_G1') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()

cairo_pdf(filename = "Volcano_PC9_663a_LRT_G1.pdf", width = 4, height = 7)
LFCs_PC9_663a %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_PC9_663a %>% filter(gene %in% LRT_G1$genes), color = "darkorange2") +
  xlim(c(-3,3)) +
  ggtitle('PC9_663a_LRT_G1') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()

cairo_pdf(filename = "Volcano_MCF7_663a_LRT_G1.pdf", width = 4, height = 7)
LFCs_MCF7_663a %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_MCF7_663a %>% filter(gene %in% LRT_G1$genes), color = "darkorange2") +
  xlim(c(-3,3)) +
  ggtitle('MCF7_663a_LRT_G1') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()

cairo_pdf(filename = "Volcano_HT29_663a_LRT_G1.pdf", width = 4, height = 7)
LFCs_HT29_663a %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_HT29_663a %>% filter(gene %in% LRT_G1$genes), color = "darkorange2") +
  xlim(c(-3,3)) +
  ggtitle('HT29_663a_LRT_G1') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()

cairo_pdf(filename = "Volcano_HL60_663a_LRT_G1.pdf", width = 4, height = 7)
LFCs_HL60_663a %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_HL60_663a %>% filter(gene %in% LRT_G1$genes), color = "darkorange2") +
  xlim(c(-3,3)) +
  ggtitle('HL60_663a_LRT_G1') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()


#G2
cairo_pdf(filename = "Volcano_HL60_483_LRT_G2.pdf", width = 4, height = 7)
LFCs_HL60_483 %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_HL60_483 %>% filter(gene %in% LRT_G2$genes), color = "darkmagenta") +
  xlim(c(-3,3)) +
  ggtitle('HL60_483_LRT_G2') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()

cairo_pdf(filename = "Volcano_HT29_483_LRT_G2.pdf", width = 4, height = 7)
LFCs_HT29_483 %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_HT29_483 %>% filter(gene %in% LRT_G2$genes), color = "darkmagenta") +
  xlim(c(-3,3)) +
  ggtitle('HT29_483_LRT_G2') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()

cairo_pdf(filename = "Volcano_PC9_483_LRT_G2.pdf", width = 4, height = 7)
LFCs_PC9_483 %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_PC9_483 %>% filter(gene %in% LRT_G2$genes), color = "darkmagenta") +
  xlim(c(-3,3)) +
  ggtitle('PC9_483_LRT_G2') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()

cairo_pdf(filename = "Volcano_MCF7_483_LRT_G2.pdf", width = 4, height = 7)
LFCs_MCF7_483 %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_MCF7_483 %>% filter(gene %in% LRT_G2$genes), color = "darkmagenta") +
  xlim(c(-3,3)) +
  ggtitle('MCF7_483_LRT_G2') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()

cairo_pdf(filename = "Volcano_MCF7_663a_LRT_G2.pdf", width = 4, height = 7)
LFCs_MCF7_663a %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_MCF7_663a %>% filter(gene %in% LRT_G2$genes), color = "darkmagenta") +
  xlim(c(-3,3)) +
  ggtitle('MCF7_663a_LRT_G2') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()

cairo_pdf(filename = "Volcano_PC9_663a_LRT_G2.pdf", width = 4, height = 7)
LFCs_PC9_663a %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_PC9_663a %>% filter(gene %in% LRT_G2$genes), color = "darkmagenta") +
  xlim(c(-3,3)) +
  ggtitle('PC9_663a_LRT_G2') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()

cairo_pdf(filename = "Volcano_HL60_663a_LRT_G2.pdf", width = 4, height = 7)
LFCs_HL60_663a %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_HL60_663a %>% filter(gene %in% LRT_G2$genes), color = "darkmagenta") +
  xlim(c(-3,3)) +
  ggtitle('HL60_663a_LRT_G2') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()

cairo_pdf(filename = "Volcano_HT29_663a_LRT_G2.pdf", width = 4, height = 7)
LFCs_HT29_663a %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_HT29_663a %>% filter(gene %in% LRT_G2$genes), color = "darkmagenta") +
  xlim(c(-3,3)) +
  ggtitle('HT29_663a_LRT_G2') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()

#G4
cairo_pdf(filename = "Volcano_HT29_483_LRT_G4.pdf", width = 4, height = 7)
LFCs_HT29_483 %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_HT29_483 %>% filter(gene %in% LRT_G4$genes), color = "darkmagenta") +
  xlim(c(-3,3)) +
  ggtitle('HT29_483_LRT_G4') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()

cairo_pdf(filename = "Volcano_HT29_663a_LRT_G4.pdf", width = 4, height = 7)
LFCs_HT29_663a %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_HT29_663a %>% filter(gene %in% LRT_G4$genes), color = "darkorange2") +
  xlim(c(-3,3)) +
  ggtitle('HT29_663a_LRT_G4') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()

cairo_pdf(filename = "Volcano_MCF7_483_LRT_G4.pdf", width = 4, height = 7)
LFCs_MCF7_483 %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_MCF7_483 %>% filter(gene %in% LRT_G4$genes), color = "darkmagenta") +
  xlim(c(-3,3)) +
  ggtitle('MCF7_483_LRT_G4') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()

cairo_pdf(filename = "Volcano_MCF7_663a_LRT_G4.pdf", width = 4, height = 7)
LFCs_MCF7_663a %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_MCF7_663a %>% filter(gene %in% LRT_G4$genes), color = "darkorange2") +
  xlim(c(-3,3)) +
  ylim(c(0,25))+
  ggtitle('MCF7_663a_LRT_G4') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()

cairo_pdf(filename = "Volcano_PC9_483_LRT_G4.pdf", width = 4, height = 7)
LFCs_PC9_483 %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_PC9_483 %>% filter(gene %in% LRT_G4$genes), color = "darkmagenta") +
  xlim(c(-3,3)) +
  ggtitle('PC9_483_LRT_G4') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()

cairo_pdf(filename = "Volcano_PC9_663a_LRT_G4.pdf", width = 4, height = 7)
LFCs_PC9_663a %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_PC9_663a %>% filter(gene %in% LRT_G4$genes), color = "darkorange2") +
  xlim(c(-3,3)) +
  ylim(c(0,25))+
  ggtitle('PC9_663a_LRT_G4') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()

###volcanos all groups at once
getwd()
setwd("/Users/lab/Documents/Daniel/ETMR_miRNA/RNAseq/results")

cairo_pdf(filename = "Volcano_MCF7_483_LRT_G1-4.pdf", width = 7, height = 7)
LFCs_MCF7_483 %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_MCF7_483 %>% filter(gene %in% LRT_G1$genes), size=3,shape=15,color = "magenta3")+
  geom_point(data = LFCs_MCF7_483 %>% filter(gene %in% LRT_G2$genes), size=3,shape=17,color = "dodgerblue1")+
  geom_point(data = LFCs_MCF7_483 %>% filter(gene %in% LRT_G3$genes), size=4,shape=18,color = "darkorange2") +
  geom_point(data = LFCs_MCF7_483 %>% filter(gene %in% LRT_G4$genes), size=4,shape=19,color = "chartreuse2")+
  xlim(c(-3,3)) +
  ylim(c(0,25))+
  ggtitle('MCF7_483_LRT_G1-4') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()


cairo_pdf(filename = "Volcano_MCF7_663a_LRT_G1-4.pdf", width = 7, height = 7)
LFCs_MCF7_663a %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_MCF7_663a %>% filter(gene %in% LRT_G1$genes), size=3,shape=15,color = "magenta3")+
  geom_point(data = LFCs_MCF7_663a %>% filter(gene %in% LRT_G4$genes), size=3,shape=19,color = "chartreuse2")+
  geom_point(data = LFCs_MCF7_663a %>% filter(gene %in% LRT_G2$genes), size=3,shape=17,color= "dodgerblue1")+
  geom_point(data = LFCs_MCF7_663a %>% filter(gene %in% LRT_G3$genes), size=4,shape=18,color = "darkorange2") +
  xlim(c(-3,3)) +
  ylim(c(0,25))+
  ggtitle('MCF7_663a_LRT_G1-4') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()

cairo_pdf(filename = "Volcano_HT29_483_LRT_G1-4.pdf", width = 7, height = 7)
LFCs_HT29_483 %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_HT29_483 %>% filter(gene %in% LRT_G1$genes), size=3,shape=15,color = "magenta3")+
  geom_point(data = LFCs_HT29_483 %>% filter(gene %in% LRT_G2$genes), size=3,shape=17,color = "dodgerblue1")+
  geom_point(data = LFCs_HT29_483 %>% filter(gene %in% LRT_G3$genes), size=4,shape=18, color = "darkorange2") +
  geom_point(data = LFCs_HT29_483 %>% filter(gene %in% LRT_G4$genes), size=3,shape=19,color = "chartreuse2")+
  xlim(c(-3,3)) +
  ylim(c(0,25))+
  ggtitle('HT29_483_LRT_G1-4') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()

cairo_pdf(filename = "Volcano_HT29_663a_LRT_G1-4.pdf", width = 7, height = 7)
LFCs_HT29_663a %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point()+
  geom_point(data = LFCs_HT29_663a %>% filter(gene %in% LRT_G1$genes), shape=15, size=3,color = "magenta3")+
  geom_point(data = LFCs_HT29_663a %>% filter(gene %in% LRT_G2$genes), shape=17,size=3,color = "dodgerblue1")+
  geom_point(data = LFCs_HT29_663a %>% filter(gene %in% LRT_G3$genes), shape=18, size = 4,color = "darkorange2") +
  geom_point(data = LFCs_HT29_663a %>% filter(gene %in% LRT_G4$genes), shape=19,size=3,color = "chartreuse2")+
  xlim(c(-3,3)) +
  ylim(c(0,25))+
  ggtitle('HT29_663a_LRT_G1-4') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()




###functional analysis for G1 to G4 using clusterProfiler

setwd("/Users/lab/Documents/Daniel/ETMR_miRNA/RNAseq/clusterProfiler")
LRT_groups_func = read.table("LRT_groups_for_function.txt", header = T)

LRT_groups_comp_GOBP <- compareCluster(gene~group, data=LRT_groups_func, fun = enricher,
                                       TERM2GENE=gmt_GOBP)
write.table(LRT_groups_comp_GOBP@compareClusterResult, file = "cluster_G1-G4_GOBP.txt")
LRT_group_GOBP = dotplot(LRT_groups_comp_GOBP, x="group",showCategory=5,font.size = 5) + 
  aes(x=fct_relevel(group, c('G1', 'G2', 'G3', 'G4'))) + xlab(NULL) +
  scale_color_gradientn(colours=c("#371ea3", "#46bac2", "#b3eebe"),
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  guides(size = guide_legend(override.aes=list(shape=1))) +
  theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
        panel.grid.major.x = element_blank())


LRT_groups_comp_REACT <- compareCluster(gene~group, data=LRT_groups_func, fun = enricher,
                                       TERM2GENE=gmt_REACTOME)
write.table(LRT_groups_comp_REACT@compareClusterResult, file = "cluster_G1-G4_REACTOME.txt")
LRT_group_REACT = dotplot(LRT_groups_comp_REACT, x="group", showCategory=5,font.size = 5) + 
  aes(x=fct_relevel(group, c('G1', 'G2', 'G3', 'G4'))) + xlab(NULL) +
  scale_color_gradientn(colours=c("#371ea3", "#46bac2", "#b3eebe"),
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  guides(size = guide_legend(override.aes=list(shape=1))) +
  theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
        panel.grid.major.x = element_blank())

LRT_groups_comp_KEGG <- compareCluster(gene~group, data=LRT_groups_func, fun = enricher,
                                        TERM2GENE=gmt_KEGG)
write.table(LRT_groups_comp_KEGG@compareClusterResult, file = "cluster_G1-G4_KEGG.txt")
LRT_group_KEGG = dotplot(LRT_groups_comp_KEGG, x="group", showCategory=5, font.size = 5) + 
  aes(x=fct_relevel(group, c('G1', 'G2', 'G3', 'G4'))) + xlab(NULL) +
  scale_color_gradientn(colours=c("#371ea3", "#46bac2", "#b3eebe"),
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  guides(size = guide_legend(override.aes=list(shape=1))) +
  theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
        panel.grid.major.x = element_blank())

LRT_groups_comp_HALLMARK <- compareCluster(gene~group, data=LRT_groups_func, fun = enricher,
                                       TERM2GENE=gmt_HALLMARK)
write.table(LRT_groups_comp_HALLMARK@compareClusterResult, file = "cluster_G1-G4_HALLMARK.txt")
LRT_group_HALLMARK = dotplot(LRT_groups_comp_HALLMARK, x="group", showCategory=5, font.size = 5) + 
  aes(x=fct_relevel(group, c('G1', 'G2', 'G3', 'G4'))) + xlab(NULL) +
  scale_color_gradientn(colours=c("#371ea3", "#46bac2", "#b3eebe"),
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  guides(size = guide_legend(override.aes=list(shape=1))) +
  theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
        panel.grid.major.x = element_blank())
print(LRT_group_HALLMARK)


ggarrange(LRT_group_HALLMARK, LRT_group_KEGG, LRT_group_GOBP, nrow=1, align=c("h"), common.legend = TRUE)



###get per condition LFCs for correlation studies and GSEA
dds = DESeq(dds)
resultsNames(dds)
res_table_all_483 = results(dds, name ="condition_mir483_vs_empty")
summary(res_table_all_483,alpha=0.05)
res_table_all_663a = results(dds, name ="condition_mir663a_vs_empty")
summary(res_table_all_663a,alpha=0.05)
res_all_483_df <- data.frame(res_table_all_483)
res_all_663a_df <- data.frame(res_table_all_663a)
write.table(res_all_483_df, file="DEG_all_mir483.txt")
write.table(res_all_663a_df, file="DEG_all_mir663a.txt")

#make correlation analysis
cor_483_663a = read.csv("Correlation_DEGs_all_lines_483_vs_663a.csv", header = T)

cairo_pdf(filename = "Correlation_all_lines_483_vs_663a_effect.pdf", width = 5, height = 5)
ggplot(data = cor_483_663a, mapping = aes(x=LFC.483, y=LFC.663a))+
  geom_point(shape=21, fill = "maroon2", color = "white", size=3)+
  xlim(c(-3,4))+
  ylim(c(-3,4))+
  sm_statCorr(color = "black", linetype ="dashed")+
  theme_classic()
dev.off()

####functional analysis of regulated genes using clusterProfiler
#get data
getwd()
setwd("/Users/lab/Documents/Daniel/ETMR_miRNA/RNAseq/results")
MCF7_483_up_ENTREZ = bitr(MCF7_483_up, fromType = "SYMBOL", toType = "ENTREZID",
                          OrgDb = "org.Hs.eg.db")
write.table(MCF7_483_up_ENTREZ, file = "MCF7_483_up_ENTREZ.txt")
MCF7_483_down_ENTREZ = bitr(MCF7_483_down, fromType = "SYMBOL", toType = "ENTREZID",
                          OrgDb = "org.Hs.eg.db")
write.table(MCF7_483_down_ENTREZ, file = "MCF7_483_down_ENTREZ.txt")

MCF7_663a_up_ENTREZ = bitr(MCF7_663a_up, fromType = "SYMBOL", toType = "ENTREZID",
                          OrgDb = "org.Hs.eg.db")
write.table(MCF7_663a_up_ENTREZ, file = "MCF7_663a_up_ENTREZ.txt")
MCF7_663a_down_ENTREZ = bitr(MCF7_663a_down, fromType = "SYMBOL", toType = "ENTREZID",
                            OrgDb = "org.Hs.eg.db")
write.table(MCF7_663a_down_ENTREZ, file = "MCF7_663a_down_ENTREZ.txt")


HL60_483_up_ENTREZ = bitr(HL60_483_up, fromType = "SYMBOL", toType = "ENTREZID",
                          OrgDb = "org.Hs.eg.db")
write.table(HL60_483_up_ENTREZ, file = "HL60_483_up_ENTREZ.txt")
MCF7_HL60_down_ENTREZ = bitr(HL60_483_down, fromType = "SYMBOL", toType = "ENTREZID",
                            OrgDb = "org.Hs.eg.db")
write.table(HL60_483_down_ENTREZ, file = "HL60_483_down_ENTREZ.txt")

HL60_663a_up_ENTREZ = bitr(HL60_663a_up, fromType = "SYMBOL", toType = "ENTREZID",
                           OrgDb = "org.Hs.eg.db")
write.table(HL60_663a_up_ENTREZ, file = "HL60_663a_up_ENTREZ.txt")
HL60_663a_down_ENTREZ = bitr(HL60_663a_down, fromType = "SYMBOL", toType = "ENTREZID",
                             OrgDb = "org.Hs.eg.db")
write.table(HL60_663a_down_ENTREZ, file = "HL60_663a_down_ENTREZ.txt")


PC9_483_up_ENTREZ = bitr(PC9_483_up, fromType = "SYMBOL", toType = "ENTREZID",
                          OrgDb = "org.Hs.eg.db")
write.table(PC9_483_up_ENTREZ, file = "PC9_483_up_ENTREZ.txt")
PC9_483_down_ENTREZ = bitr(PC9_483_down, fromType = "SYMBOL", toType = "ENTREZID",
                             OrgDb = "org.Hs.eg.db")
write.table(PC9_483_down_ENTREZ, file = "PC9_483_down_ENTREZ.txt")

PC9_663a_up_ENTREZ = bitr(PC9_663a_up, fromType = "SYMBOL", toType = "ENTREZID",
                           OrgDb = "org.Hs.eg.db")
write.table(PC9_663a_up_ENTREZ, file = "PC9_663a_up_ENTREZ.txt")
PC9_663a_down_ENTREZ = bitr(PC9_663a_down, fromType = "SYMBOL", toType = "ENTREZID",
                             OrgDb = "org.Hs.eg.db")
write.table(PC9_663a_down_ENTREZ, file = "PC9_663a_down_ENTREZ.txt")


HT29_483_up_ENTREZ = bitr(HT29_483_up, fromType = "SYMBOL", toType = "ENTREZID",
                         OrgDb = "org.Hs.eg.db")
write.table(HT29_483_up_ENTREZ, file = "HT29_483_up_ENTREZ.txt")
HT29_483_down_ENTREZ = bitr(HT29_483_down, fromType = "SYMBOL", toType = "ENTREZID",
                           OrgDb = "org.Hs.eg.db")
write.table(HT29_483_down_ENTREZ, file = "HT29_483_down_ENTREZ.txt")

HT29_663a_up_ENTREZ = bitr(HT29_663a_up, fromType = "SYMBOL", toType = "ENTREZID",
                          OrgDb = "org.Hs.eg.db")
write.table(HT29_663a_up_ENTREZ, file = "HT29_663a_up_ENTREZ.txt")
HT29_663a_down_ENTREZ = bitr(HT29_663a_down, fromType = "SYMBOL", toType = "ENTREZID",
                            OrgDb = "org.Hs.eg.db")
write.table(HT29_663a_down_ENTREZ, file = "HT29_663a_down_ENTREZ.txt")

#get dataframe with DEG upregulated in all cell lines with KO of 483 and 663a, and gene sets
getwd()
setwd("/Users/lab/Documents/Daniel/ETMR_miRNA/RNAseq/clusterProfiler")
DEG_all_lines_483_663a_up = read.table("DEG_ENTREZ_all_lines_483_663a_up.txt", header = TRUE)

gmt_GOBP = read.gmt("c5.go.bp.v2023.1.Hs.entrez.gmt")
gmt_REACTOME = read.gmt("c2.cp.reactome.v2023.1.Hs.entrez.gmt")
gmt_BIOCARTA = read.gmt("c2.cp.biocarta.v2023.1.Hs.entrez.gmt")
gmt_KEGG = read.gmt("c2.cp.kegg.v2023.1.Hs.entrez.gmt")
gmt_HALLMARK = read.gmt("h.all.v2023.1.Hs.entrez.gmt")
gmt_GOCC = read.gmt("c5.go.cc.v2023.1.Hs.entrez.gmt")
gmt_GOMF = read.gmt("c5.go.mf.v2023.1.Hs.entrez.gmt")

setwd("/Users/lab/Documents/Daniel/ETMR_miRNA/RNAseq/clusterProfiler")
cluster_all_lines_up_HALLMARK <- compareCluster(Gene~cell_line+condition, data=DEG_all_lines_483_663a_up, fun = enricher,
                     TERM2GENE=gmt_HALLMARK)
write.table(cluster_all_lines_up_HALLMARK@compareClusterResult, file = "cluster_483_663_up_HALLMARK_1.txt")
cluster_plot_all_lines_up_HALLMARK <- dotplot(cluster_all_lines_up_HALLMARK, x="cell_line",
                                              showCategory=3) + facet_grid(~condition) +
  aes(x=fct_relevel(cell_line, c('HL60', 'HT29', 'MCF7', 'PC9'))) + xlab(NULL) +
  scale_color_gradientn(colours=c("#371ea3", "#46bac2", "#b3eebe"),
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  guides(size = guide_legend(override.aes=list(shape=1))) +
  theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
        panel.grid.major.x = element_blank())
print(cluster_plot_all_lines_up_HALLMARK)
cluster_plot_all_lines_up_HALLMARK <- dotplot(cluster_all_lines_up_HALLMARK, x="cell_line",
                                              showCategory=c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
                                                             "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                                                             "HALLMARK_INFLAMMATORY_RESPONSE",
                                                             "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                                                             "HALLMARK_APOPTOSIS",
                                                             "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                                                             "HALLMARK_HYPOXIA",
                                                             "HALLMARK_P53_PATHWAY")) + facet_grid(~condition) +
  aes(x=fct_relevel(cell_line, c('HL60', 'HT29', 'MCF7', 'PC9'))) + xlab(NULL) +
  scale_color_gradientn(colours=c("#371ea3", "#46bac2", "#b3eebe"),
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  guides(size = guide_legend(override.aes=list(shape=1))) +
  theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
        panel.grid.major.x = element_blank())
print(cluster_plot_all_lines_up_HALLMARK)


cluster_all_lines_up_KEGG <- compareCluster(Gene~cell_line+condition, data=DEG_all_lines_483_663a_up, fun = enricher,
                                                TERM2GENE=gmt_KEGG)
write.table(cluster_all_lines_up_KEGG@compareClusterResult, file = "cluster_483_663_up_KEGG_1.txt")
cluster_plot_all_lines_up_KEGG <- dotplot(cluster_all_lines_up_KEGG, x="cell_line",
                                              showCategory=3) + facet_grid(~condition) +
  aes(x=fct_relevel(cell_line, c('HL60', 'HT29', 'MCF7', 'PC9'))) + xlab(NULL) +
  scale_color_gradientn(colours=c("#371ea3", "#46bac2", "#b3eebe"),
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  guides(size = guide_legend(override.aes=list(shape=1))) +
  theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
        panel.grid.major.x = element_blank())
print(cluster_plot_all_lines_up_KEGG)

cluster_all_lines_up_REACT <- compareCluster(Gene~cell_line+condition, data=DEG_all_lines_483_663a_up, fun = enricher,
                                            TERM2GENE=gmt_REACTOME)
write.table(cluster_all_lines_up_REACT@compareClusterResult, file = "cluster_483_663_up_REACTOME_1.txt")
cluster_plot_all_lines_up_REACT <- dotplot(cluster_all_lines_up_REACT, x="cell_line",
                                           showCategory=c("REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS",
                                                          "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
                                                          "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
                                                          "REACTOME_INTERFERON_SIGNALING",
                                                          "REACTOME_SIGNALING_BY_INTERLEUKINS",
                                                          "REACTOME_INTERFERON_GAMMA_SIGNALING",
                                                          "REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES",
                                                          "REACTOME_INTERLEUKIN_10_SIGNALING")) + facet_grid(~condition) +
  aes(x=fct_relevel(cell_line, c('HL60', 'HT29', 'MCF7', 'PC9'))) + xlab(NULL) +
  scale_color_gradientn(colours=c("#371ea3", "#46bac2", "#b3eebe"),
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  guides(size = guide_legend(override.aes=list(shape=1))) +
  theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
        panel.grid.major.x = element_blank())
print(cluster_plot_all_lines_up_REACT)

ggarrange(cluster_plot_all_lines_up_HALLMARK, cluster_plot_all_lines_up_REACT, nrow=1, align=c("h"), common.legend = TRUE)



###GSEA analysis
#get data
getwd()
LFC_483_all_ENTREZ = bitr(row.names(res_all_483_df), fromType = "SYMBOL", toType = "ENTREZID",
                          OrgDb = "org.Hs.eg.db")
write.table(LFC_483_all_ENTREZ, file="LFCs_ENTREZ_IDs.txt")

LFC_all_483_rank = read.table("LFC_ENTREZ_all_483.txt", header = T)
LFC_all_663a_rank = read.table("LFC_ENTREZ_all_663a.txt", header = T)
LFC_all_483_rank_vec = as.numeric(LFC_all_483_rank$log2FoldChange)
names(LFC_all_483_rank_vec) = LFC_all_483_rank$ENTREZID
LFC_all_663a_rank_vec = as.numeric(LFC_all_663a_rank$log2FoldChange)
names(LFC_all_663a_rank_vec) = LFC_all_663a_rank$ENTREZID

#KEGG GSEA
KEGG_all_483 = gseKEGG(geneList = LFC_all_483_rank_vec, organism = "hsa")
KEGG_all_663a = gseKEGG(geneList = LFC_all_663a_rank_vec, organism = "hsa")
write.table(KEGG_all_483@result, file = "KEGG_GSEA_all_483.txt")
write.table(KEGG_all_663a@result, file = "KEGG_GSEA_all_663a.txt")

## sorted by absolute values of NES
KEGG_all_483_sort <- arrange(KEGG_all_483, desc(abs(NES)))
KEGG_all_663a_sort <- arrange(KEGG_all_663a, desc(abs(NES)))

#visualize distinct sets
color <- c("#f7ca64", "#43a5bf", "#86c697", "#a670d6", "#ef998a")

#show depleted KEGG gene sets
gseaplot2(KEGG_all_483_sort, geneSetID = c("hsa00970", "hsa03030", "hsa00020", "hsa03020"),
          color = c("#000066", "#0000CC", "#3366FF", "#66CCFF"),subplots = 1:2, pvalue_table=F, base_size=20)
gseaplot2(KEGG_all_663a_sort, geneSetID = c("hsa00970", "hsa03030", "hsa00020", "hsa00620"),
          color = c("#000066", "#0000CC", "#3366FF", "#66CCFF"),subplots = 1:2, pvalue_table=F, base_size=20)

#show enriched KEGG gene sets
gseaplot2(KEGG_all_483_sort, geneSetID = c("hsa04060","hsa04668", "hsa04064" , "hsa04514"),
          color = c("#990000", "#CC3333", "#FF3333", "#FF0066"),subplots = 1:2, pvalue_table=F, base_size=20)
gseaplot2(KEGG_all_663a_sort, geneSetID = c("hsa04514","hsa04064", "hsa04060" , "hsa04668"),
          color = c("#990000", "#CC3333", "#FF3333", "#FF0066"),subplots = 1:2, pvalue_table=F, base_size=20)


















