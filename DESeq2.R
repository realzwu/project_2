setwd("/Users/zd821/Documents/project2/mRNA")
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("edgeR")
# BiocManager::install("DESeq2")

library(openxlsx)
library(DESeq2)
library(ggplot2)

counts <- read.xlsx("featureCounts_sorted.xlsx", colNames = T, rowNames = T)

# PCA plots

colData0 <- data.frame(row.names = c("RNS11A","RNS11B","RNS11C","RPCSTCA","RPCSTCB","RPCSTCC","RPCSTDA","RPCSTDB","RPCSTDC","RPSM3A","RPSM3B","RPSM3C","RSH6A","RSH6B","RSH6C","RWSH8A","RWSH8B","RWSH8C"),
                       condition = factor(c("RWPE1 control","RWPE1 control","RWPE1 control","PC3 control","PC3 control","PC3 control","PC3 overexpressed","PC3 overexpressed","PC3 overexpressed","WPMY1 control","WPMY1 control","WPMY1 control","RWPE1 silenced","RWPE1 silenced","RWPE1 silenced","WPMY1 silenced","WPMY1 silenced","WPMY1 silenced"),levels = c("RWPE1 silenced","RWPE1 control","PC3 overexpressed","PC3 control","WPMY1 silenced","WPMY1 control")))

dds0 = DESeqDataSetFromMatrix(countData = counts[c(1:9,16:24)], colData = colData0, design = ~ condition)
dds0 = DESeq(dds0)
vsd <- vst(dds0, blind=FALSE)
pcaData = plotPCA(vsd, intgroup="condition", returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=2,alpha = 0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() 

ggsave("TCGA_survival_DKK3.png",width=6.5,heigh=4.8,dpi=600)

### Differentially expressed genes in RWPE1

colData <- data.frame(row.names = c("RNS11A","RNS11B","RNS11C","RSH6A","RSH6B","RSH6C"), condition = factor(c("control","control","control","case","case","case"),levels = c("control","case")))

dds = DESeqDataSetFromMatrix(countData = counts[c(1:3,19:21)], colData = colData, design = ~ condition)
dds = DESeq(dds)
sizeFactors(dds)
res = results(dds)
class(res)
res <- as.data.frame(res)
res = cbind(rownames(res), res)
colnames(res) = c("gene_id","baseMean","log2FoldChange","lfcSE","stat","pval","padj")

write.xlsx(res, "DESeq2-all.gene-epi.xlsx")

resSig <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1),]

resSig[which(resSig$log2FoldChange > 0), "up_down"] <- "Up"
resSig[which(resSig$log2FoldChange < 0), "up_down"] <- "Down"

write.xlsx(resSig, "DESeq2-Pval0.05-FC2-epi.xlsx")


# Differentially expressed genes in WPMY1

colData <- data.frame(row.names = c("RPSM3A","RPSM3B","RPSM3C","RWSH8A","RWSH8B","RWSH8C"), condition = factor(c("control","control","control","case","case","case"),levels = c("control","case")))

dds = DESeqDataSetFromMatrix(countData = counts[c(16:18,22:24)], colData = colData, design = ~ condition)
dds = DESeq(dds)
sizeFactors(dds)
res = results(dds)
class(res)
res <- as.data.frame(res)
res = cbind(rownames(res), res)
colnames(res) = c("gene_id","baseMean","log2FoldChange","lfcSE","stat","pval","padj")

write.xlsx(res, "DESeq2-all.gene-str.xlsx")

resSig <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1),]

resSig[which(resSig$log2FoldChange > 0), "up_down"] <- "Up"
resSig[which(resSig$log2FoldChange < 0), "up_down"] <- "Down"

write.xlsx(resSig, "DESeq2-Pval0.05-FC2-str.xlsx")


# Differentially expressed genes in pc3

colData <- data.frame(row.names = c("RPCSTCA","RPCSTCB","RPCSTCC","RPCSTDA","RPCSTDB","RPCSTDC"), condition = factor(c("control","control","control","case","case","case"),levels = c("control","case")))

dds = DESeqDataSetFromMatrix(countData = counts[c(4:9)], colData = colData, design = ~ condition)
dds = DESeq(dds)
sizeFactors(dds)
res = results(dds)
class(res)
res <- as.data.frame(res)
res = cbind(rownames(res), res)
colnames(res) = c("gene_id","baseMean","log2FoldChange","lfcSE","stat","pval","padj")

write.xlsx(res, "DESeq2-all-pc3.gene.xlsx")

resSig <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1),]

resSig[which(resSig$log2FoldChange > 0), "up_down"] <- "Up"
resSig[which(resSig$log2FoldChange < 0), "up_down"] <- "Down"

write.xlsx(resSig, "DESeq2-Pval0.05-FC2-pc3.xlsx")

### Volcano plots

pc3 <- read.xlsx("Stable/DESeq2-all-pc3.gene.xlsx", colNames = T, rowNames = T)
str <- read.xlsx("Stromal/DESeq2-all.gene-str.xlsx", colNames = T, rowNames = T)
epi <- read.xlsx("Epithelial/DESeq2-all.gene-epi.xlsx", colNames = T, rowNames = T)

pc3 = na.omit(pc3)
str = na.omit(str)
epi = na.omit(epi)

pc3[which(pc3$log2FoldChange > 1), "Expression"] <- "Upregulated"
pc3[which(pc3$log2FoldChange < -1), "Expression"] <- "Downregulated"
str[which(str$log2FoldChange > 1), "Expression"] <- "Upregulated"
str[which(str$log2FoldChange < -1), "Expression"] <- "Downregulated"
epi[which(epi$log2FoldChange > 1), "Expression"] <- "Upregulated"
epi[which(epi$log2FoldChange < -1), "Expression"] <- "Downregulated"

pc3[which(pc3$log2FoldChange > 10), "log2FoldChange"] <- 10
pc3[which(pc3$log2FoldChange < -10), "log2FoldChange"] <- -10
pc3[which(pc3$padj < 1e-150), "padj"] <- 1e-150

str[which(str$log2FoldChange > 10), "log2FoldChange"] <- 10
str[which(str$log2FoldChange < -10), "log2FoldChange"] <- -10
str[which(str$padj < 1e-150), "padj"] <- 1e-150

epi[which(epi$log2FoldChange > 10), "log2FoldChange"] <- 10
epi[which(epi$log2FoldChange < -10), "log2FoldChange"] <- -10
epi[which(epi$padj < 1e-150), "padj"] <- 1e-150

ggplot(data=pc3, aes(x=log2FoldChange, y=-log10(padj), col=Expression)) + 
  geom_point(alpha = 0.7) + xlim(-10,10) + ylim(0,150) +
  theme_minimal()

ggsave("PC3_volcano.png",width=6,heigh=4.8,dpi=600)

ggplot(data=str, aes(x=log2FoldChange, y=-log10(padj), col=Expression)) + 
  geom_point(alpha = 0.7) + xlim(-10,10) + ylim(0,150) +
  theme_minimal()

ggsave("str_volcano.png",width=6,heigh=4.8,dpi=600)

ggplot(data=epi, aes(x=log2FoldChange, y=-log10(padj), col=Expression)) + 
  geom_point(alpha = 0.7) + xlim(-10,10) + ylim(0,150) +
  theme_minimal()

ggsave("epi_volcano.png",width=6,heigh=4.8,dpi=600)

