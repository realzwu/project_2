## DEG analysis

setwd("")
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("edgeR")
# BiocManager::install("DESeq2")

library(openxlsx)
library(edgeR)

counts <- read.xlsx("featureCounts_sorted.xlsx", colNames = T, rowNames = T)

group <- c(rep(1,3),rep(2,3))

y = DGEList(counts=counts[c(16:18,22:24)], group=group)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

y <- estimateDisp(y)

et <- exactTest(y)

et <- topTags(et, n=100000)

et <- as.data.frame(et)

et <- cbind(rownames(et),et)

colnames(et) <- c("gene_id", "log2FoldChange", "log2CPM", "PValue", "FDR")

write.xlsx(et, "DEG-all.gene.xlsx")

etSig <- et[which(et$PValue < 0.01 & abs(et$log2FoldChange) > 1),]

etSig[which(etSig$log2FoldChange > 0), "up_down"] <- "Up"
etSig[which(etSig$log2FoldChange < 0), "up_down"] <- "Down"

write.xlsx(etSig, "Pval0.01-FC2-edgeR.xlsx")

