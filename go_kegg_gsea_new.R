## install.packages("devtools")
# devtools::install_github("YuLab-SMU/enrichplot")
# sink("output.txt", append = TRUE, type = c("output", "message"),split = FALSE)

setwd('')

library(openxlsx)
library(ggplot2)
library(stringr)
require(enrichplot)
library(clusterProfiler)
library(DOSE)
library(ggnewscale)
library(topGO)
library(ComplexHeatmap)
library(ggridges)
library(pathview)


# Read DEGs from edgeR

info <- read.xlsx( "./Stromal/DEG-all.gene_14573.xlsx", rowNames = F,colNames = T)
GO_database <- 'org.Hs.eg.db'
KEGG_database <- 'hsa'

gene <- bitr(info$gene_id,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
gene2 = merge(gene,info,by.y = "gene_id", by.x = 'SYMBOL')
gene2 = gene2[,-4]
gene3 <- gene2[which(gene2$PValue < 0.05 & abs(gene2$log2FoldChange) > 1),]

# Over-representative analysis GO
GO<-enrichGO( gene3$ENTREZID, OrgDb = GO_database, keyType = "ENTREZID", ont = "ALL", pvalueCutoff = 0.01, readable = T)
dotplot(GO,showCategory=10, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")


# GSEA GO
names(gene2) <- c('SYMBOL','ENTREZID','Log2FoldChange','pvalue','FDR')
GSEA_input <- gene2$Log2FoldChange
names(GSEA_input) = gene2$ENTREZID
GSEA_input = sort(GSEA_input, decreasing = TRUE)

GSEA_GO2 <- gseGO(GSEA_input, OrgDb = org.Hs.eg.db, ont = "CC", pvalueCutoff = 0.2)

# GSEA plot with Gene annotation

which(GSEA_GO2$Description == "collagen-containing extracellular matrix")
GSEA_GO2[11,]

GSEA_GO3 = setReadable(GSEA_GO2,GO_database,"ENTREZID") 
id = GSEA_GO3$ID[1]
GSEA_GO3[[id]]
g2 = sample(GSEA_GO3[[id]],57)
g2[c(1:47)] = "."
p1 = gseaplot2(GSEA_GO3, geneSetID = 11, title = GSEA_GO3$Description[11])
p1[[1]] = p1[[1]] + geom_gsea_gene(g2, geom=ggrepel::geom_text_repel)
p1

# ridgeplot(GSEA_GO2,10) 
# goplot(GSEA_GO2,showCategory = 20)



# KEGG GSEA

GSEA_KEGG <- gseKEGG(GSEA_input, organism = KEGG_database, pvalueCutoff = 0.01)#GSEA富集分析

head(GSEA_KEGG,20)
ridgeplot(GSEA_KEGG) 
gseaplot2(GSEA_KEGG,1)
gseaplot2(GSEA_KEGG,2)
gseaplot2(GSEA_KEGG,3)
gseaplot2(GSEA_KEGG,4)
gseaplot2(GSEA_KEGG,6)
gseaplot2(GSEA_KEGG,9)

hsa03010 <- pathview(gene.data  = GSEA_input,
                     pathway.id = "hsa03010",
                     species    = "hsa",
                     limit      = list(gene=max(abs(GSEA_input)), cpd=1))



