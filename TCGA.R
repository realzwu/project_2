setwd("/Users/zd821/Documents/project2/TCGA")

library(openxlsx)
library(biomaRt)
library(dplyr)
library(DESeq2)
library(ggplot2)

mart <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')

### TCGA

rna.data <- read.table("TCGA.PRAD.sampleMap_HiSeqV2",sep="\t",head=T,row.names = 1)

rna.data = round(2^(rna.data)-1)
rna.data = cbind(rownames(rna.data),rna.data)

Ensembl <- rna.data[,1]

BM = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', "start_position", "end_position"), filters = 'hgnc_symbol', values = Ensembl, mart = mart)
BM =BM %>% dplyr::distinct(hgnc_symbol, .keep_all = TRUE)

eff_length <- abs(BM$end_position - BM$start_position)
BM <- cbind(BM,eff_length)

rna.data2 = merge(rna.data,BM,by.x = "rownames(rna.data)",by.y = "hgnc_symbol")
eff_list = rna.data2[,555]
rownames(rna.data2) = rna.data2[,1]
rna.data2 = rna.data2[,c(1:551)]
write.xlsx(rna.data2, "TCGA_counts.xlsx")
rna.data2 = rna.data2[,-1]

# TPM

rna.data3 <- rna.data2 / eff_list
rna.data3 = t( t(rna.data3) / colSums(rna.data3) ) * 1e6 
rna.data3 = log2(rna.data3 + 1)
rna.data3 = as.data.frame(cbind(rownames(rna.data3),rna.data3))
write.xlsx(rna.data3, "TCGA_tpm.xlsx")

### gleason score and prognosis

rna.data <- read.xlsx("TCGA_counts.xlsx",rowNames = T, colNames = F)
# rna.data <- read.xlsx("TCGA_tpm.xlsx",rowNames = T, colNames = F)
rna.data = t(rna.data)

clin.data <- read.table("TCGA.PRAD.sampleMap_PRAD_clinicalMatrix", sep="\t", header=T)
clin.data[,1]<-gsub(clin.data[,1], pattern="-", replace=".")
clin.data = clin.data[,c(1,38)]

merged_table = merge(clin.data,rna.data,by.y = "rownames(rna.data)",by.x = "sampleID")

counts <- read.table("survival%2FPRAD_survival.txt",sep="\t",head=T)

counts = counts[,c(1,9,10)]
counts[,1]<-gsub(counts[,1], pattern="-", replace=".")

counts = merge(counts,merged_table,by.y = "sampleID",by.x = "sample")

counts[grep(".01",counts[,1]),1] = "tumour"
counts[grep(".11",counts[,1]),1] = "benign"
counts[grep(".06",counts[,1]),1] = "tumour"

counts = as.data.frame(counts)
write.xlsx(counts, "TCGA_counts_new.xlsx")
# write.xlsx(counts, "TCGA_tpm_new.xlsx")


### GTEx

gtex_data<-read.table("gene_reads_2017-06-05_v8_prostate", skip = 2, header = TRUE, sep = "\t",row.names = 1)

library(stringr)
Ensembl <- gtex_data[,1]

splitEnsembl <- function(Ensembl){
  return(str_split(Ensembl[1],'[.]',simplify = T)[1])}
Ensembl <- sapply(Ensembl,splitEnsembl,simplify = T)

gtex_data[,1] <- Ensembl

BM = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', "start_position", "end_position"), filters = 'ensembl_gene_id', values = Ensembl, mart = mart)
BM =BM %>% dplyr::distinct(hgnc_symbol, .keep_all = TRUE)

eff_length <- abs(BM$end_position - BM$start_position)
BM <- cbind(BM,eff_length)

gtex_data2 = merge(gtex_data,BM,by.x = "Name",by.y = "ensembl_gene_id")
name_list = gtex_data2[,248]
eff_list = gtex_data2[,251]

write.xlsx(gtex_data2[,c(2:247)], "GTEx_counts.xlsx")
gtex_data2 = gtex_data2[,c(3:247)]

gtex_data4 <- as.matrix(gtex_data2) / eff_list
gtex_data4 = t( t(gtex_data4) / colSums(gtex_data4) ) * 1e6 
gtex_data4 = log2(gtex_data4 + 1)
gtex_data4 = as.data.frame(cbind(name_list,gtex_data4))

write.xlsx(gtex_data4, "GTEx_tpm.xlsx")


# PCA 
rna.data1 <- read.xlsx("TCGA_counts.xlsx",rowNames = F, colNames = T)
rna.data2 <- read.xlsx("GTEx_counts.xlsx",rowNames = F, colNames = T)
combinedata = merge(rna.data1,rna.data2,by.x = "rownames(rna.data)",by.y = "Description")
# combinedata = na.omit(combinedata)

combinedata = rbind(colnames(combinedata),combinedata)
combinedata[1,grep("GTEX.",combinedata[1,])] = "normal"
combinedata[1,grep(".01",combinedata[1,])] = "tumour"
combinedata[1,grep(".11",combinedata[1,])] = "benign"
combinedata[1,grep(".06",combinedata[1,])] = "tumour"
write.xlsx(combinedata, "combined.xlsx")

combinedata = combinedata[,-1]
colData0 <- data.frame(row.names = colnames(combinedata), condition = factor(combinedata[1,],levels = c("tumour","benign")))
combinedata = combinedata[-1,]
combinedata1 = apply(combinedata, 2, as.numeric)

dds0 = DESeqDataSetFromMatrix(countData = combinedata1, colData = colData0, design = ~ condition)
dds0 = DESeq(dds0)
vsd <- vst(dds0, blind=FALSE)
pcaData = plotPCA(vsd, intgroup="condition", returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=1.5, alpha=0.8) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+ stat_ellipse(type = "norm")

ggsave("pca_TCGA.png",width=5,heigh=4,dpi=600)


### boxplot3

rna.data <- read.xlsx("TCGA_boxplot3.xlsx",rowNames = F, colNames = T)

rna.data[,c(5,6)] = apply(rna.data[,c(5,6)], 2, as.numeric)
rna.data[,4] = as.character(rna.data[,4])
rna.data$DKK3 = (rna.data$DKK3 - mean(rna.data$DKK3)) / sd(rna.data$DKK3)
rna.data$Gleason_Score[rna.data$Sample == "benign"] = "Non-cancer"
rna.data$Gleason_Score = factor(rna.data$Gleason_Score,levels=c("Non-cancer", '6','7','8','9'))


library(ggprism)
library(ggbeeswarm)
library(rstatix)
library(forcats)

rna.data2 <- rna.data %>%
  rstatix::t_test(DKK3 ~ Gleason_Score, p.adjust.method = "BH", var.equal = T) %>% 
  rstatix::add_x_position(x = "Gleason_Score") 

rna.data2 = rna.data2[c(1,7),]
rna.data2$label = c("***", "p = 0.83")

# new
p <- ggplot(rna.data, aes(x = Gleason_Score, y = DKK3))
p <- p + ggbeeswarm::geom_beeswarm(aes(fill = Gleason_Score), cex = 1.0, size=1.2, alpha=0.9, shape = 21) 
p <- p + stat_summary(geom = "crossbar", aes(fill = Gleason_Score), fun = mean, position = position_dodge(0.9),
                      colour = "red", size = 0.4, width = 0.5, show.legend = FALSE)
p <- p + add_pvalue(rna.data2, y = 4.2, tip.length = 0, bracket.shorten = c(0,1),
                    fontface = "italic", lineend = "round",label.size = 4, fontface = "bold", bracket.size = 0.5)
print(p)

# old
p <- ggplot(rna.data, aes(x = Gleason_Score, y = DKK3))
p <- p + ggbeeswarm::geom_beeswarm(aes(fill = forcats::fct_rev(Sample)), dodge.width = 0.9, cex = 1, size=1.6, alpha=0.9, shape =21) 
p <- p + stat_summary(geom = "crossbar", aes(fill = forcats::fct_rev(Sample)), fun = mean, position = position_dodge(0.9),
  colour = "red", size = 0.4, width = 0.7, show.legend = FALSE)
p <- p + add_pvalue(rna.data2, y = 4.5, xmin = "xmin", xmax = "xmax", tip.length = 0, 
  fontface = "italic", lineend = "round",label.size = 4, fontface = "bold", bracket.size = 0.5)
p

ggsave("TCGA_stageplot_DKK3.png",width=7.5,heigh=4.8,dpi=600)



rna.data2 = rna.data[rna.data[,1] == "tumour",]

res_aov <- aov(DKK3 ~ Gleason_Score, data = rna.data2)
# ANOVA p-value DKK3 = 0.883
summary(res_aov)[[1]][["Pr(>F)"]][1]

### correlation

rna.data <- read.xlsx("TCGA_boxplot3.xlsx",rowNames = F, colNames = T)
rna.data = rna.data[rna.data$Sample == "tumour",]
rna.data[,5] = as.numeric(rna.data[,5])
rna.data[,6] = as.numeric(rna.data[,6])
rna.data[,4] = as.character(rna.data[,4])
ggplot(rna.data, aes(x=AR, y=DKK3)) + geom_point(alpha = 0.7) + geom_smooth(method = lm ,alpha = 0.2)


ggplot(rna.data, aes(x=AR, y=DKK3, color = Gleason_Score)) + geom_point(alpha = 0.7) + geom_smooth(method = lm ,alpha = 0.2)


### survival

library("survival")
library("party")
library("survminer")
library("gridExtra")

rna.data <- read.xlsx("TCGA_boxplot4.xlsx",rowNames = F, colNames = T)
rna.data$DKK3 = as.numeric(rna.data$DKK3)
rna.data$PFI = as.numeric(rna.data$PFI)

# scatter plot
rna.data$status[rna.data$PFI == 0] = "No Recurrence"
rna.data$status[rna.data$PFI != 0] = "Recurrence"

rna.data$DKK3 = (rna.data$DKK3 - mean(rna.data$DKK3)) / sd(rna.data$DKK3)

scatterPlot <- ggplot(rna.data,aes(PFI.year, DKK3, color=status)) + geom_point() + 
  scale_color_manual(values = c('#3DB7E4','#FF8849')) + 
  theme(legend.position=c(0.6,1), legend.justification=c(0,1))
scatterPlot
# Marginal density plot of x (top panel)
xdensity <- ggplot(rna.data, aes(PFI.year, fill=status)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#3DB7E4','#FF8849')) + 
  theme(legend.position = "none")
xdensity
# Marginal density plot of y (right panel)
ydensity <- ggplot(rna.data, aes(DKK3, fill=status)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#3DB7E4','#FF8849')) + 
  theme(legend.position = "none")
ydensity

blankPlot <- ggplot()+geom_blank(aes(1,1))+ theme(plot.background = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), panel.border = element_blank(),panel.background = element_blank(),axis.title.x = element_blank(),
  axis.title.y = element_blank(),axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank())

grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, ncol=2, nrow=2, widths=c(3.4, 2), heights=c(2, 3.4))



rna.data$geneexp_cp <- rna.data$DKK3 <= mean(rna.data$DKK3)

fit <- survfit(Surv(PFI.year, PFI) ~ geneexp_cp, data = rna.data)

ggsurvplot(fit,
           data = rna.data,
           size = 1,                 # 线条大小
           xlim = c(0,10),
           ylim = c(0,1),
           palette = c("#4575B4","#992320"),
           conf.int = TRUE, 
           pval = T,             
           censor.shape=".",
           risk.table = TRUE,        
           legend.labs = c("High DKK3", "Low DKK3"), 
           xlab='Time to Biochemical Recurrence (years)', 
           ylab='Probability of Recurrence Free',
           break.time.by = 1,
           risk.table.height = 0.25, 
           ggtheme = theme_bw())


ggsave("TCGA_survive.png",width=6,heigh=4.8,dpi=600)



### LASSO
library(glmnet)
ecm <- read.xlsx("GO0062023.xlsx",rowNames = F, colNames = F)
rna.data <- read.xlsx("TCGA_tpm_new.xlsx",rowNames = F, colNames = F)

rna.data2 = rna.data[rna.data[,1] == "tumour" & rna.data[,2] == 1,]
pfitime = as.numeric(rna.data2$X3) / 365.25
rna.data2 = rbind(rna.data[1,],rna.data2)
rna.data2 = t(rna.data2)
rna.data2 = t(merge(ecm,rna.data2,by.x = "X1", by.y = "1"))
rna.data2 = as.data.frame(rna.data2)
colnames(rna.data2) = rna.data2[1,]
rna.data2 = rna.data2[-1,]
rna.data2 = apply(rna.data2, 2, as.numeric)
rna.data4 = apply(rna.data2, 2, mean)
rna.data5 = apply(rna.data2, 2, sd)
for(i in 1:129){rna.data2[,i] = (rna.data2[,i] - rna.data4[i])/ rna.data5[i] }

rna.data3 = as.matrix(rna.data2[,c("COL1A1","FREM2","DPT","MATN3","COL10A1")])

cvfit <- cv.glmnet(rna.data3, pfitime, family = "gaussian")
plot(cvfit)
cvfit$lambda.min

fit <- glmnet(rna.data3, pfitime, family = "gaussian",intercept = T)
plot(fit, label = T)
plot(fit, xvar = "lambda", label = TRUE)
fit2 = coef(fit, s = 0.1)
fit2

## calculation of risk score
colnames(rna.data) = rna.data[1,]

rna.data6 = rna.data[-1 & rna.data$sample == "tumour",c("PFI","PFI.time","COL1A1","FREM2","DPT","MATN3","COL10A1")]
rna.data6[,c(3:7)] = apply(rna.data6[,c(3:7)], 2, as.numeric)

rna.data7 = apply(rna.data6[,c(3:7)], 2, mean)
rna.data8 = apply(rna.data6[,c(3:7)], 2, sd)
for(i in 1:5){rna.data6[,i+2] = (rna.data6[,i+2] - rna.data7[i]) / rna.data8[i] }

rna.data6$risk = -0.175 * as.numeric(rna.data6$COL1A1) - 0.129 * as.numeric(rna.data6$FREM2) +0.407 * as.numeric(rna.data6$DPT) -0.104 * as.numeric(rna.data6$MATN3) +0.431 * as.numeric(rna.data6$COL10A1)
rna.data6$PFI.year = as.numeric(rna.data6$PFI.time) / 365.25
rna.data6$PFI= as.numeric(rna.data6$PFI)
ctree_xfs <- ctree(Surv(PFI.year, PFI) ~ risk, data = rna.data6)
pvalue <- 1 - ctree_xfs@tree$criterion$maxcriterion
ps2 <- party:::cutpoints_list(ctree_xfs@tree, variableID=1)
ps  <- signif(ps2[1], digits = 3)
rna.data6$risk = rna.data6$risk <= ps
fit <- survfit(Surv(PFI.year, PFI) ~ risk, data = rna.data6)

ggsurvplot(
  fit,
  data = rna.data6,
  size = 1,                 # 线条大小
  # xlim = c(0,10),
  ylim = c(0,1),
  palette = c("#F8776D", "#02C1C5"),# 分组颜色
  conf.int = TRUE,          # 是否添加置信区间
  pval = T,              # 是否添加p值
  censor.shape=".",
  risk.table = TRUE,        # 是否添加风险表
  risk.table.col = "strata",# 分线表颜色
  legend.labs =
    c("High risk", "Low risk"),    # 图例标签
  break.time.by = 1,
  risk.table.height = 0.25, # 生存曲线图下所有生存表的高度，数值0-1之间
  ggtheme = theme_bw())

