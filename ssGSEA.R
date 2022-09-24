setwd("")

library(openxlsx)
library(ggpubr)
library(corrplot)
library(Hmisc)
library(GSVA)
library(reshape2)

cleaned <- read.xlsx("microenvironment_zscore.xlsx",sheet = 1,rowNames = F,colNames = T)

# Unique DKK3

which(cleaned$namelist == "DKK3")
cleaned = cleaned[-2819,]
which(cleaned$namelist == "DKK3")

# Read ECM-related gene list

myfile <- "ADAMTS10.txt"
corrlist <- read.delim(myfile)[,1]

### calculating ssGSEA

gsea = cleaned[!duplicated(cleaned[,1]),]
roi = gsea[3,]
gsea = gsea[c(-2,-3),]
rownames(gsea) = gsea[,1]
namelist = gsea[-1,1]
gsea_T = gsea[,roi == "T"]
gsea_P = gsea[,roi == "P"]
gsea_B = gsea[,roi == "B"]
gsea_sB = gsea[,roi == "sB"]
gsea_sP = gsea[,roi == "sP"]
gsea_sT = gsea[,roi == "sT"]


gleasonlist = gsea_T[1,]

gsea_B = apply(gsea_B[-1,],2,as.numeric)
gsea_P = apply(gsea_P[-1,],2,as.numeric)
gsea_T = apply(gsea_T[-1,],2,as.numeric)

gsea_sB = apply(gsea_sB[-1,],2,as.numeric)
gsea_sP = apply(gsea_sP[-1,],2,as.numeric)
gsea_sT = apply(gsea_sT[-1,],2,as.numeric)


rownames(gsea_B) = namelist
rownames(gsea_P) = namelist
rownames(gsea_T) = namelist
rownames(gsea_sB) = namelist
rownames(gsea_sP) = namelist
rownames(gsea_sT) = namelist


dkkB = gsea_B[rownames(gsea_B) == "DKK3",]
dkkP = gsea_P[rownames(gsea_P) == "DKK3",]
dkkT = gsea_T[rownames(gsea_T) == "DKK3",]
dkksB = gsea_sB[rownames(gsea_sB) == "DKK3",]
dkksP = gsea_sP[rownames(gsea_sP) == "DKK3",]
dkksT = gsea_sT[rownames(gsea_sT) == "DKK3",]

dkkB[dkkB < median(dkkB)] = "Low DKK3"
dkkB[dkkB != "Low DKK3"] = "High DKK3"
dkkP[dkkP < median(dkkP)] = "Low DKK3"
dkkP[dkkP != "Low DKK3"] = "High DKK3"
dkkT[dkkT < median(dkkT)] = "Low DKK3"
dkkT[dkkT != "Low DKK3"] = "High DKK3"

dkksB[dkksB < median(dkksB)] = "Low DKK3"
dkksB[dkksB != "Low DKK3"] = "High DKK3"
dkksP[dkksP < median(dkksP)] = "Low DKK3"
dkksP[dkksP != "Low DKK3"] = "High DKK3"
dkksT[dkksT < median(dkksT)] = "Low DKK3"
dkksT[dkksT != "Low DKK3"] = "High DKK3"

DKK3 = rbind(as.matrix(dkkB),as.matrix(dkkP),as.matrix(dkkT),as.matrix(dkksB),as.matrix(dkksP),as.matrix(dkksT))

gs <- as.list(sample(0, size=1, replace=TRUE))
gs[[1]] = corrlist

ssgsea_score3_B = gsva(as.matrix(gsea_B), as.list(gs), method = "ssgsea", ssgsea.norm = T, verbose = T)   # signature 'matrix,list'
ssgsea_score3_P = gsva(as.matrix(gsea_P), as.list(gs), method = "ssgsea", ssgsea.norm = T, verbose = T)   
ssgsea_score3_T = gsva(as.matrix(gsea_T), as.list(gs), method = "ssgsea", ssgsea.norm = T, verbose = T)   
ssgsea_score3_sB = gsva(as.matrix(gsea_sB), as.list(gs), method = "ssgsea", ssgsea.norm = T, verbose = T)  
ssgsea_score3_sP = gsva(as.matrix(gsea_sP), as.list(gs), method = "ssgsea", ssgsea.norm = T, verbose = T)  
ssgsea_score3_sT = gsva(as.matrix(gsea_sT), as.list(gs), method = "ssgsea", ssgsea.norm = T, verbose = T)  

ssgsea_score = rbind(ssgsea_score3_B,ssgsea_score3_P,ssgsea_score3_T,ssgsea_score3_sB,ssgsea_score3_sP,ssgsea_score3_sT)
rownames(ssgsea_score) = c("B","P","T","sB","sP","sT")


ssgsea_score2 = melt(t(ssgsea_score))
ssgsea_score2[,1] = DKK3
colnames(ssgsea_score2) = c("DKK3_level","Sample","ssGSEA")

library(ggprism)
library(ggbeeswarm)
library(rstatix)
library(forcats)
library(reshape)

# Draw jitter plot for ssGSEA Enrichent score

rna.data2 <- ssgsea_score2 %>%  group_by(Sample) %>%
  rstatix::t_test(ssGSEA ~ DKK3_level, p.adjust.method = "BH", 
                  var.equal = T, ref.group = "Low DKK3") %>% 
  rstatix::add_x_position(x = "Sample", dodge = 0.7) %>% 
  mutate(label = c("**", "*", "**","***", "**", "***"))

rna.label2 <- ssgsea_score2 %>%
  rstatix::t_test(ssGSEA ~ Sample, p.adjust.method = "BH", 
                  var.equal = T) %>% 
  rstatix::add_x_position(x = "Sample")

rna.label2 = rna.label2[c(2,14),]
rna.label2$label = c("**", "***")

p <- ggplot(ssgsea_score2, aes(x = Sample, y = ssGSEA))
p <- p + ggbeeswarm::geom_beeswarm(aes(fill = DKK3_level), dodge.width = 0.7, cex = 1, size=2.2, alpha=0.8, shape = 21) 
p <- p + stat_summary(geom = "crossbar", aes(fill = DKK3_level), fun = mean, position = position_dodge(0.9),
                      colour = "red", size = 0.4, width = 0.4, show.legend = FALSE)
p <- p + add_pvalue(rna.data2, y = 1.5, xmin = "xmin", xmax = "xmax", tip.length = 0, 
                    fontface = "italic", lineend = "round",label.size = 4.2, fontface = "bold", bracket.size = 0.5)
p <- p + add_pvalue(rna.label2, y = 1.7, xmin = "xmin", xmax = "xmax", tip.length = 0, 
                    fontface = "italic", lineend = "round",label.size = 4.2, fontface = "bold", bracket.size = 0.5)
p <- p + theme_prism(base_fontface = "plain", base_line_size = 0.7, base_family = "Arial")
p <- p + scale_x_discrete(guide = guide_prism_bracket(width = 0.07), labels = scales::wrap_format(5)) 
p <- p + scale_y_continuous(limits = c(-1, 2), expand = c(0, 0), breaks = seq(-1, 1.5, 0.5), guide = "prism_offset") + labs(y = "ssGSEA Enrichment Score")

p

# ggsave("TCGA_stageplot_DKK3.png",width=7.5,heigh=4.8,dpi=600)

