setwd("/Users/zd821/Documents/project2/mRNA/weihuanjing")

library(openxlsx)
library(biomaRt)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(openxlsx)
library(ggprism)
library(ggbeeswarm)
library(rstatix)
library(forcats)
library(GGally)

# jitter plot

rna.data <- read.xlsx("microenvironment_zscore.xlsx",rowNames = F, colNames = T)

# rna.DKK3 = rna.data[rna.data$namelist == "DKK3",]
# rna.DKK3 = t(rbind(rna.data[c(1,3),-1],rna.DKK3[2,-1],rna.DKK3[2,-1]))
# colnames(rna.DKK3) = c("Gleason", "Sample", "DKK3")

# rna.DKK3 = as.data.frame(rna.DKK3)
# rna.DKK3[,1] = as.character(rna.DKK3[,1])
# rna.DKK3[,3] = as.numeric(rna.DKK3[,3])

rna.DKK3 = rna.data[rna.data$namelist == "DKK3",][2,]
rna.FBLN1 = rna.data[rna.data$namelist == "FBLN1",]
rna.ANGPT1 = rna.data[rna.data$namelist == "ANGPT1",]
rna.DPT = rna.data[rna.data$namelist == "DPT",]
rna.FREM2 = rna.data[rna.data$namelist == "FREM2",]

rna.all = t(rbind(rna.data[c(1,3),-1],rna.DKK3[1,-1],rna.FBLN1[1,-1],rna.ANGPT1[1,-1],rna.DPT[1,-1],rna.FREM2[1,-1]))
colnames(rna.all) = c("Gleason", "Sample", "DKK3", "FBLN1", "ANGPT1", "DPT", "FREM2")

rna.all = as.data.frame(rna.all)
rna.all[,1] = as.character(rna.all[,1])
rna.all[,c(3:7)] = apply(rna.all[,c(3:7)],2,as.numeric)
rna.all$Gleason[rna.all$Gleason == "6"] = "= 6"
rna.all$Gleason[rna.all$Gleason == "≥8"] = "≥ 8"
                
rna.sT = rna.all[rna.all$Sample == "sT",]
rna.T = rna.all[rna.all$Sample == "T",]
rna.B = rna.all[rna.all$Sample == "B",]

ggpairs(rna.all[,c(1,3:7)],columns = 2:6, ggplot2::aes(colour=forcats::fct_rev(Gleason)))
ggpairs(rna.sT[,c(1,3:7)],columns = 2:6, ggplot2::aes(colour=forcats::fct_rev(Gleason)))
ggpairs(rna.B[,c(1,3:7)],columns = 2:6, ggplot2::aes(colour=forcats::fct_rev(Gleason)))
ggpairs(rna.T[,c(1,3:7)],columns = 2:6, ggplot2::aes(colour= forcats::fct_rev(Gleason)))


rna.all$Sample = factor(rna.all$Sample,levels=c("B", 'P','T','sB','sP',"sT"))

# t test
rna.label <- rna.all %>% group_by(Sample) %>%
  rstatix::t_test(DKK3 ~ Gleason, p.adjust.method = "BH", 
                  var.equal = T) %>% 
  rstatix::add_x_position(x = "Sample", dodge = 0.9)
rna.label$label = c("p = 0.25", "p = 0.88","p = 0.35","p = 0.09", "p = 0.82", "*")

rna.label2 <- rna.all %>%
  rstatix::t_test(DKK3 ~ Sample, p.adjust.method = "BH", 
                  var.equal = T) %>% 
  rstatix::add_x_position(x = "Sample") %>% rstatix::add_significance(p.col = "p")

# label t-test
rna.label2 = rna.label2[c(2,14),]
rna.label2$label = c("**", "*")


p <- ggplot(rna.all, aes(x = Sample, y = DKK3))
p <- p + ggbeeswarm::geom_beeswarm(aes(fill = forcats::fct_rev(Gleason)), dodge.width = 0.9, cex = 1, size=2, alpha=0.9, shape =21) 
p <- p + stat_summary(geom = "crossbar", aes(fill = forcats::fct_rev(Gleason)), fun = mean, position = position_dodge(0.9),
                      colour = "red", size = 0.4, width = 0.7, show.legend = FALSE)
p <- p + add_pvalue(rna.label, xmin = "xmin",xmax = "xmax",y = 3, tip.length = 0,
                    fontface = "italic", lineend = "round",label.size = 4, fontface = "bold", bracket.size = 0.5)                   
p <- p + add_pvalue(rna.label2, y = 3.5,xmin = "xmin",xmax = "xmax", tip.length = 0,
                    fontface = "italic", lineend = "round",label.size = 4, fontface = "bold", bracket.size = 0.5)
p

ggsave("weihuanjing_stageplot_DKK3.png",width=8,heigh=4.8,dpi=600)

