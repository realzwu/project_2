# Main coding for MSKCC cohort

setwd("")
##Load the required libraries

library(ggplot2)
library(tidyr)
library(devtools)
library(dplyr)
library(RColorBrewer)
library(party)
library(survival)
library(prostateCancerTaylor)
iclusPal <- brewer.pal(5, 'Set1')
myfile <- "xxx.txt"
genes <- read.delim(myfile)[,1]

# genes = "DKK3"

###Convert into data convenient for dplyr

data(taylor,package = 'prostateCancerTaylor')
pd_taylor <- tbl_df(pData(taylor))
fd_taylor <- tbl_df(fData(taylor))
exp_taylor <- tbl_df(data.frame(ID = as.character(featureNames(taylor)),log2(exprs(taylor))))
probes <- fd_taylor %>% filter(Gene %in% genes) %>% dplyr::select(ID) %>% unique %>% as.matrix %>%  as.character
data <- exp_taylor  %>% filter(ID %in% probes) %>% gather(geo_accession,Expression,-ID)
fd <- fd_taylor %>% mutate(Symbol = Gene)

pd <- mutate(pd_taylor,Gleason = gsub('4+3', '7', pd_taylor$Gleason,fixed=TRUE)) %>% 
mutate(Gleason = gsub('3+4', '7', Gleason,fixed=TRUE)) %>% 
mutate(Gleason = gsub('3+3', '6', Gleason,fixed=TRUE)) %>% 
mutate(Gleason = gsub('4+5', '9', Gleason,fixed=TRUE)) %>% 
mutate(Gleason = gsub('4+4', '8', Gleason,fixed=TRUE)) %>% 
mutate(Gleason = gsub('5+3', '8', Gleason,fixed=TRUE)) %>% 
mutate(Gleason=factor(Gleason,levels=c("Non-cancer", '6','7','8','9')))

summary_stats <- data %>% group_by(ID) %>% 
summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))

data <- left_join(data,summary_stats) %>% mutate(Z = (Expression - mean) / sd)

mostVarProbes <- left_join(summary_stats,fd) %>% 
arrange(Symbol,desc(iqr)) %>% 
distinct(Symbol,.keep_all=TRUE) %>% 
dplyr::select(ID) %>%  as.matrix %>%  as.character

data <- filter(data, ID %in% mostVarProbes)
data <- left_join(data, dplyr::select(fd, ID, Symbol))
data <- left_join(data, pd)
data <- mutate(data, Expression=Z)
data$Gleason[data$Sample_Group == "normal adjacent benign prostate"] = "Non-cancer"
data = na.omit(data[,c(1:11)])

## Jitter plot

library(ggprism)
library(ggbeeswarm)
library(rstatix)
library(forcats)

data$Symbol=factor(data$Symbol,levels=c("Non-cancer", 'FBLN1','ANGPT1','DPT','FREM2'))

rna.data2 <- data %>%
  group_by(Symbol) %>%
  rstatix::t_test(Z ~ Gleason, p.adjust.method = "BH", 
                  var.equal = T) %>% 
  rstatix::add_x_position(x = "Gleason")

rna.data2 = rna.data2[c(1,7,11,17,21,27,31,37),]
rna.data2$label =  c(" ", "***","***", "***"," ", "**","***", "**")

p <- ggplot(data, aes(x = Gleason, y = Expression))
p <- p + ggbeeswarm::geom_beeswarm(aes(fill = Gleason), cex = 2.0, size=1.2, alpha=0.9, shape = 21) 
p <- p + facet_wrap(  ~ Symbol, scales = "free", labeller = labeller(Symbol = c(FBLN1 = "FBLN1", ANGPT1 = "ANGPT1", DPT = "DPT", FREM2 = "FREM2")))

p <- p + stat_summary(geom = "crossbar", aes(fill = Gleason), fun = mean, position = position_dodge(0.9),
                      colour = "red", size = 0.4, width = 0.5, show.legend = FALSE)
p <- p + add_pvalue(rna.data2, y = c(3,3,4,4,4,4,3,3), tip.length = 0, bracket.shorten = c(0,1,0,1,0,1,0,1),
                    fontface = "italic", lineend = "round",label.size = 4, fontface = "bold", bracket.size = 0.5)
print(p)

ggsave("MSKCC_stage.png",width=6,heigh=4.8,dpi=600)


## K_M plot

data <- left_join(data,summary_stats) %>% mutate(Z = (Expression - mean) / sd)

mostVarProbes <- left_join(summary_stats,fd) %>% 
  arrange(Symbol,desc(iqr)) %>% 
  distinct(Symbol,.keep_all=TRUE) %>% 
  dplyr::select(ID) %>%  as.matrix %>%  as.character
data <- filter(data, ID %in% mostVarProbes)
data <- left_join(data, dplyr::select(fd, ID, Symbol))
data <- left_join(data, pd)
data <- data %>% filter(!is.na(Time) & !is.na(Event))
surv.xfs <- Surv((as.numeric(as.character(data$Time))/12), data$Event)
data$surv.xfs <- surv.xfs

ctree_xfs <- ctree(surv.xfs ~ Expression, data = data)
pvalue <- 1 - ctree_xfs@tree$criterion$maxcriterion
newPval <- signif(pvalue, digits = 2)

if(newPval < 0.05) {
  ps2 <- party:::cutpoints_list(ctree_xfs@tree, variableID=1)
  ps  <- signif(ps2[1], digits = 3)
  if(length(ps2)==1) {
    data$geneexp_cp <- data$Expression <= ps2[1]
    nt <- table(data$geneexp_cp)
    fit <- survfit(surv.xfs ~ geneexp_cp, data = data)
    plot(fit, xlab='Time to BCR (years)', ylab='Probability of Freedom from Biochemical Recurrence', main=paste(genes,', p=', newPval), col=c(2,4))
    legend('bottomleft', c(paste(genes, '>', ps, 'n=', nt[[1]]), paste(genes, '<=', ps, 'n=', nt[[2]])), col=c(2,4), lty=1, lwd=1.5, bty='n')
  }
} else {
  ps <- round(median(data$Expression),3)
  data$geneexp_cp <- data$Expression<= ps
  nt <- table(data$geneexp_cp)
  fit <- survfit( surv.xfs ~ geneexp_cp, data = data)
  test <- survdiff(surv.xfs~geneexp_cp,data=data)
  newPval <- round(pchisq(test$chisq, df = length(test$n)-1, lower.tail=FALSE),3)
  plot(fit, xlab='Time to BCR (years)', ylab='Probability of Recurrence Free', main=paste(genes,', p=', newPval), col=c(2,4))
  legend('bottomleft', c(paste(genes, '>', ps, 'n=', nt[[1]]), paste(genes, '<=', ps, 'n=', nt[[2]])), col=c(2,4), lty=1, lwd=1.5, bty='n')
}

ggsurvplot(fit,
  data = data,
  size = 1,                 # 线条大小
  xlim = c(0,10),
  ylim = c(0,1),
  palette = c("#F8776D", "#02C1C5"),
  conf.int = TRUE,          # 是否添加置信区间
  pval = T,             
  censor.shape=".",
  risk.table = TRUE,        # 是否添加风险表
  risk.table.col = "strata",# 分线表颜色
  legend.labs = c("High DKK3", "Low DKK3"),    # 图例标签
  xlab='Time to Biochemical Recurrence (years)', 
  ylab='Probability of Recurrence Free',
  break.time.by = 1,
  risk.table.height = 0.25, 
  ggtheme = theme_bw())

ggsave("MSKCC_survive.png",width=6,heigh=4.8,dpi=600)


### GGally for making correlation matrix
library(GGally)
library(forcats)

exp_taylor2 = exp_taylor %>% filter(ID %in% c("13551","6553","5779","6495","25055"))
# ID sequence:74 463 6495 6911 25055
exp_taylor2 = t(exp_taylor2[,-1])
exp_taylor4 = apply(exp_taylor2, 2, mean)
exp_taylor5 = apply(exp_taylor2, 2, sd)
for(i in 1:5){exp_taylor2[,i] = (exp_taylor2[,i] - exp_taylor4[i])/ exp_taylor5[i] }

exp_taylor2 = cbind(rownames(exp_taylor2),exp_taylor2)

data2 = unique(data[,c(2,16,17)])

exp_taylor3 = merge(exp_taylor2,data,by.x = "V1", by.y ="geo_accession")

exp_taylor4 = cbind(exp_taylor3$Gleason,exp_taylor3[,c(2:6)])
colnames(exp_taylor4) = c("Gleason","DKK3","FBLN1","ANGPT1","DPT","FREM2" )

exp_taylor4[,1] = as.character(exp_taylor4[,1])
exp_taylor4[,c(2:6)] = apply(exp_taylor4[,c(2:6)],2,as.numeric)
exp_taylor4 = exp_taylor4 %>% unique

exp_taylor4[,1][exp_taylor4[,1] == "6"] = "≤ 7"
exp_taylor4[,1][exp_taylor4[,1] == "7"] = "≤ 7"
exp_taylor4[,1][exp_taylor4[,1] == "8"] = "≥ 8"
exp_taylor4[,1][exp_taylor4[,1] == "9"] = "≥ 8"

ggpairs(exp_taylor4,columns = 2:6, ggplot2::aes(colour = forcats::fct_rev(Gleason)))
