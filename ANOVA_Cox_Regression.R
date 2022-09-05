## For Figure 6A 6B

setwd("")

### sort four favourable gene sets
library(survival)
library(openxlsx)
myfile <- "FBLN1.txt"
genes <- read.delim(myfile)[,1]
## genes = "DKK3" "FBLN1" "ANGPT1" "DPT" "FREM2"

rna.data_ <- read.xlsx("TCGA_tpm.xlsx",rowNames = T, colNames = F)
clinical = rna.data_[,c(1:3)]
rna.data3 = t(rna.data_) 
rna.data3 = rna.data3 %>% as.data.frame %>% filter(sample %in% genes)
rna.data3 = cbind(clinical,t(rna.data3))
rna.data3 = as.data.frame(cbind(rna.data3[,3],rna.data3[,5],rna.data3[,7],rna.data3[,4],rna.data3[,6],rna.data3[,8]))
colnames(rna.data3) = rna.data3[1,]
rna.data3 = rna.data3[-1,]

### draw correlation matrix
library(GGally)
library(forcats)
rna.data3[,1] = as.character(rna.data3[,1])
rna.data3[,c(2:6)] = apply(rna.data3[,c(2:6)],2,as.numeric)
rna.data4 = rna.data3
rna.data4$gleason_score[rna.data4$gleason_score == "6"] = "≤ 7"
rna.data4$gleason_score[rna.data4$gleason_score == "7"] = "≤ 7"
rna.data4$gleason_score[rna.data4$gleason_score == "8"] = "≥ 8"
rna.data4$gleason_score[rna.data4$gleason_score == "9"] = "≥ 8"
rna.data4$gleason_score[rna.data4$gleason_score == "10"] = "≥ 8"

ggpairs(rna.data4,columns = 2:6, ggplot2::aes(colour= forcats::fct_rev(gleason_score)))


## ANOVA

df = as.data.frame(matrix(nrow=1,ncol=0))

res_aov <- aov(ADAM19 ~ gleason_score, data = rna.data3)
df$ADAM19 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov(ADAMTS15 ~ gleason_score, data = rna.data3)
df$ADAMTS15 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov(AEBP1 ~ gleason_score, data = rna.data3)
df$AEBP1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov(AHSG ~ gleason_score, data = rna.data3)
df$AHSG = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov(ANG ~ gleason_score, data = rna.data3)
df$ANG = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov(ANGPT1 ~ gleason_score, data = rna.data3)
df$ANGPT1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov(ANGPTL2 ~ gleason_score, data = rna.data3)
df$ANGPTL2 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov(ANGPT4 ~ gleason_score, data = rna.data3)
df$ANGPT4 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov(ANXA2P2 ~ gleason_score, data = rna.data3)
df$ANXA2P2 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov(ANXA6 ~ gleason_score, data = rna.data3)
df$ANXA6 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov(ANXA8 ~ gleason_score, data = rna.data3)
df$ANXA8 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov(APLP1 ~ gleason_score, data = rna.data3)
df$APLP1 = summary(res_aov)[[1]][["Pr(>F)"]][1]

 
res_aov <- aov(BMP7 ~ gleason_score, data = rna.data3)
df$BMP7 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov(CCDC80 ~ gleason_score, data = rna.data3)
df$CCDC80 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov(CDH13 ~ gleason_score, data = rna.data3)
df$CDH13 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov(CDH2 ~ gleason_score, data = rna.data3)
df$CDH2 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov(COCH ~ gleason_score, data = rna.data3)
df$COCH = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov(COL10A1 ~ gleason_score, data = rna.data3)
df$COL10A1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov(COL12A1 ~ gleason_score, data = rna.data3)
df$COL12A1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov(COL14A1 ~ gleason_score, data = rna.data3)
df$COL14A1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov(COL16A1 ~ gleason_score, data = rna.data3)
df$COL16A1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov(COL1A1 ~ gleason_score, data = rna.data3)
df$COL1A1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov(COL23A1 ~ gleason_score, data = rna.data3)
df$COL23A1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov(COL25A1 ~ gleason_score, data = rna.data3)
df$COL25A1 = summary(res_aov)[[1]][["Pr(>F)"]][1]

res_aov <- aov( COL2A1 ~ gleason_score, data = rna.data3)
df$COL2A1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( COL4A1 ~ gleason_score, data = rna.data3)
df$COL4A1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( COL4A2 ~ gleason_score, data = rna.data3)
df$COL4A2 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( COL4A4 ~ gleason_score, data = rna.data3)
df$COL4A4 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( COL4A5 ~ gleason_score, data = rna.data3)
df$COL4A5 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( COL4A6 ~ gleason_score, data = rna.data3)
df$COL4A6 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( COL5A1 ~ gleason_score, data = rna.data3)
df$COL5A1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( COL6A2 ~ gleason_score, data = rna.data3)
df$COL6A2 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( COL6A3 ~ gleason_score, data = rna.data3)
df$COL6A3 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( COL7A1 ~ gleason_score, data = rna.data3)
df$COL7A1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( COL8A1 ~ gleason_score, data = rna.data3)
df$COL8A1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( COL9A2 ~ gleason_score, data = rna.data3)
df$COL9A2 = summary(res_aov)[[1]][["Pr(>F)"]][1]


res_aov <- aov( CPA3 ~ gleason_score, data = rna.data3)
df$CPA3 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( CSPG4 ~ gleason_score, data = rna.data3)
df$CSPG4 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( CTSB ~ gleason_score, data = rna.data3)
df$CTSB = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( DCN ~ gleason_score, data = rna.data3)
df$DCN = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( DPT ~ gleason_score, data = rna.data3)
df$DPT = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( DST ~ gleason_score, data = rna.data3)
df$DST = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( EDIL3 ~ gleason_score, data = rna.data3)
df$EDIL3 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( EFEMP1 ~ gleason_score, data = rna.data3)
df$EFEMP1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( EFEMP2 ~ gleason_score, data = rna.data3)
df$EFEMP2 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( EMILIN1 ~ gleason_score, data = rna.data3)
df$EMILIN1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( ENTPD2 ~ gleason_score, data = rna.data3)
df$ENTPD2 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( FBLN1 ~ gleason_score, data = rna.data3)
df$FBLN1 = summary(res_aov)[[1]][["Pr(>F)"]][1]


res_aov <- aov( FBLN2 ~ gleason_score, data = rna.data3)
df$FBLN2 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( FBN1 ~ gleason_score, data = rna.data3)
df$FBN1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( FBN2 ~ gleason_score, data = rna.data3)
df$FBN2 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( FGFBP3 ~ gleason_score, data = rna.data3)
df$FGFBP3 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( FGFR2 ~ gleason_score, data = rna.data3)
df$FGFR2 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( FLG ~ gleason_score, data = rna.data3)
df$FLG = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( FN1 ~ gleason_score, data = rna.data3)
df$FN1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( FRAS1 ~ gleason_score, data = rna.data3)
df$FRAS1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( FREM2 ~ gleason_score, data = rna.data3)
df$FREM2 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( GDF15 ~ gleason_score, data = rna.data3)
df$GDF15 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( GPC5 ~ gleason_score, data = rna.data3)
df$GPC5 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( GPC6 ~ gleason_score, data = rna.data3)
df$GPC6 = summary(res_aov)[[1]][["Pr(>F)"]][1]


res_aov <- aov( HMCN1 ~ gleason_score, data = rna.data3)
df$HMCN1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( KRT1 ~ gleason_score, data = rna.data3)
df$KRT1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( L1CAM ~ gleason_score, data = rna.data3)
df$L1CAM = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( LAMA1 ~ gleason_score, data = rna.data3)
df$LAMA1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( LAMA2 ~ gleason_score, data = rna.data3)
df$LAMA2 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( LAMA3 ~ gleason_score, data = rna.data3)
df$LAMA3 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( LAMB1 ~ gleason_score, data = rna.data3)
df$LAMB1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( LAMB2 ~ gleason_score, data = rna.data3)
df$LAMB2 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( LAMB4 ~ gleason_score, data = rna.data3)
df$LAMB4 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( LAMC3 ~ gleason_score, data = rna.data3)
df$LAMC3 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( LGALS1 ~ gleason_score, data = rna.data3)
df$LGALS1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( LOXL2 ~ gleason_score, data = rna.data3)
df$LOXL2 = summary(res_aov)[[1]][["Pr(>F)"]][1]


res_aov <- aov( LRRC15 ~ gleason_score, data = rna.data3)
df$LRRC15 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( LTBP3 ~ gleason_score, data = rna.data3)
df$LTBP3 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( LUM ~ gleason_score, data = rna.data3)
df$LUM = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( MATN3 ~ gleason_score, data = rna.data3)
df$MATN3 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( MDK ~ gleason_score, data = rna.data3)
df$MDK = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( MFAP2 ~ gleason_score, data = rna.data3)
df$MFAP2 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( MFAP5 ~ gleason_score, data = rna.data3)
df$MFAP5 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( MGP ~ gleason_score, data = rna.data3)
df$MGP = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( MMP2 ~ gleason_score, data = rna.data3)
df$MMP2 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( MMP9 ~ gleason_score, data = rna.data3)
df$MMP9 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( MXRA5 ~ gleason_score, data = rna.data3)
df$v = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( NAV2 ~ gleason_score, data = rna.data3)
df$NAV2 = summary(res_aov)[[1]][["Pr(>F)"]][1]

res_aov <- aov( NCAM1 ~ gleason_score, data = rna.data3)
df$NCAM1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( NID1 ~ gleason_score, data = rna.data3)
df$NID1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( NID2 ~ gleason_score, data = rna.data3)
df$NID2 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( NPNT ~ gleason_score, data = rna.data3)
df$NPNT = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( NTN1 ~ gleason_score, data = rna.data3)
df$NTN1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( NTN4 ~ gleason_score, data = rna.data3)
df$NTN4 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( NTNG1 ~ gleason_score, data = rna.data3)
df$NTNG1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( NTNG2 ~ gleason_score, data = rna.data3)
df$NTNG2 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( PCOLCE ~ gleason_score, data = rna.data3)
df$PCOLCE = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( PDGFB ~ gleason_score, data = rna.data3)
df$PDGFB = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( POSTN ~ gleason_score, data = rna.data3)
df$POSTN = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( RTBDN ~ gleason_score, data = rna.data3)
df$RTBDN = summary(res_aov)[[1]][["Pr(>F)"]][1]


res_aov <- aov( S100A4 ~ gleason_score, data = rna.data3)
df$S100A4 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( S100A6 ~ gleason_score, data = rna.data3)
df$S100A6 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( S100A7 ~ gleason_score, data = rna.data3)
df$S100A7 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( S100A8 ~ gleason_score, data = rna.data3)
df$S100A8 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( S100A9 ~ gleason_score, data = rna.data3)
df$S100A9 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( SDC2 ~ gleason_score, data = rna.data3)
df$SDC2 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( SEMA3B ~ gleason_score, data = rna.data3)
df$SEMA3B = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( SEMA7A ~ gleason_score, data = rna.data3)
df$SEMA7A = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( SERPINB12 ~ gleason_score, data = rna.data3)
df$SERPINB12 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( SERPINE1 ~ gleason_score, data = rna.data3)
df$SERPINE1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( SERPINE2 ~ gleason_score, data = rna.data3)
df$SERPINE2 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( SERPINF1 ~ gleason_score, data = rna.data3)
df$SERPINF1 = summary(res_aov)[[1]][["Pr(>F)"]][1]


res_aov <- aov( SERPING1 ~ gleason_score, data = rna.data3)
df$SERPING1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( SFRP1 ~ gleason_score, data = rna.data3)
df$SFRP1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( SHH ~ gleason_score, data = rna.data3)
df$SHH = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( SMOC1 ~ gleason_score, data = rna.data3)
df$SMOC1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( SMOC2 ~ gleason_score, data = rna.data3)
df$SMOC2 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( SOD3 ~ gleason_score, data = rna.data3)
df$SOD3 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( SPARC ~ gleason_score, data = rna.data3)
df$SPARC = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( SRPX ~ gleason_score, data = rna.data3)
df$SRPX = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( SRPX2 ~ gleason_score, data = rna.data3)
df$SRPX2 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( SSC5D ~ gleason_score, data = rna.data3)
df$SSC5D = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( TGFB1I1 ~ gleason_score, data = rna.data3)
df$TGFB1I1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( TGFB2 ~ gleason_score, data = rna.data3)
df$TGFB2 = summary(res_aov)[[1]][["Pr(>F)"]][1]


res_aov <- aov( TGM2 ~ gleason_score, data = rna.data3)
df$TGM2 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( TGM4 ~ gleason_score, data = rna.data3)
df$TGM4 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( THBS1 ~ gleason_score, data = rna.data3)
df$THBS1 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( THBS2 ~ gleason_score, data = rna.data3)
df$THBS2 = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( TNC ~ gleason_score, data = rna.data3)
df$TNC = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( VCAN ~ gleason_score, data = rna.data3)
df$VCAN = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( VIT ~ gleason_score, data = rna.data3)
df$VIT = summary(res_aov)[[1]][["Pr(>F)"]][1]
res_aov <- aov( VWF ~ gleason_score, data = rna.data3)
df$VWF = summary(res_aov)[[1]][["Pr(>F)"]][1]

df= as.data.frame(t(df))
df = cbind(rownames(df),df)
write.xlsx(df,"ANOVA_TCGA.xlsx",)

rna.data3 = apply(rna.data3, 2, as.numeric)
rna.data3 = as.data.frame(rna.data3)


## Cox Regression

df = as.data.frame(matrix(nrow=4,ncol=0))

cox_result = coxph(Surv(PFI.time, PFI) ~ ADAM19, rna.data3)
                                      df$ADAM19 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ ADAMTS15, rna.data3)
                                      df$ADAMTS15 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ AEBP1, rna.data3)
                                      df$AEBP1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ AHSG, rna.data3)
                                      df$AHSG = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ ANG, rna.data3)
                                      df$ANG = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ ANGPT1, rna.data3)
                                      df$ANGPT1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ ANGPT4, rna.data3)
                                      df$ANGPT4 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ ANGPTL2, rna.data3)
                                      df$ANGPTL2 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ ANXA2P2, rna.data3)
                                      df$ANXA2P2 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ ANXA6, rna.data3)
                                      df$ANXA6 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ ANXA8, rna.data3)
                                      df$ANXA8 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ APLP1, rna.data3)
                                      df$APLP1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
                                    
cox_result = coxph(Surv(PFI.time, PFI) ~ BMP7, rna.data3)
df$BMP7 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ CCDC80, rna.data3)
df$CCDC80 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ CDH13, rna.data3)
df$CDH13 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ CDH2, rna.data3)
df$CDH2 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ COCH, rna.data3)
df$COCH = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ COL10A1, rna.data3)
df$COL10A1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ COL12A1, rna.data3)
df$COL12A1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ COL14A1, rna.data3)
df$COL14A1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ COL16A1, rna.data3)
df$COL16A1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ COL1A1, rna.data3)
df$COL1A1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ COL23A1, rna.data3)
df$COL23A1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ COL25A1, rna.data3)
df$COL25A1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])


cox_result = coxph(Surv(PFI.time, PFI) ~ COL2A1, rna.data3)
df$COL2A1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ COL4A1, rna.data3)
df$COL4A1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ COL4A2, rna.data3)
df$COL4A2 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ COL4A4, rna.data3)
df$COL4A4 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ COL4A5, rna.data3)
df$COL4A5 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ COL4A6, rna.data3)
df$COL4A6 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ COL5A1, rna.data3)
df$COL5A1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ COL6A2, rna.data3)
df$COL6A2 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ COL6A3, rna.data3)
df$COL6A3 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ COL7A1, rna.data3)
df$COL7A1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ COL8A1, rna.data3)
df$COL8A1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ COL9A2, rna.data3)
df$COL9A2 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])


cox_result = coxph(Surv(PFI.time, PFI) ~ CPA3, rna.data3)
df$CPA3 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ CSPG4, rna.data3)
df$CSPG4 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ CTSB, rna.data3)
df$CTSB = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ DCN, rna.data3)
df$DCN = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ DPT, rna.data3)
df$DPT = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ DST, rna.data3)
df$DST = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ EDIL3, rna.data3)
df$EDIL3 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ EFEMP1, rna.data3)
df$EFEMP1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ EFEMP2, rna.data3)
df$EFEMP2 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ EMILIN1, rna.data3)
df$EMILIN1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ ENTPD2, rna.data3)
df$ENTPD2 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ FBLN1, rna.data3)
df$FBLN1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])


cox_result = coxph(Surv(PFI.time, PFI) ~ FBLN2, rna.data3)
df$FBLN2 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ FBN1, rna.data3)
df$FBN1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ FBN2, rna.data3)
df$FBN2 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ FGFBP3, rna.data3)
df$FGFBP3 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ FGFR2, rna.data3)
df$FGFR2 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ FLG, rna.data3)
df$FLG = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ FN1, rna.data3)
df$FN1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ FRAS1, rna.data3)
df$FRAS1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ FREM2, rna.data3)
df$FREM2 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ GDF15, rna.data3)
df$GDF15 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ GPC5, rna.data3)
df$GPC5 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ GPC6, rna.data3)
df$GPC6 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])


cox_result = coxph(Surv(PFI.time, PFI) ~ HMCN1, rna.data3)
df$HMCN1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ KRT1, rna.data3)
df$KRT1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ L1CAM, rna.data3)
df$L1CAM = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ LAMA1, rna.data3)
df$LAMA1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ LAMA2, rna.data3)
df$LAMA2 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ LAMA3, rna.data3)
df$LAMA3 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ LAMB1, rna.data3)
df$LAMB1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ LAMB2, rna.data3)
df$LAMB2 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ LAMB4, rna.data3)
df$LAMB4 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ LAMC3, rna.data3)
df$LAMC3 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ LGALS1, rna.data3)
df$LGALS1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ LOXL2, rna.data3)
df$LOXL2 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])


cox_result = coxph(Surv(PFI.time, PFI) ~ LRRC15, rna.data3)
df$LRRC15 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ LTBP3, rna.data3)
df$LTBP3 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ LUM, rna.data3)
df$LUM = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ MATN3, rna.data3)
df$MATN3 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ MDK, rna.data3)
df$MDK = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ MFAP2, rna.data3)
df$MFAP2 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ MFAP5, rna.data3)
df$MFAP5 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ MGP, rna.data3)
df$MGP = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ MMP2, rna.data3)
df$MMP2 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ MMP9, rna.data3)
df$MMP9 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ MXRA5, rna.data3)
df$MXRA5 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ NAV2, rna.data3)
df$NAV2 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])


cox_result = coxph(Surv(PFI.time, PFI) ~ NCAM1, rna.data3)
df$NCAM1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ NID1, rna.data3)
df$NID1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ NID2, rna.data3)
df$NID2 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ NPNT, rna.data3)
df$NPNT = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ NTN1, rna.data3)
df$NTN1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ NTN4, rna.data3)
df$NTN4 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ NTNG1, rna.data3)
df$NTNG1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ NTNG2, rna.data3)
df$NTNG2 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ PCOLCE, rna.data3)
df$PCOLCE = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ PDGFB, rna.data3)
df$PDGFB = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ POSTN, rna.data3)
df$POSTN = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ RTBDN, rna.data3)
df$RTBDN = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])


cox_result = coxph(Surv(PFI.time, PFI) ~ S100A4, rna.data3)
df$S100A4 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ S100A6, rna.data3)
df$S100A6 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ S100A7, rna.data3)
df$S100A7 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ S100A8, rna.data3)
df$S100A8 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ S100A9, rna.data3)
df$S100A9 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ SDC2, rna.data3)
df$SDC2 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ SEMA3B, rna.data3)
df$SEMA3B = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ SEMA7A, rna.data3)
df$SEMA7A = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ SERPINB12, rna.data3)
df$SERPINB12 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ SERPINE1, rna.data3)
df$SERPINE1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ SERPINE2, rna.data3)
df$SERPINE2 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ SERPINF1, rna.data3)
df$SERPINF1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])


cox_result = coxph(Surv(PFI.time, PFI) ~ SERPING1, rna.data3)
df$SERPING1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ SFRP1, rna.data3)
df$SFRP1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ SHH, rna.data3)
df$SHH = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ SMOC1, rna.data3)
df$SMOC1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ SMOC2, rna.data3)
df$SMOC2 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ SOD3, rna.data3)
df$SOD3 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ SPARC, rna.data3)
df$SPARC = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ SRPX, rna.data3)
df$SRPX = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ SRPX2, rna.data3)
df$SRPX2 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ SSC5D, rna.data3)
df$SSC5D = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ TGFB1I1, rna.data3)
df$TGFB1I1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ TGFB2, rna.data3)
df$TGFB2 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])


cox_result = coxph(Surv(PFI.time, PFI) ~ TGM2, rna.data3)
df$TGM2 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ TGM4, rna.data3)
df$TGM4 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ THBS1, rna.data3)
df$THBS1 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ THBS2, rna.data3)
df$THBS2 = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ TNC, rna.data3)
df$TNC = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ VCAN, rna.data3)
df$VCAN = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ VIT, rna.data3)
df$VIT = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])
cox_result = coxph(Surv(PFI.time, PFI) ~ VWF, rna.data3)
df$VWF = c(summary(cox_result)[['coefficients']][c(2,5)],summary(cox_result)[['conf.int']][c(3,4)])

df= as.data.frame(t(df))
df = cbind(rownames(df),df)
write.xlsx(df,"coxph_TCGA.xlsx",)


cox_result = coxph(Surv(PFI.time, PFI) ~ COL1A1+FREM2+DPT+MATN3+COL10A1, rna.data3)
summary(cox_result)[['coefficients']][,c(2,5)]
summary(cox_result)[['conf.int']][,c(3,4)]

