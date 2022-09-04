setwd("/Users/zd821/Documents/project2/mRNA/weihuanjing")

library(openxlsx)

uncleaned = read.xlsx("microenvironment_id_transformed.xlsx",sheet = 1,rowNames = F,colNames = T)

unclean_symbol = uncleaned[,2]

uncleaned2 = apply(uncleaned[,c(-1,-2)],2,as.numeric)

average = apply(uncleaned2, 1, mean)
standard = apply(uncleaned2, 1, sd)
for(i in 1:19869){uncleaned2[i,] = (uncleaned2[i,] - average[i]) / standard[i] }

cleaned = as.data.frame(cbind(unclean_symbol,uncleaned2))

write.xlsx(cleaned, "microenvironment_zscore.xlsx")

cleaned <- read.xlsx("microenvironment_zscore.xlsx",sheet = 1,rowNames = F,colNames = T)
cleaned = cleaned[-which(cleaned[,1] == "DKK3")[1],]

cleaned2 <- cleaned[!duplicated(cleaned[,1], fromLast=TRUE), ] 

write.xlsx(cleaned2[-2,], "microenvironment_zscore_cleaned.xlsx")

namelist = cleaned2[,1]
cleaned2 = cleaned2[,-1]
cleaned2 = cleaned2[,order(as.matrix(cleaned2[3,]),as.matrix(cleaned2[1,]))]
cleaned3 = apply(cleaned3[c(-1:-3),], 2 , as.numeric)

# PCA
# B P sB sP sT T
library(limma)

plotMDS(cleaned3, col = c(rep("#F5A700",25), rep("#F72000",25),rep("#0DC14B",25), rep("#8480CC",25),rep("#24207C",25), rep("#AA003C",25)), labels = c(rep("•",150)))



### DEG

rownames(cleaned2) = namelist
cleaned_B = cleaned2[c(-1:-3),cleaned2[3,] == "B"]
cleaned_P = cleaned2[c(-1:-3),cleaned2[3,] == "P"]
cleaned_T = cleaned2[c(-1:-3),cleaned2[3,] == "T"]
cleaned_sB = cleaned2[c(-1:-3),cleaned2[3,] == "sB"]
cleaned_sP = cleaned2[c(-1:-3),cleaned2[3,] == "sP"]
cleaned_sT = cleaned2[c(-1:-3),cleaned2[3,] == "sT"]

cleaned_B[,c(1:25)] = apply(cleaned_B, 2, as.numeric)
cleaned_P[,c(1:25)] = apply(cleaned_P, 2, as.numeric)
cleaned_T[,c(1:25)] = apply(cleaned_T, 2, as.numeric)
cleaned_sB[,c(1:25)] = apply(cleaned_sB, 2, as.numeric)
cleaned_sP[,c(1:25)] = apply(cleaned_sP, 2, as.numeric)
cleaned_sT[,c(1:25)] = apply(cleaned_sT, 2, as.numeric)

pairinfo = factor(c(rep("Eight",13),rep("Six",12)))
  
design <- model.matrix(~0+pairinfo)

colnames(design)=levels(pairinfo)

fit_B <- lmFit(cleaned_B, design)
fit_P <- lmFit(cleaned_P, design)
fit_T <- lmFit(cleaned_T, design)
fit_sB <- lmFit(cleaned_sB, design)
fit_sP <- lmFit(cleaned_sP, design)
fit_sT <- lmFit(cleaned_sT, design)

contrast.matrix<-makeContrasts(paste0(unique(pairinfo),collapse = "-"),levels = design)

fit2_B <- eBayes(contrasts.fit(fit_B, contrast.matrix)) ##这一步很重要，大家可以自行看看效果
fit2_P <- eBayes(contrasts.fit(fit_P, contrast.matrix))
fit2_T <- eBayes(contrasts.fit(fit_T, contrast.matrix))
fit2_sB <- eBayes(contrasts.fit(fit_sB, contrast.matrix))
fit2_sP <- eBayes(contrasts.fit(fit_sP, contrast.matrix))
fit2_sT <- eBayes(contrasts.fit(fit_sT, contrast.matrix))


tempOutput_B = topTable(fit2_B, coef=1, n=Inf, adjust.method = "BH")
tempOutput_P = topTable(fit2_P, coef=1, n=Inf, adjust.method = "BH")
tempOutput_T = topTable(fit2_T, coef=1, n=Inf, adjust.method = "BH")
tempOutput_sB = topTable(fit2_sB, coef=1, n=Inf, adjust.method = "BH")
tempOutput_sP = topTable(fit2_sP, coef=1, n=Inf, adjust.method = "BH")
tempOutput_sT = topTable(fit2_sT, coef=1, n=Inf, adjust.method = "BH")

nrDEG_B = na.omit(tempOutput_B) 
nrDEG_P = na.omit(tempOutput_P) 
nrDEG_T = na.omit(tempOutput_T) 
nrDEG_sB = na.omit(tempOutput_sB) 
nrDEG_sP = na.omit(tempOutput_sP) 
nrDEG_sT = na.omit(tempOutput_sT) 

nrDEG[,7] = rownames(nrDEG)
nrDEG[,7] = rownames(nrDEG)
nrDEG[,7] = rownames(nrDEG)
nrDEG[,7] = rownames(nrDEG)
nrDEG[,7] = rownames(nrDEG)
nrDEG[,7] = rownames(nrDEG)

Diff_B = topTable(fit2_B, coef=1, n = Inf,p.value = 0.05, adjust.method = "none")
Diff_P = topTable(fit2_P, coef=1, n = Inf,p.value = 0.05, adjust.method = "none")
Diff_T = topTable(fit2_T, coef=1, n = Inf,p.value = 0.05, adjust.method = "none")
Diff_sB = topTable(fit2_sB, coef=1, n = Inf,p.value = 0.05, adjust.method = "none")
Diff_sP = topTable(fit2_sP, coef=1, n = Inf,p.value = 0.05, adjust.method = "none")
Diff_sT = topTable(fit2_sT, coef=1, n = Inf,p.value = 0.05, adjust.method = "none")

Diff_B = cbind(rownames(Diff_B), Diff_B)
Diff_P = cbind(rownames(Diff_P), Diff_P)
Diff_T = cbind(rownames(Diff_T), Diff_T)
Diff_sB = cbind(rownames(Diff_sB), Diff_sB)
Diff_sP = cbind(rownames(Diff_sP), Diff_sP)
Diff_sT = cbind(rownames(Diff_sT), Diff_sT)

write.xlsx(Diff_B,"DEG_diffb.xlsx")
write.xlsx(Diff_P,"DEG_diffp.xlsx")
write.xlsx(Diff_T,"DEG_difft.xlsx")
write.xlsx(Diff_sB,"DEG_diffsb.xlsx")
write.xlsx(Diff_sP,"DEG_diffsp.xlsx")
write.xlsx(Diff_sT,"DEG_diffst.xlsx")

### merge ECM_list and DEG
corrlist = read.xlsx("corrplot_list.xlsx",rowNames = T,colNames = F)
list.sB = rownames(Diff_sB)[rownames(Diff_sB) %in% corrlist[1,]]
list.sP = rownames(Diff_sP)[rownames(Diff_sP) %in% corrlist[1,]]
list.sT = rownames(Diff_sT)[rownames(Diff_sT) %in% corrlist[1,]]
list.B = rownames(Diff_B)[rownames(Diff_B) %in% corrlist[2,]]
list.P1 = rownames(Diff_P)[rownames(Diff_P) %in% corrlist[2,]]
list.P2 = rownames(Diff_P)[rownames(Diff_P) %in% corrlist[3,]]
list.T = rownames(Diff_T)[rownames(Diff_T) %in% corrlist[3,]]

