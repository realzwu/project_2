# Figure 4 Heatmap of normalized expression profile of PC3, RWPE1 qnd WPMY1

setwd("")

library(openxlsx)
library(ComplexHeatmap)
library(rstatix)
library(circlize)

myfile <- "file/133genes.txt" # see attached files
genes <- read.delim(myfile)[,1]

raw.data = read.xlsx("featureCounts_tpm.xlsx",rowNames = F,colNames = T)
raw.data2 = raw.data[raw.data[,1]%in% genes,]
rownames(raw.data2) = raw.data2[,1]
raw.data2 = raw.data2[,-1]
# raw.data2 = raw.data2[,c(22:24,16:18,19:21,1:3,4:6,7:9)]

average = apply(raw.data2, 1, mean)
sd = apply(raw.data2, 1, sd)

# Normalization
for(i in 1:132) {raw.data2[i,] = (raw.data2[i,] - average[i])/ sd[i]}

raw.data3 = as.matrix(raw.data2)

## extract up/downregulation for annotation
corrlist = read.xlsx("files/133genes_separated",rowNames = T,colNames = F)
strlist = na.omit(as.character(corrlist[1,]))
epilist = na.omit(as.character(corrlist[2,]))
pc3list = na.omit(as.character(corrlist[3,]))

WPMY1 = rownames(raw.data2) 
RWPE1 = rownames(raw.data2) 
PC3 = rownames(raw.data2) 
for(i in 1:132){if(WPMY1[i] %in% strlist){if(mean(raw.data3[i,c(1:3)]) > mean(raw.data3[i,c(4:6)])){WPMY1[i] = "Up"} else{WPMY1[i] = "Down"}} else{WPMY1[i] = " "}}
for(i in 1:132){if(RWPE1[i] %in% epilist){if(mean(raw.data3[i,c(7:9)]) > mean(raw.data3[i,c(10:12)])){RWPE1[i] = "Up"} else{RWPE1[i] = "Down"}} else{RWPE1[i] = " "}}
for(i in 1:132){if(PC3[i] %in% pc3list){if(mean(raw.data3[i,c(13:15)]) < mean(raw.data3[i,c(16:18)])){PC3[i] = "Up"} else{PC3[i] = "Down"}} else{PC3[i] = " "}}

mean.value = apply(raw.data2,2,mean)
  
ha_column1 = rowAnnotation( PC3 = PC3, RWPE1 = RWPE1,  WPMY1 = WPMY1,col = list(WPMY1 = c("Up" = "#FF8849", "Down" = "#3DB7E4"), RWPE1 = c("Up" = "#FF8849", "Down" = "#3DB7E4"), PC3 = c("Up" = "#FF8849", "Down" = "#3DB7E4")))

## annotate cell lines
Cell = c(rep("WPMY1", 6), rep("RWPE1", 6), rep("PC3", 6))
# DKK3 = c(rep(c(rep("Low",3),rep("High",3)),3))
Sample = c(rep("DKK3 Silenced",3),rep("Control",3),rep("DKK3 Silenced",3),rep("Control",3),rep("Control",3),rep("DKK3 Inducted",3))
ha_column2 = HeatmapAnnotation( Cell = Cell, Sample = Sample,  col = list(Sample = c("DKK3 Inducted" = "#3DB7E4", "DKK3 Silenced" = "#FF8849", "Control" = "#69BE28"), Cell = c("WPMY1" = "#20639B", "RWPE1" = "#3CAEA3", "PC3" = "#F6D55C")))

# Dotplots
ha_column3 = HeatmapAnnotation(. = anno_points(apply(raw.data2,2,mean), height = unit(1, "cm")))

# annotate specific names (FBLN1, etc.)
ha_column4 = rowAnnotation(foo = anno_mark(at = c(6,12,19,23,42,49,58,71,74,77,78,96,126,128), labels = rownames(raw.data3)[c(6,12,19,23,42,49,58,71,74,77,78,96,126,128)]),annotation_name_gp= gpar(fontsize = 4))
  
col_fun = colorRamp2(c(-2.5, 0, 2.5), c("#40A5F4", "#000000", "#F97340"))
# fa = cluster_within_group(raw.data3,cell)

pdf("heatmap_pc3_2.pdf")

Heatmap(raw.data3, name = "Expression", row_title = "ECM-related Gene Sets", row_names_side = "right", row_dend_side = "left", 
        col = col_fun, column_title = " ", right_annotation = c(ha_column1,ha_column4), bottom_annotation = ha_column3 ,top_annotation = ha_column2, cluster_columns = F, column_dend_reorder = F, 
        row_dend_reorder = T, show_column_dend = F,column_km = 3, show_row_names = F, show_column_names = F)

dev.off()
