### query gene expression dataset from GEO


setwd("/Users/zd821/Documents/project2/mRNA/weihuanjing")
# BiocManager::install("hugene10sttranscriptcluster.db")
# BiocManager::install("GEOquery")
# BiocManager::install("tidyverse")

rm(list=ls())
options(stringsAsFactors = F)
library(GEOquery)
library(dplyr)
library(tidyverse)
library(hugene10sttranscriptcluster.db)
library(openxlsx)

gset1 <- getGEO('GSE97284',destdir = '.',getGPL = F)
class(gset1)

gset <- gset1[[1]]


expr <- gset@assayData[["exprs"]]
exp <- as.data.frame(expr)
colnames(exp)

ids <- toTable(hugene10sttranscriptclusterSYMBOL) 

exp <- exp[ids$probe_id,]



## remove duplicate (should be careful)

exp <- cbind(ids,exp)
exp1 <- distinct(exp,symbol,.keep_all = T)

write.xlsx(exp,"weihuanjing_newer_id_transformed.xlsx")

