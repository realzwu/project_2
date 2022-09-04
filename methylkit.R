setwd("")
library(methylKit)

## CpG calling / control batch effect / DMR detect 
file.list=list("bismark_extractor/DWSH8A_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DWSH8B_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DWSH8C_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DPSM3A_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DPSM3B_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DPSM3C_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DSH6A_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DSH6B_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DSH6C_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DNS11A_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DNS11B_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DNS11C_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DPCSTDA_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DPCSTDB_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DPCSTDC_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DPCSTCA_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DPCSTCB_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DPCSTCC_1_trimmed_bismark_hisat2.bismark.cov.gz")

file.list=list("bismark_extractor/DWSH8A_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DWSH8B_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DWSH8C_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DPSM3A_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DPSM3B_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DPSM3C_1_trimmed_bismark_hisat2.bismark.cov.gz")

file.list=list("bismark_extractor/DSH6A_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DSH6B_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DSH6C_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DNS11A_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DNS11B_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DNS11C_1_trimmed_bismark_hisat2.bismark.cov.gz")

file.list=list("bismark_extractor/DPCSTDA_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DPCSTDB_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DPCSTDC_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DPCSTCA_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DPCSTCB_1_trimmed_bismark_hisat2.bismark.cov.gz",
               "bismark_extractor/DPCSTCC_1_trimmed_bismark_hisat2.bismark.cov.gz")


# read the files to a methylRawList object: myobj
# myobj=methRead(file.list, sample.id=list("WPMY1-1","WPMY1-2","WPMY1-3","WPMY1-C1","WPMY1-C2","WPMY1-C3"),
#               assembly="hg38", treatment=c(rep(1,3),rep(0,3)), context="CpG", mincov = 10, pipeline = "bismarkCoverage")

myobj=methRead(file.list, sample.id=list("RWPE1-1","RWPE1-2","RWPE1-3","RWPE1-C1","RWPE1-C2","RWPE1-C3"),
               assembly="hg38", treatment=c(rep(1,3),rep(0,3)), context="CpG", mincov = 10, pipeline = "bismarkCoverage")

# myobj=methRead(file.list, sample.id=list("PC3-1","PC3-2","PC3-3","PC3-C1","PC3-C2","PC3-C3"),
#                assembly="hg38", treatment=c(rep(1,3),rep(0,3)), context="CpG", mincov = 10, pipeline = "bismarkCoverage")

for(i in 1:6){myobj[[i]]$chr = paste("chr",myobj[[i]]$chr,sep = "")}

# getMethylationStats(myobj[[2]], plot=TRUE, both.strands=FALSE)
# getCoverageStats(myobj[[2]], plot=TRUE, both.strands=FALSE)

# View(myobj[[2]])

filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)

meth=unite(filtered.myobj, destrand=FALSE)

# perc.meth=percMethylation(meth)
# mean(perc.meth[,1])
# mean(perc.meth[,2])

# meth.min=unite(myobj,min.per.group=1L)

# getCorrelation(meth,plot=FALSE)

# pdf("PCA_all.pdf")
clusterSamples(meth, dist="correlation", method="ward.D2", plot=T)
# dev.off()

PCASamples(meth, screeplot=TRUE)

PCASamples(meth, screeplot=FALSE)

# make some batch data frame
# this is a bogus data frame
# we don't have batch information for the example data

# sampleAnnotation=data.frame(batch_id=c("a","a","a","b","b","b"), age=c(19,34,19,19,23,40))

# as=assocComp(mBase=meth,sampleAnnotation)

# construct a new object by removing the first pricipal component
# from percent methylation value matrix
# newObj=removeComp(meth,comp=1)
# mat=percMethylation(meth)

# do some changes in the matrix
# this is just a toy example
# ideally you want to correct the matrix
# for batch effects
# mat[mat==100]=80

# reconstruct the methylBase from the corrected matrix
# newobj=reconstruct(mat,meth)

# myobj.low=methRead(file.list, sample.id=list("Stable1","Stable2","Stable3","Control1","Control2","Control3"),
#               assembly="hg38", treatment=c(1,1,1,0,0,0), context="CpG", mincov = 3, pipeline = "bismarkCoverage")

# tiles = tileMethylCounts(myobj.low,win.size=1000,step.size=1000,cov.bases = 10)
# head(tiles[[1]],3)

myDiff=calculateDiffMeth(meth,overdispersion="MN",test="Chisq",mc.cores=2)

# get all differentially methylated bases
myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")

# myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)

# diffMethPerChr(myDiff,plot=TRUE,qvalue.cutoff=0.01, meth.cutoff=25)




### Gene annotation

library(genomation)

gene.obj=genomation::readTranscriptFeatures("../ref_data/hg38_GENCODE.v38.bed")

genomation::annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)

cpg.obj=readFeatureFlank("../ref_data/CpGislands.Hsapiens.hg38.UCSC.bed", feature.flank.name=c("CpGi","shores"))

diffCpGann=annotateWithFeatureFlank(as(myDiff25p,"GRanges"), cpg.obj$CpGi,cpg.obj$shores,feature.name="CpGi",flank.name="shores")

# promoters=regionCounts(filtered.myobj,gene.obj$promoters)

# promoter1=promoters[[1]]
# head(promoter1)
# promoter2=promoters[[2]]
# head(promoter2)

diffAnn=annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)
# DiffTss=getAssociationWithTSS(diffAnn)

diffAnn1=annotateWithGeneParts(as(myDiff25p.hyper,"GRanges"),gene.obj)
DiffTss1=getAssociationWithTSS(diffAnn1)

diffAnn2=annotateWithGeneParts(as(myDiff25p.hypo,"GRanges"),gene.obj)
DiffTss2=getAssociationWithTSS(diffAnn2)

getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)
plotTargetAnnotation(diffAnn,precedence=TRUE, main="Differential Methylation of PC3")
getFeatsWithTargetsStats(diffAnn,percentage=TRUE)

plotTargetAnnotation(diffCpGann,col=c("#20639B","#3CAEA3","#F6D55C"), main=" Differential Methylation of PC3")




myDiff25p.hyper = cbind(target.row = c(1:length(rownames(myDiff25p.hyper))), myDiff25p.hyper)
myDiff25p.hyper = merge(DiffTss1, myDiff25p.hyper, by.x = "target.row", by.y = "target.row")

myDiff25p.hypo = cbind(target.row = c(1:length(rownames(myDiff25p.hypo))), myDiff25p.hypo)
myDiff25p.hypo = merge(DiffTss2, myDiff25p.hypo, by.x = "target.row", by.y = "target.row")

library(stringr)
splitEnsembl <- function(x){
  return(str_split(x[1],'[.]',simplify = T)[1])}
myDiff25p.hyper$feature.name <- sapply(myDiff25p.hyper$feature.name,splitEnsembl,simplify = T)
myDiff25p.hypo$feature.name <- sapply(myDiff25p.hypo$feature.name,splitEnsembl,simplify = T)

library(biomaRt)

mart <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')

BM1 = getBM(attributes = c('ensembl_transcript_id', 'hgnc_symbol'), filters = 'ensembl_transcript_id', values = myDiff25p.hyper$feature.name, mart = mart)
BM2 = getBM(attributes = c('ensembl_transcript_id', 'hgnc_symbol'), filters = 'ensembl_transcript_id', values = myDiff25p.hypo$feature.name, mart = mart)

counts_new1 = merge(myDiff25p.hyper,BM1,by.x = "feature.name",by.y = "ensembl_transcript_id")
counts_new2 = merge(myDiff25p.hypo,BM2,by.x = "feature.name",by.y = "ensembl_transcript_id")

library(openxlsx)
write.xlsx(counts_new1,"RWPE1_DEM_hyper.xlsx")
write.xlsx(counts_new2,"RWPE1_DEM_hypo.xlsx")
