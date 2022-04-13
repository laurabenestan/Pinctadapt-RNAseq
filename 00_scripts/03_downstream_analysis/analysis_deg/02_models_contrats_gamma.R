#!/usr/bin/Rscript

ls()
rm(list=ls())
ls()

# Installing packages
#

library("plot3D")
# Loading packages

library('DESeq2')
library("ggplot2")
library("pheatmap")
library("RColorBrewer")
library(gridExtra)

#Global variables
setwd("analysis_deg/")



# Loading data
coldata<-read.table("../../01_info_files/design.txt", header=T,check.names=F)
coldata<-coldata[coldata$Time != 0,]
coldata$condition <- paste(coldata$Time,coldata$Genotype,coldata$Temperature,sep='_')

rownames(coldata)<-coldata$Ind
coldata$Ind<-NULL
coldata$Temperature_Time<-NULL
coldata <- coldata[order(rownames(coldata)),] 


#counts
cts<-read.table("join_global_genome_110621.txt", header=T,check.names=F)
rownames(cts)<-cts$genes
cts$genes<-NULL

cts<-cts[,colnames(cts) %in% rownames(coldata) ]
cts <- cts[,order(colnames(cts))] 

#check

rownames(coldata)
colnames(cts)
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
str(cts)
str(coldata)
#Create dds
#coldata$Temperature<-as.factor(coldata$Temperature) Polynomial
#coldata$Genotype<-NULL
#coldata$Temperature<-NULL
#coldata$Time<-NULL
coldata$condition<-as.factor(coldata$condition)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
### NEsted anova like https://support.bioconductor.org/p/107564/
name_filtered <-rownames(dds)

summary(colSums(cts, na.rm = FALSE, dims = 1))
colSums(cts, na.rm = FALSE, dims = 1)
mean(colSums(cts, na.rm = FALSE, dims = 1))

# filtering
dds <- estimateSizeFactors(dds)
dds$sizeFactor
plot(dds$sizeFactor)
idx <- rowSums(counts(dds,normalized=TRUE) >= 7 ) >= 8 # check use normalize
dds <- dds[idx,]
dim(dds)
# filtering
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)


#matrix log
norm.counts.host <- counts(dds, normalized=TRUE)
log.norm.counts.host <- log2(norm.counts.host + 2)

vsd.fast <- vst(dds, fitType='local',blind=FALSE)
vsd.assay<-assay(vst(dds, blind=TRUE))
mat.tot<-(as.matrix(vsd.assay))





## Perform analysis with contrasts
ddsTest<-dds
ddsTest$condition <- relevel(ddsTest$condition, ref = "1_Gambier_34")
ddsTest<- DESeq(ddsTest)
resultsNames(ddsTest)


T1_MAR_23vs34<-results(ddsTest, contrast=c("condition","1_Marquesas_23","1_Marquesas_34"))
T1_MAR_23vs34<-subset(T1_MAR_23vs34,padj<0.01)
T1_MAR_23vs34<-subset(T1_MAR_23vs34,abs(log2FoldChange)>1)
str(T1_MAR_23vs34@rownames)
write.table(T1_MAR_23vs34,file="contrasts/T1_MAR_23vs34.txt",quote=F)

ddsTest$condition <- relevel(ddsTest$condition, ref = "60_Marquesas_34")
ddsTest<- DESeq(ddsTest)
resultsNames(ddsTest)
resLFC<- lfcShrink(ddsTest, coef="condition_60_Marquesas_27_vs_60_Marquesas_34", type="apeglm")
resLFC_60_MAR_27_34<-resLFC
resLFC<-subset(resLFC,abs(log2FoldChange)>1)
resLFC<-subset(resLFC,padj<0.01)
write.table(resLFC,file="~/Desktop/temp.FC2.Mr.60_27vs34.txt",quote=F,sep="\t")

ddsTest$condition <- relevel(ddsTest$condition, ref = "60_Gambier_34")
ddsTest<- DESeq(ddsTest)
resultsNames(ddsTest)
resLFC<- lfcShrink(ddsTest, coef="condition_60_Gambier_27_vs_60_Gambier_34", type="apeglm")
resLFC_60_GAM_27_34<-resLFC
resLFC<-subset(resLFC,abs(log2FoldChange)>1)
resLFC<-subset(resLFC,padj<0.01)
write.table(resLFC,file="~/Desktop/temp.FC2.GAM.60_27vs34.txt",quote=F,sep="\t")

ddsTest$condition <- relevel(ddsTest$condition, ref = "1_Gambier_34")
ddsTest<- DESeq(ddsTest)
resultsNames(ddsTest)
resLFC<- lfcShrink(ddsTest, coef="condition_1_Gambier_27_vs_1_Gambier_34", type="apeglm")
resLFC_1_GAM_27_34<-resLFC
write.table(resLFC,file="~/Desktop/Gamma/revision_JAE/analysis_deg/contrasts/TOTAL.GAM.1_27vs34.txt",quote=F,sep="\t")
resLFC<-subset(resLFC,abs(log2FoldChange)>1)
resLFC<-subset(resLFC,padj<0.01)
summary(resLFC)
write.table(resLFC,file="~/Desktop/temp.FC2.GAM.1_27vs34.txt",quote=F,sep="\t")

## fix 1_Marquesas_34

ddsTest$condition <- relevel(ddsTest$condition, ref = "1_Marquesas_34")
ddsTest<- DESeq(ddsTest)
resultsNames(ddsTest)

resLFC<- lfcShrink(ddsTest, coef="condition_1_Marquesas_27_vs_1_Marquesas_34", type="apeglm")
resLFC_1_MAR_27_34<-resLFC
resLFC<-subset(resLFC,abs(log2FoldChange)>1)
resLFC<-subset(resLFC,padj<0.01)
summary(resLFC)
write.table(resLFC,file="~/Desktop/temp.FC2.1.27vs34.MAR.txt",quote=F,sep="\t")

T1_MAR_23vs34<-results(ddsTest, contrast=list("condition1_Marquesas_23","condition1_Marquesas_34"))
T1_MAR_23vs34<-subset(T1_MAR_23vs34,padj<0.01)
T1_MAR_23vs34<-subset(T1_MAR_23vs34,abs(log2FoldChange)>1)
str(T1_MAR_23vs34@rownames)

T1_MAR_23vs32<-results(ddsTest, contrast=c("condition","1_Marquesas_23","1_Marquesas_32"))
T1_MAR_23vs32<-subset(T1_MAR_23vs32,padj<0.01)
T1_MAR_23vs32<-subset(T1_MAR_23vs32,abs(log2FoldChange)>1)
str(T1_MAR_23vs32@rownames)
write.table(T1_MAR_23vs32,file="contrasts/T1_MAR_23vs32.txt",quote=F)

T1_MAR_23vs27<-results(ddsTest, contrast=c("condition","1_Marquesas_23","1_Marquesas_27"))
T1_MAR_23vs27<-subset(T1_MAR_23vs27,padj<0.01)
T1_MAR_23vs27<-subset(T1_MAR_23vs27,abs(log2FoldChange)>1)
write.table(T1_MAR_23vs27,file="contrasts/T1_MAR_23vs27.txt",quote=F)

T1_MAR_27vs32<-results(ddsTest, contrast=c("condition","1_Marquesas_27","1_Marquesas_32"))
T1_MAR_27vs32<-subset(T1_MAR_27vs32,padj<0.01)
T1_MAR_27vs32<-subset(T1_MAR_27vs32,abs(log2FoldChange)>1)
write.table(T1_MAR_27vs32,file="contrasts/T1_MAR_27vs32.txt",quote=F)

T1_MAR_27vs34<-results(ddsTest, contrast=c("condition","1_Marquesas_27","1_Marquesas_34"))
resLFC_1_MAR_27_34<-T1_MAR_27vs34
T1_MAR_27vs34<-subset(T1_MAR_27vs34,padj<0.01)
T1_MAR_27vs34<-subset(T1_MAR_27vs34,abs(log2FoldChange)>1)
write.table(T1_MAR_27vs34,file="contrasts/T1_MAR_27vs34.txt",quote=F)
write.table(T1_MAR_27vs34,file="~/Desktop/T1_MAR_27vs34_fdr.txt",quote=F)


                  
T1_MAR_32vs34<-results(ddsTest, contrast=c("condition","1_Marquesas_32","1_Marquesas_34"))
T1_MAR_32vs34<-subset(T1_MAR_32vs34,padj<0.01)
T1_MAR_32vs34<-subset(T1_MAR_32vs34,abs(log2FoldChange)>1)
write.table(T1_MAR_32vs34,file="contrasts/T1_MAR_32vs34.txt",quote=F)

## T1 for gambier
T1_GAM_23vs34<-results(ddsTest, contrast=c("condition","1_Gambier_23","1_Gambier_34"))
T1_GAM_23vs34<-subset(T1_GAM_23vs34,padj<0.01)
T1_GAM_23vs34<-subset(T1_GAM_23vs34,abs(log2FoldChange)>1)
write.table(T1_GAM_23vs34,file="contrasts/T1_GAM_23vs34.txt",quote=F)

T1_GAM_23vs32<-results(ddsTest, contrast=c("condition","1_Gambier_23","1_Gambier_32"))
T1_GAM_23vs32<-subset(T1_GAM_23vs32,padj<0.01)
T1_GAM_23vs32<-subset(T1_GAM_23vs32,abs(log2FoldChange)>1)
write.table(T1_GAM_23vs32,file="contrasts/T1_GAM_23vs32.txt",quote=F)

T1_GAM_23vs27<-results(ddsTest, contrast=c("condition","1_Gambier_23","1_Gambier_27"))
T1_GAM_23vs27<-subset(T1_GAM_23vs27,padj<0.01)
T1_GAM_23vs27<-subset(T1_GAM_23vs27,abs(log2FoldChange)>1)
write.table(T1_GAM_23vs27,file="contrasts/T1_GAM_23vs27.txt",quote=F)

T1_GAM_27vs32<-results(ddsTest, contrast=c("condition","1_Gambier_27","1_Gambier_32"))
T1_GAM_27vs32<-subset(T1_GAM_27vs32,padj<0.01)
T1_GAM_27vs32<-subset(T1_GAM_27vs32,abs(log2FoldChange)>1)
write.table(T1_GAM_27vs32,file="contrasts/T1_GAM_27vs32.txt",quote=F) ###

T1_GAM_27vs34<-results(ddsTest, contrast=c("condition","1_Gambier_27","1_Gambier_34")) #Gambier 2d 27°C vs 34°C
T1_GAM_27vs34<-subset(T1_GAM_27vs34,padj<0.01)
T1_GAM_27vs34<-subset(T1_GAM_27vs34,abs(log2FoldChange)>1)
write.table(T1_GAM_27vs34,file="~/Desktop/T1_GAM_27vs34_fdr.txt",quote=F)


T1_GAM_32vs34<-results(ddsTest, contrast=c("condition","1_Gambier_32","1_Gambier_34"))
T1_GAM_32vs34<-subset(T1_GAM_32vs34,padj<0.01)
T1_GAM_32vs34<-subset(T1_GAM_32vs34,abs(log2FoldChange)>1)
write.table(T1_GAM_32vs34,file="contrasts/T1_GAM_32vs34.txt",quote=F)

## T60 MAR
T60_MAR_23vs34<-results(ddsTest, contrast=c("condition","60_Marquesas_23","60_Marquesas_34"))
T60_MAR_23vs34<-subset(T60_MAR_23vs34,padj<0.01)
T60_MAR_23vs34<-subset(T60_MAR_23vs34,abs(log2FoldChange)>1)
write.table(T60_MAR_23vs34,file="contrasts/T60_MAR_23vs34.txt",quote=F)

T60_MAR_23vs32<-results(ddsTest, contrast=c("condition","60_Marquesas_23","60_Marquesas_32"))
T60_MAR_23vs32<-subset(T60_MAR_23vs32,padj<0.01)
T60_MAR_23vs32<-subset(T60_MAR_23vs32,abs(log2FoldChange)>1)
write.table(T60_MAR_23vs32,file="contrasts/T60_MAR_23vs32.txt",quote=F)

T60_MAR_23vs27<-results(ddsTest, contrast=c("condition","60_Marquesas_23","60_Marquesas_27"))
T60_MAR_23vs27<-subset(T60_MAR_23vs27,padj<0.01)
T60_MAR_23vs27<-subset(T60_MAR_23vs27,abs(log2FoldChange)>1)
write.table(T60_MAR_23vs27,file="contrasts/T60_MAR_23vs27.txt",quote=F)

T60_MAR_27vs32<-results(ddsTest, contrast=c("condition","60_Marquesas_27","60_Marquesas_32"))
T60_MAR_27vs32<-subset(T60_MAR_27vs32,padj<0.01)
T60_MAR_27vs32<-subset(T60_MAR_27vs32,abs(log2FoldChange)>1)
write.table(T60_MAR_27vs32,file="contrasts/T60_MAR_27vs32.txt",quote=F)

T60_MAR_27vs34<-results(ddsTest, contrast=c("condition","60_Marquesas_27","60_Marquesas_34"))
T60_MAR_27vs34<-subset(T60_MAR_27vs34,padj<0.01)
T60_MAR_27vs34<-subset(T60_MAR_27vs34,abs(log2FoldChange)>1)
str(T60_MAR_27vs34@rownames)
write.table(T60_MAR_27vs34,file="contrasts/T60_MAR_27vs34_full.txt",quote=F) 
write.table(T60_MAR_27vs34,file="~/Desktop/T60_MAR_27vs34_fdr.txt",quote=F)

T60_MAR_32vs34<-results(ddsTest, contrast=c("condition","60_Marquesas_32","60_Marquesas_34"))
T60_MAR_32vs34<-subset(T60_MAR_32vs34,padj<0.01)
T60_MAR_32vs34<-subset(T60_MAR_32vs34,abs(log2FoldChange)>1)
write.table(T60_MAR_32vs34,file="contrasts/T60_MAR_32vs34.txt",quote=F)

## T60 for gambier
T60_GAM_23vs34<-results(ddsTest, contrast=c("condition","60_Gambier_23","60_Gambier_34"))
T60_GAM_23vs34<-subset(T60_GAM_23vs34,padj<0.01)
T60_GAM_23vs34<-subset(T60_GAM_23vs34,abs(log2FoldChange)>1)
write.table(T60_GAM_23vs34,file="contrasts/T60_GAM_23vs34.txt",quote=F)

T60_GAM_23vs32<-results(ddsTest, contrast=c("condition","60_Gambier_23","60_Gambier_32"))
T60_GAM_23vs32<-subset(T60_GAM_23vs32,padj<0.01)
T60_GAM_23vs32<-subset(T60_GAM_23vs32,abs(log2FoldChange)>1)
write.table(T60_GAM_23vs32,file="contrasts/T60_GAM_23vs32.txt",quote=F)

T60_GAM_23vs27<-results(ddsTest, contrast=c("condition","60_Gambier_23","60_Gambier_27"))
T60_GAM_23vs27<-subset(T60_GAM_23vs27,padj<0.01)
T60_GAM_23vs27<-subset(T60_GAM_23vs27,abs(log2FoldChange)>1)
write.table(T60_GAM_23vs27,file="contrasts/T60_GAM_23vs27.txt",quote=F)

T60_GAM_27vs32<-results(ddsTest, contrast=c("condition","60_Gambier_27","60_Gambier_32"))
T60_GAM_27vs32<-subset(T60_GAM_27vs32,padj<0.01)
T60_GAM_27vs32<-subset(T60_GAM_27vs32,abs(log2FoldChange)>1)
write.table(T60_GAM_27vs32,file="contrasts/T60_GAM_27vs32.txt",quote=F)

T60_GAM_27vs34<-results(ddsTest, contrast=c("condition","60_Gambier_27","60_Gambier_34"))
T60_GAM_27vs34<-subset(T60_GAM_27vs34,padj<0.01)
T60_GAM_27vs34<-subset(T60_GAM_27vs34,abs(log2FoldChange)>1)
write.table(T60_GAM_27vs34,file="~/Desktop/T60_GAM_27vs34_fdr.txt",quote=F)

T60_GAM_32vs34<-results(ddsTest, contrast=c("condition","60_Gambier_32","60_Gambier_34"))
T60_GAM_32vs34<-subset(T60_GAM_32vs34,padj<0.01)
T60_GAM_32vs34<-subset(T60_GAM_32vs34,abs(log2FoldChange)>1)
write.table(T60_GAM_32vs34,file="contrasts/T60_GAM_32vs34.txt",quote=F)


### Fix marquesas Gambier ## fix 1_Marquesas_34

ddsTest$condition <- relevel(ddsTest$condition, ref = "1_Marquesas_27")
ddsTest<- DESeq(ddsTest)
resultsNames(ddsTest)
## Combination t1 vs T1 MAR- GAM 27°C
T1vs1_GAM_MAR_27vs27<-results(ddsTest, contrast=c("condition","1_Marquesas_27","1_Gambier_27"))
T1vs1_GAM_MAR_27vs27<-subset(T1vs1_GAM_MAR_27vs27,padj<0.01)
T1vs1_GAM_MAR_27vs27<-subset(T1vs1_GAM_MAR_27vs27,abs(log2FoldChange)>1)
write.table(T1vs1_GAM_MAR_27vs27,file="contrasts/T1vs1_GAM_MAR_27vs27.total.txt",quote=F)

## Combination t1 vs T60 MAR
T1vs60_MAR_23vs23<-results(ddsTest, contrast=c("condition","1_Marquesas_23","60_Marquesas_23"))
T1vs60_MAR_23vs23<-subset(T1vs60_MAR_23vs23,padj<0.01)
T1vs60_MAR_23vs23<-subset(T1vs60_MAR_23vs23,abs(log2FoldChange)>1)
write.table(T1vs60_MAR_23vs23,file="contrasts/T1vs60_MAR_23vs23.txt",quote=F)

T1vs60_MAR_34vs34<-results(ddsTest, contrast=c("condition","1_Marquesas_34","60_Marquesas_34"))
T1vs60_MAR_34vs34<-subset(T1vs60_MAR_34vs34,padj<0.01)
T1vs60_MAR_34vs34<-subset(T1vs60_MAR_34vs34,abs(log2FoldChange)>1)
write.table(T1vs60_MAR_34vs34,file="contrasts/T1vs60_MAR_34vs34.txt",quote=F)
T1vs60_MAR_34vs34.up<-subset(T1vs60_MAR_34vs34,log2FoldChange > 0)
list_T1vs60_MAR_34vs34.up<-T1vs60_MAR_34vs34.up@rownames
T1vs60_MAR_34vs34.down<-subset(T1vs60_MAR_34vs34,log2FoldChange < 0)
list_T1vs60_MAR_34vs34.down<-T1vs60_MAR_34vs34.down@rownames    
write.table(list_T1vs60_MAR_34vs34.up,file="contrasts/T1vs60_MAR_34vs34.up.txt",quote=F,col.names = F,row.names = F)
write.table(list_T1vs60_MAR_34vs34.down,file="contrasts/T1vs60_MAR_34vs34.down.txt",quote=F,col.names = F,row.names = F)


T1vs60_MAR_32vs32<-results(ddsTest, contrast=c("condition","1_Marquesas_32","60_Marquesas_32"))
T1vs60_MAR_32vs32<-subset(T1vs60_MAR_32vs32,padj<0.01)
T1vs60_MAR_32vs32<-subset(T1vs60_MAR_32vs32,abs(log2FoldChange)>1)
write.table(T1_MAR_23vs34,file="contrasts/T1_MAR_23vs34.txt",quote=F)

T1vs60_MAR_27vs27<-results(ddsTest, contrast=c("condition","1_Marquesas_27","60_Marquesas_27"))
T1vs60_MAR_27vs27<-subset(T1vs60_MAR_27vs27,padj<0.01)
T1vs60_MAR_27vs27<-subset(T1vs60_MAR_27vs27,abs(log2FoldChange)>1)
write.table(T1vs60_MAR_27vs27,file="contrasts/T1vs60_MAR_27vs27.txt",quote=F)


T1vs60_MAR_23vs27<-results(ddsTest, contrast=c("condition","1_Marquesas_23","60_Marquesas_27"))
T1vs60_MAR_23vs27<-subset(T1vs60_MAR_23vs27,padj<0.01)
T1vs60_MAR_23vs27<-subset(T1vs60_MAR_23vs27,abs(log2FoldChange)>1)
write.table(T1vs60_MAR_23vs27,file="contrasts/T1vs60_MAR_23vs27.txt",quote=F)

T1vs60_MAR_23vs32<-results(ddsTest, contrast=c("condition","1_Marquesas_23","60_Marquesas_32"))
T1vs60_MAR_23vs32<-subset(T1vs60_MAR_23vs32,padj<0.01)
T1vs60_MAR_23vs32<-subset(T1vs60_MAR_23vs32,abs(log2FoldChange)>1)
write.table(T1vs60_MAR_23vs32,file="contrasts/T1vs60_MAR_23vs32.txt",quote=F)

T1vs60_MAR_23vs34<-results(ddsTest, contrast=c("condition","1_Marquesas_23","60_Marquesas_34"))
T1vs60_MAR_23vs34<-subset(T1vs60_MAR_23vs34,padj<0.01)
T1vs60_MAR_23vs34<-subset(T1vs60_MAR_23vs34,abs(log2FoldChange)>1)
write.table(T1vs60_MAR_23vs34,file="contrasts/T1vs60_MAR_23vs34.txt",quote=F)




T1vs60_MAR_27vs32<-results(ddsTest, contrast=c("condition","1_Marquesas_27","60_Marquesas_32"))
T1vs60_MAR_27vs32<-subset(T1vs60_MAR_27vs32,padj<0.01)
T1vs60_MAR_27vs32<-subset(T1vs60_MAR_27vs32,abs(log2FoldChange)>1)
write.table(T1vs60_MAR_27vs32,file="contrasts/T1vs60_MAR_27vs32.txt",quote=F)

T1vs60_MAR_27vs34<-results(ddsTest, contrast=c("condition","1_Marquesas_27","60_Marquesas_34"))
T1vs60_MAR_27vs34<-subset(T1vs60_MAR_27vs34,padj<0.01)
T1vs60_MAR_27vs34<-subset(T1vs60_MAR_27vs34,abs(log2FoldChange)>1)
write.table(T1vs60_MAR_27vs34,file="contrasts/T1vs60_MAR_27vs34.txt",quote=F)

T1vs60_MAR_27vs23<-results(ddsTest, contrast=c("condition","1_Marquesas_27","60_Marquesas_23"))
T1vs60_MAR_27vs23<-subset(T1vs60_MAR_27vs23,padj<0.01)
T1vs60_MAR_27vs23<-subset(T1vs60_MAR_27vs23,abs(log2FoldChange)>1)
write.table(T1vs60_MAR_27vs23,file="contrasts/T1vs60_MAR_27vs23.txt",quote=F)

T1vs60_MAR_32vs23<-results(ddsTest, contrast=c("condition","1_Marquesas_32","60_Marquesas_23"))
T1vs60_MAR_32vs23<-subset(T1vs60_MAR_32vs23,padj<0.01)
T1vs60_MAR_32vs23<-subset(T1vs60_MAR_32vs23,abs(log2FoldChange)>1)
write.table(T1vs60_MAR_32vs23,file="contrasts/T1vs60_MAR_32vs23.txt",quote=F)

T1vs60_MAR_32vs27<-results(ddsTest, contrast=c("condition","1_Marquesas_32","60_Marquesas_27"))
T1vs60_MAR_32vs27<-subset(T1vs60_MAR_32vs27,padj<0.01)
T1vs60_MAR_32vs27<-subset(T1vs60_MAR_32vs27,abs(log2FoldChange)>1)
write.table(T1vs60_MAR_32vs27,file="contrasts/T1vs60_MAR_32vs27.txt",quote=F)

T1vs60_MAR_32vs34<-results(ddsTest, contrast=c("condition","1_Marquesas_32","60_Marquesas_34"))
T1vs60_MAR_32vs34<-subset(T1vs60_MAR_32vs34,padj<0.01)
T1vs60_MAR_32vs34<-subset(T1vs60_MAR_32vs34,abs(log2FoldChange)>1)
write.table(T1vs60_MAR_32vs34,file="contrasts/T1vs60_MAR_32vs34.txt",quote=F)

T1vs60_MAR_32vs32<-results(ddsTest, contrast=c("condition","1_Marquesas_32","60_Marquesas_32"))
T1vs60_MAR_32vs32<-subset(T1vs60_MAR_32vs32,padj<0.01)
T1vs60_MAR_32vs32<-subset(T1vs60_MAR_32vs32,abs(log2FoldChange)>1)
write.table(T1vs60_MAR_32vs32,file="contrasts/T1vs60_MAR_32vs32.txt",quote=F)

T1vs60_MAR_34vs23<-results(ddsTest, contrast=c("condition","1_Marquesas_34","60_Marquesas_23"))
T1vs60_MAR_34vs23<-subset(T1vs60_MAR_34vs23,padj<0.01)
T1vs60_MAR_34vs23<-subset(T1vs60_MAR_34vs23,abs(log2FoldChange)>1)
write.table(T1vs60_MAR_34vs23,file="contrasts/T1vs60_MAR_34vs23.txt",quote=F)

T1vs60_MAR_34vs27<-results(ddsTest, contrast=c("condition","1_Marquesas_34","60_Marquesas_27"))
T1vs60_MAR_34vs27<-subset(T1vs60_MAR_34vs27,padj<0.01)
T1vs60_MAR_34vs27<-subset(T1vs60_MAR_34vs27,abs(log2FoldChange)>1)
write.table(T1vs60_MAR_34vs27,file="contrasts/T1vs60_MAR_34vs27.txt",quote=F)

T1vs60_MAR_34vs32<-results(ddsTest, contrast=c("condition","1_Marquesas_34","60_Marquesas_32"))
T1vs60_MAR_34vs32<-subset(T1vs60_MAR_34vs32,padj<0.01)
T1vs60_MAR_34vs32<-subset(T1vs60_MAR_34vs32,abs(log2FoldChange)>1)
write.table(T1vs60_MAR_34vs32,file="contrasts/T1vs60_MAR_34vs32.txt",quote=F)


## T1 vs 60 Gambier
T1vs60_GAM_23vs23<-results(ddsTest, contrast=c("condition","1_Gambier_23","60_Gambier_23"))
T1vs60_GAM_23vs23<-subset(T1vs60_GAM_23vs23,padj<0.01)
T1vs60_GAM_23vs23<-subset(T1vs60_GAM_23vs23,abs(log2FoldChange)>1)
write.table(T1vs60_GAM_23vs23,file="contrasts/T1vs60_GAM_23vs23.txt",quote=F)

T1vs60_GAM_34vs34<-results(ddsTest, contrast=c("condition","1_Gambier_34","60_Gambier_34"))
T1vs60_GAM_34vs34<-subset(T1vs60_GAM_34vs34,padj<0.01)
T1vs60_GAM_34vs34<-subset(T1vs60_GAM_34vs34,abs(log2FoldChange)>1)
write.table(T1vs60_GAM_34vs34,file="contrasts/T1vs60_GAM_34vs34.txt",quote=F)
T1vs60_GAM_34vs34.up<-subset(T1vs60_GAM_34vs34,log2FoldChange > 0)
list_T1vs60_GAM_34vs34.up<-T1vs60_GAM_34vs34.up@rownames
T1vs60_GAM_34vs34.down<-subset(T1vs60_GAM_34vs34,log2FoldChange < 0)
list_T1vs60_GAM_34vs34.down<-T1vs60_GAM_34vs34.down@rownames                            
write.table(list_T1vs60_GAM_34vs34.up,file="contrasts/T1vs60_GAM_34vs34.up.txt",quote=F,col.names = F,row.names = F)
write.table(list_T1vs60_GAM_34vs34.down,file="contrasts/T1vs60_GAM_34vs34.down.txt",quote=F,col.names = F,row.names = F)


T1vs60_GAM_32vs32<-results(ddsTest, contrast=c("condition","1_Gambier_32","60_Gambier_32"))
T1vs60_GAM_32vs32<-subset(T1vs60_GAM_32vs32,padj<0.01)
T1vs60_GAM_32vs32<-subset(T1vs60_GAM_32vs32,abs(log2FoldChange)>1)
write.table(T1vs60_GAM_32vs32,file="contrasts/T1vs60_GAM_32vs32.txt",quote=F)

T1vs60_GAM_27vs27<-results(ddsTest, contrast=c("condition","1_Gambier_27","60_Gambier_27"))
T1vs60_GAM_27vs27<-subset(T1vs60_GAM_27vs27,padj<0.01)
T1vs60_GAM_27vs27<-subset(T1vs60_GAM_27vs27,abs(log2FoldChange)>1)
write.table(T1vs60_GAM_27vs27,file="contrasts/T1vs60_GAM_27vs27.txt",quote=F)


T1vs60_GAM_23vs27<-results(ddsTest, contrast=c("condition","1_Gambier_23","60_Gambier_27"))
T1vs60_GAM_23vs27<-subset(T1vs60_GAM_23vs27,padj<0.01)
T1vs60_GAM_23vs27<-subset(T1vs60_GAM_23vs27,abs(log2FoldChange)>1)
write.table(T1vs60_GAM_23vs27,file="contrasts/T1vs60_GAM_23vs27.txt",quote=F)

T1vs60_GAM_23vs32<-results(ddsTest, contrast=c("condition","1_Gambier_23","60_Gambier_32"))
T1vs60_GAM_23vs32<-subset(T1vs60_GAM_23vs32,padj<0.01)
T1vs60_GAM_23vs32<-subset(T1vs60_GAM_23vs32,abs(log2FoldChange)>1)
write.table(T1vs60_GAM_23vs32,file="contrasts/T1vs60_GAM_23vs32.txt",quote=F)

T1vs60_GAM_23vs34<-results(ddsTest, contrast=c("condition","1_Gambier_23","60_Gambier_34"))
T1vs60_GAM_23vs34<-subset(T1vs60_GAM_23vs34,padj<0.01)
T1vs60_GAM_23vs34<-subset(T1vs60_GAM_23vs34,abs(log2FoldChange)>1)
write.table(T1vs60_GAM_23vs34,file="contrasts/T1vs60_GAM_23vs34.txt",quote=F)



T1vs60_GAM_27vs32<-results(ddsTest, contrast=c("condition","1_Gambier_27","60_Gambier_32"))
T1vs60_GAM_27vs32<-subset(T1vs60_GAM_27vs32,padj<0.01)
T1vs60_GAM_27vs32<-subset(T1vs60_GAM_27vs32,abs(log2FoldChange)>1)
write.table(T1vs60_GAM_27vs32,file="contrasts/T1vs60_GAM_27vs32.txt",quote=F)

T1vs60_GAM_27vs34<-results(ddsTest, contrast=c("condition","1_Gambier_27","60_Gambier_34"))
T1vs60_GAM_27vs34<-subset(T1vs60_GAM_27vs34,padj<0.01)
T1vs60_GAM_27vs34<-subset(T1vs60_GAM_27vs34,abs(log2FoldChange)>1)
write.table(T1vs60_GAM_27vs34,file="contrasts/T1vs60_GAM_27vs34.txt",quote=F)

T1vs60_GAM_27vs23<-results(ddsTest, contrast=c("condition","1_Gambier_27","60_Gambier_23"))
T1vs60_GAM_27vs23<-subset(T1vs60_GAM_27vs23,padj<0.01)
T1vs60_GAM_27vs23<-subset(T1vs60_GAM_27vs23,abs(log2FoldChange)>1)
write.table(T1vs60_GAM_27vs23,file="contrasts/T1vs60_GAM_27vs23.txt",quote=F)


T1vs60_GAM_32vs23<-results(ddsTest, contrast=c("condition","1_Gambier_32","60_Gambier_23"))
T1vs60_GAM_32vs23<-subset(T1vs60_GAM_32vs23,padj<0.01)
T1vs60_GAM_32vs23<-subset(T1vs60_GAM_32vs23,abs(log2FoldChange)>1)
write.table(T1vs60_GAM_32vs23,file="contrasts/T1vs60_GAM_32vs23.txt",quote=F)

T1vs60_GAM_32vs27<-results(ddsTest, contrast=c("condition","1_Gambier_32","60_Gambier_27"))
T1vs60_GAM_32vs27<-subset(T1vs60_GAM_32vs27,padj<0.01)
T1vs60_GAM_32vs27<-subset(T1vs60_GAM_32vs27,abs(log2FoldChange)>1)
write.table(T1vs60_GAM_32vs27,file="contrasts/T1vs60_GAM_32vs27.txt",quote=F)

T1vs60_GAM_32vs34<-results(ddsTest, contrast=c("condition","1_Gambier_32","60_Gambier_34"))
T1vs60_GAM_32vs34<-subset(T1vs60_GAM_32vs34,padj<0.01)
T1vs60_GAM_32vs34<-subset(T1vs60_GAM_32vs34,abs(log2FoldChange)>1)
write.table(T1vs60_GAM_32vs34,file="contrasts/T1vs60_GAM_32vs34.txt",quote=F)


T1vs60_GAM_34vs23<-results(ddsTest, contrast=c("condition","1_Gambier_34","60_Gambier_23"))
T1vs60_GAM_34vs23<-subset(T1vs60_GAM_34vs23,padj<0.01)
T1vs60_GAM_34vs23<-subset(T1vs60_GAM_34vs23,abs(log2FoldChange)>1)
write.table(T1vs60_GAM_34vs23,file="contrasts/T1vs60_GAM_34vs23.txt",quote=F)

T1vs60_GAM_34vs27<-results(ddsTest, contrast=c("condition","1_Gambier_34","60_Gambier_27"))
T1vs60_GAM_34vs27<-subset(T1vs60_GAM_342vs27,padj<0.01)
T1vs60_GAM_34vs27<-subset(T1vs60_GAM_34vs27,abs(log2FoldChange)>1)
write.table(T1vs60_GAM_34vs27,file="contrasts/T1vs60_GAM_34vs27.txt",quote=F)

T1vs60_GAM_34vs32<-results(ddsTest, contrast=c("condition","1_Gambier_34","60_Gambier_32"))
T1vs60_GAM_34vs32<-subset(T1vs60_GAM_34vs32,padj<0.01)
T1vs60_GAM_34vs32<-subset(T1vs60_GAM_34vs32,abs(log2FoldChange)>1)
write.table(T1vs60_GAM_34vs32,file="contrasts/T1vs60_GAM_34vs32.txt",quote=F)


T1_27_GAM_vs_MAR<-results(ddsTest, contrast=c("condition","1_Gambier_27","1_Marquesas_27"))
T1_27_GAM_vs_MAR<-subset(T1_27_GAM_vs_MAR,padj<0.05)
T1_27_GAM_vs_MAR<-subset(T1_27_GAM_vs_MAR,abs(log2FoldChange)>1)
write.table(T1_27_GAM_vs_MAR,file="contrasts/T1_27_GAM_vs_MAR.txt",quote=F)

T1_34_GAM_vs_MAR<-results(ddsTest, contrast=c("condition","1_Gambier_34","1_Marquesas_34"))
T1_34_GAM_vs_MAR<-subset(T1_34_GAM_vs_MAR,padj<0.05)
T1_34_GAM_vs_MAR<-subset(T1_34_GAM_vs_MAR,abs(log2FoldChange)>1)
write.table(T1_34_GAM_vs_MAR,file="contrasts/T1_34_GAM_vs_MAR.txt",quote=F)



#exctraction GEI by time using LRT
## Time 1
# Loading data
coldata<-read.table("design.txt", header=T,check.names=F)
coldata<-coldata[coldata$Time == 1,]

rownames(coldata)<-coldata$Ind
coldata$Ind<-NULL
coldata$Temperature_Time<-NULL
coldata <- coldata[order(rownames(coldata)),] 

#counts
cts<-read.table("join_gamma.txt", header=T,check.names=F)
rownames(cts)<-cts$genes
cts$genes<-NULL

cts<-cts[,colnames(cts) %in% rownames(coldata) ]
cts <- cts[,order(colnames(cts))] 

#clean coldata
vectorname<-colnames(cts)
coldata<-coldata[rownames(coldata) %in% vectorname, ]
conta=as.character(c("TRINITY_DN63165_c1_g2_i6", "TRINITY_DN88404_c13_g1_i8"))
cts<-cts[!rownames(cts) %in% conta, ]
# Sequences major proportion

#check

rownames(coldata)
colnames(cts)
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
str(cts)
str(coldata)


coldata$Genotype<-as.factor(coldata$Genotype)
coldata$Temperature<-as.factor(coldata$Temperature)
coldata$Time<-NULL


ddsLRT <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Genotype+Temperature+Genotype*Temperature)
name_filtered <-rownames(dds)


# filtering
ddsLRT <- estimateSizeFactors(ddsLRT)
ddsLRT$sizeFactor
plot(ddsLRT$sizeFactor)
idx <- rowSums(counts(ddsLRT,normalized=TRUE) >= 7 ) >= 4 # check use normalize
ddsLRT <- ddsLRT[idx,]
dim(ddsLRT)
# filtering
ddsLRT <- estimateSizeFactors(ddsLRT)
ddsLRT <- estimateDispersions(ddsLRT)

#test against reduced model
ddsLRT <- DESeq(ddsLRT, test="LRT", reduced=~Genotype + Temperature)
resGEIT1 <- results(ddsLRT)
resGEIT1<-subset(resGEIT1,padj<0.05)
write.table(resGEIT1,file="contrasts/T1_GEI.txt",quote=F)

## Time 60
# Loading data
coldata<-read.table("design.txt", header=T,check.names=F)
coldata<-coldata[coldata$Time == 60,]

rownames(coldata)<-coldata$Ind
coldata$Ind<-NULL
coldata$Temperature_Time<-NULL
coldata <- coldata[order(rownames(coldata)),] 

#counts
cts<-read.table("join_gamma.txt", header=T,check.names=F)
rownames(cts)<-cts$genes
cts$genes<-NULL

cts<-cts[,colnames(cts) %in% rownames(coldata) ]
cts <- cts[,order(colnames(cts))] 

#clean coldata
vectorname<-colnames(cts)
coldata<-coldata[rownames(coldata) %in% vectorname, ]
conta=as.character(c("TRINITY_DN63165_c1_g2_i6", "TRINITY_DN88404_c13_g1_i8"))
cts<-cts[!rownames(cts) %in% conta, ]
# Sequences major proportion

#check

rownames(coldata)
colnames(cts)
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
str(cts)
str(coldata)


coldata$Genotype<-as.factor(coldata$Genotype)
coldata$Temperature<-as.factor(coldata$Temperature)
coldata$Time<-NULL


ddsLRT <- DESeqDataSetFromMatrix(countData = cts,
                                 colData = coldata,
                                 design = ~ Genotype+Temperature+Genotype*Temperature)
name_filtered <-rownames(dds)


# filtering
ddsLRT <- estimateSizeFactors(ddsLRT)
ddsLRT$sizeFactor
plot(ddsLRT$sizeFactor)
idx <- rowSums(counts(ddsLRT,normalized=TRUE) >= 7 ) >= 4 # check use normalize
ddsLRT <- ddsLRT[idx,]
dim(ddsLRT)
# filtering
ddsLRT <- estimateSizeFactors(ddsLRT)
ddsLRT <- estimateDispersions(ddsLRT)

#test against reduced model
ddsLRT <- DESeq(ddsLRT, test="LRT", reduced=~Genotype + Temperature)
resGEIT60 <- results(ddsLRT)
resGEIT60<-subset(resGEIT60,padj<0.05)
write.table(resGEIT60,file="contrasts/T60_GEI.txt",quote=F)


#plot outliers
#Counts outliers
W <- res$stat
maxCooks <- apply(assays(dds)[["cooks"]],1,max)
idx <- !is.na(W)
plot(rank(W[idx]), maxCooks[idx], xlab="rank of Wald statistic", 
     ylab="maximum Cook's distance per gene",
     ylim=c(0,5), cex=.4, col=rgb(0,0,0,.3))
m <- ncol(dds)
p <- 3
abline(h=qf(.99, p, m - p))

### Plot specific genes
coldata<-read.table("design.txt", header=T,check.names=F)
cts<-read.table("join_gamma.txt", header=T,check.names=F)
rownames(cts)<-cts$genes
cts$genes<-NULL

cts<-cts[,colnames(cts) %in% rownames(coldata) ]
cts <- cts[,order(colnames(cts))] 

#clean coldata
vectorname<-colnames(cts)
coldata<-coldata[rownames(coldata) %in% vectorname, ]
conta=as.character(c("TRINITY_DN63165_c1_g2_i6", "TRINITY_DN88404_c13_g1_i8"))
cts<-cts[!rownames(cts) %in% conta, ]

rownames(coldata)<-coldata$Ind
cts_t<-t(cts)
coldata<-coldata[coldata$Time == 1,]
data.merge<-merge(coldata,cts_t,by=0)

rownames(data.merge)<-data.merge$Row.names
data.merge$Row.names<-NULL

data.merge$Time<-NULL
data.merge$Genotype<-as.factor(data.merge$Genotype)
data.merge$Temperature<-as.factor(data.merge$Temperature)
summary(data.merge[,1:4])

ggplot(data.merge,aes(x=Temperature,y=TRINITY_DN23289_c0_g1_i1,color=Genotype))+geom_point()+
  facet_grid(~Genotype)
dev.off()




####### Test full model
# Loading data
coldata<-read.table("design.txt", header=T,check.names=F)
coldata<-coldata[coldata$Time != 0,]

rownames(coldata)<-coldata$Ind
coldata$Ind<-NULL
coldata$Temperature_Time<-NULL
coldata <- coldata[order(rownames(coldata)),] 

#counts
cts<-read.table("join_gamma.txt", header=T,check.names=F)
rownames(cts)<-cts$genes
cts$genes<-NULL

cts<-cts[,colnames(cts) %in% rownames(coldata) ]
cts <- cts[,order(colnames(cts))] 

#clean coldata
vectorname<-colnames(cts)
coldata<-coldata[rownames(coldata) %in% vectorname, ]
conta=as.character(c("TRINITY_DN63165_c1_g2_i6", "TRINITY_DN88404_c13_g1_i8"))
cts<-cts[!rownames(cts) %in% conta, ]
# Sequences major proportion

#check

rownames(coldata)
colnames(cts)
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
str(cts)
str(coldata)


coldata$Genotype<-as.factor(coldata$Genotype)
coldata$Temperature<-as.factor(coldata$Temperature)
coldata$Time<-as.factor(coldata$Time)


ddsFull <- DESeqDataSetFromMatrix(countData = cts,
                                 colData = coldata,
                                 design = ~ Time*Genotype*Temperature)
name_filtered <-rownames(dds)


# filtering
ddsFull <- estimateSizeFactors(ddsFull)
ddsFull$sizeFactor

# filtering
ddsFull <- estimateSizeFactors(ddsFull)
ddsFull$sizeFactor
plot(ddsFull$sizeFactor)

idx <- rowSums(counts(ddsFull,normalized=TRUE) >= 7 ) >= 4 # check use normalize
ddsFull <- ddsFull[idx,]
dim(ddsFull)
# filtering
ddsFull <- estimateSizeFactors(ddsFull)
ddsFull <- estimateDispersions(ddsFull)

#testing
ddsFull <- DESeq(ddsFull)
resultsNames(ddsFull)


#Time
FullTime<-results(ddsFull, name="Time_60_vs_1")
FullTime<-subset(FullTime,padj<0.05)
str(FullTime@rownames)
write.table(FullTime,file="models_full_fc/FullTime.txt",quote=F)


#Genotype
FullGeno<-results(ddsFull, name="Genotype_Marquesas_vs_Gambier")
FullGeno<-subset(FullGeno,padj<0.05)
str(FullGeno@rownames)
write.table(FullGeno,file="models_full_fc/FullGeno.txt",quote=F)

#Temperature
FullTemp3v1<-results(ddsFull, name="Temperature_3_vs_1")
FullTemp3v1<-subset(FullTemp3v1,padj<0.05)
str(FullTemp3v1@rownames)
write.table(FullTemp3v1,file="models_full_fc/FullTemp3v1.txt",quote=F)

FullTemp2v1<-results(ddsFull, name="Temperature_2_vs_1")
FullTemp2v1<-subset(FullTemp2v1,padj<0.05)
str(FullTemp2v1@rownames)
write.table(FullTemp2v1,file="models_full_fc/FullTemp2v1.txt",quote=F)

FullTemp4v1<-results(ddsFull, name="Temperature_4_vs_1")
FullTemp4v1<-subset(FullTemp4v1,padj<0.05)
str(FullTemp4v1@rownames)
write.table(FullTemp4v1,file="models_full_fc/FullTemp4v1.txt",quote=F)

resultsNames(ddsFull)

#GEI
FullGenot2<-results(ddsFull,name="GenotypeMarquesas.Temperature2")
FullGenot2<-subset(FullGenot2,padj<0.05)
str(FullGenot2@rownames)
write.table(FullGenot2,file="models_full_fc/FullGenot2.txt",quote=F)

FullGenot3<-results(ddsFull,name="GenotypeMarquesas.Temperature3")
FullGenot3<-subset(FullGenot3,padj<0.05)
str(FullGenot2@rownames)
write.table(FullGenot3,file="models_full_fc/FullGenot3.txt",quote=F)

FullGenot4<-results(ddsFull,name="GenotypeMarquesas.Temperature4")
FullGenot4<-subset(FullGenot4,padj<0.05)
str(FullGenot4@rownames)
write.table(FullGenot4,file="models_full_fc/FullGenot4.txt",quote=F)

#
resultsNames(ddsFull)

## missing info for timxtemp

#Three-ways interaction
FullTriplt2<-results(ddsFull, name="Time60.GenotypeMarquesas.Temperature2")
FullTriplt2<-subset(FullTriplt2,padj<0.05)
str(FullTriplt2@rownames)
write.table(FullTriplt2,file="models_full_fc/FullTriplt2.txt",quote=F)

FullTriplt3<-results(ddsFull, name="Time60.GenotypeMarquesas.Temperature3")
FullTriplt3<-subset(FullTriplt3,padj<0.05)
str(FullTriplt3@rownames)
write.table(FullTriplt3,file="models_full_fc/FullTriplt3.txt",quote=F)

FullTriplt4<-results(ddsFull, name="Time60.GenotypeMarquesas.Temperature4")
FullTriplt4<-subset(FullTriplt4,padj<0.05)
str(FullTriplt4@rownames)
write.table(FullTriplt4,file="models_full_fc/FullTriplt4.txt",quote=F)

###### plot
cts_t<-t(norm.counts.host) #test
data.merge<-merge(coldata,cts_t,by=0)
rownames(data.merge)<-data.merge$Row.names
data.merge$Row.names<-NULL
data.merge$condition<-as.factor(data.merge$condition)
head(data.merge[1:10])

library(Rmisc)

list_deg_unconsistant=c("TRINITY_DN57644_c0_g1_i2",
                        "TRINITY_DN58266_c0_g1_i1",                        
                        "TRINITY_DN59968_c0_g1_i4", 
                        "TRINITY_DN83319_c0_g2_i1",
                        "TRINITY_DN81780_c1_g2_i1",
                        "condition")


list_nadh_unconsistant=c("TRINITY_DN60556_c0_g1_i7",
  "TRINITY_DN65736_c0_g1_i1",
  "TRINITY_DN68075_c8_g1_i1",
  "TRINITY_DN71288_c0_g1_i2",
  "TRINITY_DN73358_c2_g1_i1",
  "TRINITY_DN81033_c1_g3_i1",
  "TRINITY_DN88171_c1_g1_i6",
                        "condition")

theme<-theme(panel.background = element_blank(),
             panel.border=element_rect(fill=NA),
             panel.grid.major = element_line(colour="gray"),
             axis.text.x=element_text(colour="black",size=8,angle = 90, vjust = .5),
             axis.text.y=element_text(colour="black",size=8),
             axis.title.x=element_blank(),
             axis.title.y=element_blank(),
             axis.ticks=element_line(colour="black"),
             plot.title = element_text(size=6,hjust=0),
             aspect.ratio=1,
             legend.position = "none")
summary(data.merge$condition)
keepcon<-c("1_Marquesas_34","1_Marquesas_27","1_Marquesas_32","1_Gambier_34","1_Gambier_27","1_Gambier_32",
           "60_Marquesas_34","60_Marquesas_27","60_Marquesas_32","60_Gambier_34","60_Gambier_27","60_Gambier_32")
keepcon
data.merge<-data.merge[data.merge$condition %in% keepcon,]

data.merge.hsp<-data.merge[,colnames(data.merge) %in% list_deg_unconsistant]

head(data.merge.hsp)

colNamesVec=c("condition", "HSP70_TRINITY_DN57644_c0_g1_i2",
              "Hsc70c2_TRINITY_DN58266_c0_g1_i1",
              "HSP70A1_TRINITY_DN59968_c0_g1_i4",
              "Hsc71c8_TRINITY_DN81780_c1_g2_i1",
              "HSP71_TRINITY_DN83319_c0_g2_i1")
colnames(data.merge.hsp)<-colNamesVec
head(data.merge.hsp)

#HSP
pl <- vector("list", length = ncol(data.merge.hsp)-1)
for(ii in seq_along(pl)){
  .col <- colnames(data.merge.hsp)[-1][ii]
  .p<-ggplot(data.merge.hsp,aes_string(x="condition",y=.col))+
    geom_boxplot()+
    geom_dotplot(binaxis='y', stackdir='center',bin=10, aes(fill=condition)) +
    theme +ggtitle(label=.col)+
    scale_fill_manual(values=c("yellow","orange","red","yellow","orange","red","yellow","orange","red","yellow","orange","red")) +
    scale_x_discrete(labels=c("d1G27","d1G32","d1G34","d1M27","d1M32","d1M34","d48G27","d48G32","d48G34","d48M27","d48M32","d48M34"))
 
  pl[[ii]] <- .p
}
grid.arrange(grobs=pl)
dev.off()


data.merge.nadh<-data.merge[,colnames(data.merge) %in% list_nadh_unconsistant]
#NADH
pl <- vector("list", length = ncol(data.merge.nadh)-1)
for(ii in seq_along(pl)){
  .col <- colnames(data.merge.nadh)[-1][ii]
  .p<-ggplot(data.merge.nadh,aes_string(x="condition",y=.col))+
    geom_boxplot()+
    geom_dotplot(binaxis='y', stackdir='center',bin=10, aes(fill=condition)) +
    theme +ggtitle(label=.col)+
    scale_fill_manual(values=c("yellow","orange","red","yellow","orange","red","yellow","orange","red","yellow","orange","red")) +
    scale_x_discrete(labels=c("d1G27","d1G32","d1G34","d1M27","d1M32","d1M34","d48G27","d48G32","d48G34","d48M27","d48M32","d48M34"))
  
  pl[[ii]] <- .p
}
grid.arrange(grobs=pl)

## overlap SNPs & DEG
list_overlap_unconsistant<-c("TRINITY_DN86328_c1_g2_i3",
                             "TRINITY_DN68120_c0_g2_i1",
                             "condition")
data.merge.overlap<-data.merge[,colnames(data.merge) %in% list_overlap_unconsistant]

pl <- vector("list", length = ncol(data.merge.overlap)-1)
for(ii in seq_along(pl)){
  .col <- colnames(data.merge.overlap)[-1][ii]
  .p<-ggplot(data.merge.overlap,aes_string(x="condition",y=.col))+
    geom_boxplot()+
    geom_dotplot(binaxis='y', stackdir='center',bin=10, aes(fill=condition)) +
    theme +ggtitle(label=.col)+
    scale_fill_manual(values=c("yellow","orange","red","yellow","orange","red","yellow","orange","red","yellow","orange","red")) +
    scale_x_discrete(labels=c("d1G27","d1G32","d1G34","d1M27","d1M32","d1M34","d48G27","d48G32","d48G34","d48M27","d48M32","d48M34"))
  
  pl[[ii]] <- .p
}
grid.arrange(grobs=pl)

## GST

list_gst_unconsistant=c("TRINITY_DN63624_c0_g1_i8",
                        "TRINITY_DN71566_c7_g3_i1",
                        "TRINITY_DN78912_c1_g2_i2", 
                        "TRINITY_DN79397_c0_g1_i5", 
                        "TRINITY_DN83398_c3_g1_i1",
                        "condition")
data.merge.gst<-data.merge[,colnames(data.merge) %in% list_gst_unconsistant]

pl <- vector("list", length = ncol(data.merge.gst)-1)
for(ii in seq_along(pl)){
  .col <- colnames(data.merge.gst)[-1][ii]
  .p<-ggplot(data.merge.gst,aes_string(x="condition",y=.col))+
    geom_boxplot()+
    geom_dotplot(binaxis='y', stackdir='center',bin=10, aes(fill=condition)) +
    theme +ggtitle(label=.col)+
    scale_fill_manual(values=c("yellow","orange","red","yellow","orange","red","yellow","orange","red","yellow","orange","red")) +
    scale_x_discrete(labels=c("d1G27","d1G32","d1G34","d1M27","d1M32","d1M34","d48G27","d48G32","d48G34","d48M27","d48M32","d48M34"))
  
  pl[[ii]] <- .p
}
grid.arrange(grobs=pl)

#plot single genes
ggplot(data.merge,aes_string(x="condition",y="TRINITY_DN40109_c0_g1_i1"))+
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center',bin=10, aes(fill=condition)) +
  theme



###### Plot single genes HSP
#"HSP70_TRINITY_DN57644_c0_g1_i2",
#              "Hsc70c2_TRINITY_DN58266_c0_g1_i1",
#            "HSP70A1_TRINITY_DN59968_c0_g1_i4",
#             "Hsc71c8_TRINITY_DN81780_c1_g2_i1",
#             "HSP71_TRINITY_DN83319_c0_g2_i1"

hsp<-list("TRINITY_DN57644_c0_g1_i2","TRINITY_DN58266_c0_g1_i1",
"TRINITY_DN59968_c0_g1_i4",
"TRINITY_DN81780_c1_g2_i1",
"TRINITY_DN83319_c0_g2_i1","Temperature","Time","Genotype","condition")

data.merge$Temperature<-as.factor(data.merge$Temperature)
data.merge<-data.merge[,colnames(data.merge) %in% hsp ]
data.merge <- data.merge[,order(colnames(data.merge))]
data.merge<-data.merge[data.merge$Temperature!= 32,]
data.merge<-data.merge[data.merge$Temperature!= 23,]
data.merge
shape=c(21,24)
ggplot(data.merge,aes_string(x="condition",y="TRINITY_DN84730_c1_g1_i2"))+
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center',bin=10, aes(fill=Temperature)) +
  theme +ggtitle(label="HSP70")+
  scale_fill_manual(values=c("dodgerblue","tomato2")) +
  scale_x_discrete(labels=c("27°C","34°C","27°C","34°C"))+
  scale_shape_manual(values=shape,labels=c("Gambier", "Marquesas"))+
  facet_grid(. ~ Time)


##@ Plot for test
ggplot(data.merge,aes_string(x="condition",y="TRINITY_DN80915_c0_g2_i8"))+
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center',bin=10, aes(fill=Temperature)) +
  theme


#### Plot Volcano style
#Gambier T1
list.MAR.27_34_T1<-read.table("contrasts/list_T1_MAR_27vs34.color.txt")
rownames(list.MAR.27_34_T1)<-list.MAR.27_34_T1$V1
df.resLFC_1_GAM_27_34<-as.data.frame(resLFC_1_GAM_27_34)
data.vol.GAM.T1<-merge(df.resLFC_1_GAM_27_34,list.MAR.27_34_T1,by=0,all=T)
library("dplyr")
data.vol.GAM.T1 <- data.vol.GAM.T1 %>%
  mutate(padj = coalesce(padj, 1))

list.MAR.27_34_T60<-read.table("contrasts/list_T60_MAR_27vs34.color.txt")
rownames(list.MAR.27_34_T60)<-list.MAR.27_34_T60$V1
df.resLFC_60_GAM_27_34<-as.data.frame(resLFC_60_GAM_27_34)
data.vol.GAM.T60<-merge(df.resLFC_60_GAM_27_34,list.MAR.27_34_T60,by=0,all=T)



data.vol.GAM.T60 <- data.vol.GAM.T60 %>%
  mutate(padj = coalesce(padj, 1))

theme1<-theme(panel.background = element_blank(),
             panel.border=element_rect(fill=NA),
             panel.grid.major = element_blank(),
             axis.text.x=element_text(colour="black",size=10, vjust = .5),
             axis.text.y=element_text(colour="black",size=10),
             axis.title.x=element_text(colour="black",size=12),
             axis.title.y=element_text(colour="black",size=12),
             axis.ticks=element_line(colour="black"),
             plot.title = element_text(size=6,hjust=0),
             aspect.ratio=1,
             legend.position = "none")

volcano_gambier_day2_27vs34<-ggplot(data.vol.GAM.T1,aes(x=log2FoldChange,y=(-log10(padj)),fill=V2))+geom_point(pch=21)+
  scale_fill_manual(values=c("tomato2","white"),na.value = "gray100")+ylim(0.01,40)+
  geom_hline(yintercept=2, linetype="dashed", color = "black")+
  geom_vline(xintercept=-1, linetype="dashed", color = "black")+
  geom_vline(xintercept=1, linetype="dashed", color = "black")+theme1+
  ylab("-log10(P.adj)")+xlab("Log2FC")+
  annotate("text", x = -4, y = 35, label = "Gambier 27°C vs 34°C\nday 2")


volcano_gambier_day48_27vs34<-ggplot(data.vol.GAM.T60,aes(x=log2FoldChange,y=(-log10(padj)),fill=V2))+geom_point(pch=21)+
  scale_fill_manual(values=c("tomato2","white"),na.value = "gray100")+
  geom_hline(yintercept=2, linetype="dashed", color = "black")+
  geom_vline(xintercept=-1, linetype="dashed", color = "black")+
  geom_vline(xintercept=1, linetype="dashed", color = "black")+theme1+
  ylab("-log10(P.adj)")+xlab("Log2FC")+ylim(0.01,25)+
  annotate("text", x = 4, y = 18, label = "Gambier 27°C vs 34°C\nday 48")
  
#Marquesas
list.GAM.27_34_T1<-read.table("contrasts/list_T1_GAM_27vs34.color.txt")
rownames(list.GAM.27_34_T1)<-list.GAM.27_34_T1$V1
df.resLFC_1_MAR_27_34<-as.data.frame(resLFC_1_MAR_27_34)
data.vol.MAR.T1<-merge(df.resLFC_1_MAR_27_34,list.GAM.27_34_T1,by=0,all=T)


volcano_marquesas_day2_27vs34<-ggplot(data.vol.MAR.T1,aes(x=log2FoldChange,y=(-log10(padj)),fill=V2))+geom_point(pch=21)+
  scale_fill_manual(values=c("dodgerblue","white"),na.value = "gray100")+
  geom_hline(yintercept=2, linetype="dashed", color = "black")+
  geom_vline(xintercept=-1, linetype="dashed", color = "black")+
  geom_vline(xintercept=1, linetype="dashed", color = "black")+theme1+
  ylab("-log10(P.adj)")+xlab("Log2FC")+
  annotate("text", x = -3, y = 50, label = "Marquesas 27°C vs 34°C\nday 2")



list.GAM.27_34_T60<-read.table("contrasts/list_T60_GAM_27vs34.color.txt")
rownames(list.GAM.27_34_T60)<-list.GAM.27_34_T60$V1
df.resLFC_60_MAR_27_34<-as.data.frame(resLFC_60_MAR_27_34)
data.vol.MAR.T60<-merge(df.resLFC_60_MAR_27_34,list.GAM.27_34_T60,by=0,all=T)


volcano_marquesas_day60_27vs34<-ggplot(data.vol.MAR.T60,aes(x=log2FoldChange,y=(-log10(padj)),fill=V2))+geom_point(pch=21)+
  scale_fill_manual(values=c("dodgerblue","white"),na.value = "gray100")+
  geom_hline(yintercept=2, linetype="dashed", color = "black")+
  geom_vline(xintercept=-1, linetype="dashed", color = "black")+
  geom_vline(xintercept=1, linetype="dashed", color = "black")+theme1+
  ylab("-log10(P.adj)")+xlab("Log2FC")+
  annotate("text", x = 3, y = 25, label = "Marquesas 27°C vs 34°C\nday 48")


library(ggpubr)

pdf(file="volcano_contrast_27_34.pdf",height=10,width=10)
ggarrange(volcano_gambier_day2_27vs34, volcano_gambier_day48_27vs34, 
          volcano_marquesas_day2_27vs34,volcano_marquesas_day60_27vs34,
          ncol = 2, nrow = 2)
dev.off()

### plot enhanced volcano
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

head(T1_GAM_27vs34)
data.res.2d.merge<-merge(as.T1_GAM_27vs34,T1_MAR_27vs34,by=0)
## Gambier 2d 27C vs 34°C 
Volc1dGAM27Cvs34C<-EnhancedVolcano(T1_GAM_27vs34,
                lab = rownames(T1_GAM_27vs34),
                x = 'log2FoldChange',
                y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab= bquote(~Log[10]~'FDR'),
                title = 'Gambier 2d 27°C vs 34°C',
                pCutoff = 0.01,
                FCcutoff = 2,
                pointSize = 3.0,
                labSize = 2.0)

Volc60dGAM27Cvs34C<-EnhancedVolcano(T60_GAM_27vs34,
                lab = rownames(T60_GAM_27vs34),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Gambier 60d 27°C vs 34°C',
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab= bquote(~Log[10]~'FDR'),
                pCutoff = 0.01,
                FCcutoff = 2,
                pointSize = 3.0,
                labSize = 2.0)

Volc60dMAR27Cvs34C<-EnhancedVolcano(T60_MAR_27vs34,
                lab = rownames(T60_MAR_27vs34),
                x = 'log2FoldChange',
                y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab= bquote(~Log[10]~'FDR'),
                title = 'Marquesas 60d 27°C vs 34°C',
                pCutoff = 0.01,
                FCcutoff = 2,
                pointSize = 3.0,
                labSize = 2.0)

Volc1dMAR27Cvs34C<-EnhancedVolcano(T1_MAR_27vs34,
                lab = rownames(T1_MAR_27vs34),
                x = 'log2FoldChange',
                y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab= bquote(~Log[10]~'FDR'),
                title = 'Marquesas 2d 27°C vs 34°C',
                pCutoff = 0.01,
                FCcutoff = 2,
                pointSize = 3.0,
                labSize = 2.0)

grid.arrange(Volc1dMAR27Cvs34C, Volc60dMAR27Cvs34C,Volc1dGAM27Cvs34C, Volc60dGAM27Cvs34C,
             ncol = 2,nrow=2)


## correlation log2FC
data.1d.gam<-read.table("~/Desktop/T1_GAM_27vs34.txt",header=T)
head(data.1d.gam)
data.1d.mar<-read.table("~/Desktop/T1_MAR_27vs34.txt",header=T)

temp.merge<-merge(data.1d.mar,data.1d.gam,by=0)
head(temp.merge)
ggscatter(temp.merge, x = "log2FoldChange.x", y = "log2FoldChange.y", 
         conf.int = TRUE, 
         cor.coef = TRUE, cor.method = "pearson",
         ylab = "gam", xlab = "mar")
