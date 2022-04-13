#!/usr/bin/Rscript

################ WGCNA RNA-seq Analysis project snapper ############
####################################################
ls()
rm(list=ls())
ls()

# Installing packages
#source('http://bioconductor.org/biocLite.R')
#biocLite('DESeq2')


# Loading packages
library('DESeq2')
library("ggplot2")
#install.packages("pheatmap")
library("pheatmap")
#install.packages("RColorBrewer")
library("RColorBrewer")
#install.packages("vsn") #not available
#library("vsn")
#source("https://bioconductor.org/biocLite.R")
#biocLite("IHW")


#Global variables
setwd("analysis_deg/")


# Loading data
coldata<-read.table("../../01_info_files/design.txt", header=T,check.names=F)
coldata$Temperature_Time<-NULL
coldata<-coldata[coldata$Time != 0,]

rownames(coldata)<-coldata$Ind
coldata$Ind<-NULL
coldata <- coldata[order(rownames(coldata)),] 
coldata <- coldata[(rownames(coldata)) != 16,] 

#counts
cts<-read.table("join_global_genome_110621.txt", header=T,check.names=F) 
rownames(cts)<-cts$genes
cts$genes<-NULL

cts<-cts[,colnames(cts) %in% rownames(coldata) ]
cts <- cts[,order(colnames(cts))] 


# Sequences major proportion
summary(colSums(cts, na.rm = FALSE, dims = 1))
colSums(cts, na.rm = FALSE, dims = 1)
mean(colSums(cts, na.rm = FALSE, dims = 1))
#check

rownames(coldata)
colnames(cts)
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))

#Create dds
coldata$Temperature<-as.factor(coldata$Temperature)

coldata$Genotype<-as.factor(coldata$Genotype)
coldata$Time<-as.factor(coldata$Time)



# buolt objects
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Temperature + Time + Genotype)

name_filtered <-rownames(dds)


# filtering
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

#For test
idx <- rowSums(counts(dds,normalized=TRUE) >= 6 ) >= 4  # check use normalize
dds <- dds[idx,]
dim(dds)

#matrix log
norm.counts <- counts(dds, normalized=TRUE)
log.norm.counts <- log2(norm.counts + 2)
#write.table(log.norm.counts,file="03_results/log_matrix_gamma_noT0.txt",quote=F) # transcriptome
#write.table(log.norm.counts,file="log2cpm_matrix_gamma_genome_110621.txt",quote=F)
df.log.norm.counts<-as.data.frame(log.norm.counts)
norm.counts<-as.data.frame(norm.counts)

# Check redundancy
majSequences <- function(counts, n=3, group, col=c("lightblue","orange","MediumVioletRed","SpringGreen"), outfile=TRUE){
  
  seqnames <- apply(counts, 2, function(x){x <- sort(x, decreasing=TRUE); names(x)[1:n]})
  seqnames <- unique(unlist(as.character(seqnames)))
  
  sum <- apply(counts,2,sum)
  counts <- counts[seqnames,]
  sum <- matrix(sum,nrow(counts),ncol(counts),byrow=TRUE)
  p <- round(100*counts/sum,digits=3)
  
  if (outfile) png(filename="est.png",width=min(3600,1800+800*ncol(counts)/10),height=1800,res=300)
  maj <- apply(p, 2, max)
  seqname <- rownames(p)[apply(p, 2, which.max)]
  x <- barplot(maj, col=col[as.integer(group)], main="Percentage of reads from most expressed sequence",
               ylim=c(0, max(maj)*1.2), las=2, ylab="Percentage of reads")
  legend("topright", levels(group), fill=col[1:nlevels(group)], bty="n")
  for (i in 1:length(seqname)) text(x[i], maj[i]/2, seqname[i], cex=0.8, srt=90, adj=0)
  if (outfile) dev.off()
  
  return(invisible(p))
}

majSequences(counts(dds), group=as.factor(coldata$Genotype))

# Data transformation

vsd.fast <- vst(dds, fitType='local',blind=FALSE)
write.table(assay(vsd.fast),file="../analysis_wgcna/vst_matrix_genome_110621.txt", quote=F)

## Correlation matrix individuals
sampleDists <- dist(t(assay(vsd.fast)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd.fast$Time, vsd.fast$Temperature, 
                                    vsd.fast$Genotype,sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()
### plot PCA

plotData<-plotPCA(vsd.fast,intgroup=c("Temperature","Time","Genotype"), returnData=TRUE)
percentVar <- round(100 * attr(plotData, "percentVar"))
PCA<-ggplot(plotData, aes(PC1, PC2,color=Temperature,shape=Genotype)) + geom_point(aes(size=Time),alpha=0.8,stroke=1)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed()
PCA + geom_text(aes(label = name))

theme_glob<-theme(axis.text.x=element_text(colour="black",size=12),
                  axis.text.y=element_text(colour="black",size=12),
                  axis.title.x=element_text(colour="black",size=12,face="bold"),
                  axis.title.y=element_text(colour="black",size=12,face="bold"),
                  panel.background = element_blank(),
                  panel.border=element_rect(colour = "black", fill=NA, size=1.5),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.margin=unit(c(1,1,1,1),"line"),
                  legend.title=element_blank(),
                  aspect.ratio=1,
                  legend.background = element_blank(),
                  legend.key = element_rect(fill=NA),
                  legend.text = element_text(colour="black",size=12,face="bold")) 
couleurs=c("dodgerblue2","orange2","black","red","pink") 
shape=c(21,24)
size=c(1,2,4)
#tiff(file = "03_results/rda.tiff", width = 20, height = 15, units = "cm", res = 300)
set.seed(1)
graph.pca<-PCA +theme_glob + scale_color_manual(values=couleurs,labels=c("23°C","27°C","32°C","34°C"))+
  scale_size_manual(values=size,labels=c("T0", "T1", "T3"))+
    scale_shape_manual(values=shape,labels=c("Gambier","Marquises"))+
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted")
graph.pca
#dev.off()


## extract loadings
data<-assay(vsd.fast)
head(data)
pca <- prcomp(t(data),center=T,scale=T)
pca
loadings.PC1<-data.frame(sort(abs(pca$rotation[,"PC1"]), decreasing=TRUE)[1:50])
write.table(rownames(loadings.PC1),file="03_results/loadings.PC1.txt",quote=F)
rownames(loadings.PC1)
loadings.PC2<-data.frame(sort(abs(pca$rotation[,"PC2"]), decreasing=TRUE)[1:50])
write.table(rownames(loadings.PC2),file="03_results/loadings.PC2.txt",quote=F)
rownames(loadings.PC1)
head(loadings.PC1)


####### Plot individual genes
#Plot specific genes
cts_t<-t(df.log.norm.counts)
datamerge<-merge(coldata,cts_t,by=0)
datamerge[1:5,1:5]
# plot
theme1<-theme(panel.background = element_blank(),
              panel.border=element_rect(colour = "black", fill=NA, size=1.5),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.background=element_blank(),
              axis.text.x=element_text(colour="black",size=20),
              axis.text.y=element_text(colour="black",size=20),
              axis.title.x=element_text(colour="black",size=20,face="bold"),
              axis.title.y=element_text(colour="black",size=20,face="bold"),
              axis.ticks=element_line(colour="black"),
              plot.margin=unit(c(1,1,1,1),"line"),
              aspect.ratio=1,
              legend.background = element_blank(),
              legend.text=element_text(size=18,face="bold"),
              legend.key = element_rect(fill=NA))

### RDA
##################################### RDA ###########################################


genet=assay(vsd.fast)

#rownames(genet)<-genet$gene
#genet$gene <- NULL
genet_trans=t(genet)
row.names(genet_trans)
str(genet_trans)
str(genet)
head(genet)
design=read.table("design.txt", sep="\t",header=T,na.strings="NA") # design_cell.culture.txt
rownames(design) <-design$Ind
design$Ind <-NULL
summary(design)

#subset remove T0
design<-design[ which(design$Time != 0), ]
rownames(design)

df<-merge(design, genet_trans, by=0)
summary(df)
rownames(df)=df$Row.names
df$Row.names <- NULL
rownames(df)

df <- df[order(rownames(df)),] 
rownames(df)
head(design)

#mettre dans l'ordre
df=df[ order(df$Temperature), ]
rownames(df)

df$Temperature<-as.factor(df$Temperature)
df$Time<-as.factor(df$Time)
df$Temperature_Time<-NULL
df$Genotype<-as.factor(df$Genotype)

#produire PCOA ?? partir d'une matrice euclidean
library(ape)
library(vegan)
library(cluster)
library(dplyr)

str(df)

genet.pcoa=pcoa(daisy(df[,4:ncol(df)], metric="euclidean"))
genet.pcoa$values

#Produire s??lection de variable de la db-RDA (var2 = sex, var3 = river, var4 = treatment)
ordistep(rda(genet.pcoa$vectors[,1:17]~df$Time+df$Temperature+df$Genotype, scale=F))

#Produire db-RDA avec mod??le s??lectionn?? par ordistep, soit sex + river
rda_genet=rda(formula = genet.pcoa$vectors[, 1:17] ~ as.numeric(df$Time)+ as.numeric(df$Temperature)+ as.numeric(df$Genotype), scale = F)

#test significance of global model
anova(rda_genet, step=1000)

#Test significance for each variable
anova(rda_genet, step=1000, by='margin')

#Calcul?? le pourcentage de variation expliqu?? par le mod??le global      
RsquareAdj(rda_genet)

#partial db-RDA for time alone
rda_genet2=rda(genet.pcoa$vectors[,1:17],as.numeric(df$Time),as.numeric(df$Temperature)*as.numeric(df$Genotype), scale=F)

#Test?? significativit?? du mod??le
anova(rda_genet2, step=1000)

#RsquareAdj(rda_genet3)
RsquareAdj(rda_genet2)


#partial db-RDA for genotype alone
rda_genet3=rda(genet.pcoa$vectors[,1:17],as.numeric(df$Genotype),as.numeric(df$Time)*as.numeric(df$Temperature), scale=F)

#Test?? significativit?? du mod??le
anova(rda_genet3, step=1000)

#RsquareAdj(rda_genet3)
RsquareAdj(rda_genet3)


#partial db-RDA for temperature
rda_genet4=rda(genet.pcoa$vectors[,1:17],as.numeric(df$Temperature),as.numeric(df$Time)*as.numeric(df$Time), scale=F)

#Test?? significativit?? du mod??le
anova(rda_genet4, step=1000)

#RsquareAdj(rda_genet3)
RsquareAdj(rda_genet4)


# GRAPHIQUES DE LA RDA;
# -------------------;
library(ggplot2)
sommaire = summary(rda_genet)
sommaire
df1  <- data.frame(sommaire$sites[,1:2])       # PC1 and PC2
df2  <- data.frame(sommaire$species[,1:2])
df2
str(df1)
# prepare df1 with group info
df1<-merge(df1, design, by=0, all=TRUE)
rownames(df1)=df1$Row.names
df1$Row.names <- NULL
colnames(df1)
df1$Time<-as.factor(df1$Time)
df1$Genotype<-as.factor(df1$Genotype)
df1$Temperature<-as.factor(df1$Temperature)

#mettre dans l'ordre
df1=df1[ order(df1$Time), ]

# prepare df2
head(df2)
str(df2)
df2subset<-df2[c("Axis.1","Axis.2","Axis.3","Axis.4","Axis.5","Axis.6","Axis.7","Axis.8"),]
str(df2subset)

# loadings for PC1 and PC2
title_lab=expression(bold(paste("Model adj. ",R^{2},"=0.45; P < 0.001")))

theme1<-theme(panel.background = element_blank(),
              panel.border=element_rect(colour = "black", fill=NA, size=1.5),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.background=element_blank(),
              axis.text.x=element_text(colour="black",size=20),
              axis.text.y=element_text(colour="black",size=20),
              axis.title.x=element_text(colour="black",size=20,face="bold"),
              axis.title.y=element_text(colour="black",size=20,face="bold"),
              axis.ticks=element_line(colour="black"),
              plot.margin=unit(c(1,1,1,1),"line"),
              legend.title=element_blank(),
              aspect.ratio=1,
              legend.background = element_blank(),
              legend.text=element_text(size=14,face="bold"),
              legend.key = element_rect(fill=NA))

couleurs=c("dodgerblue","yellow","orange","red") 
shape=c(21,24)
size=c(3,6)

graph.rda.genet <- ggplot(df1, aes(x=RDA1, y=RDA2,shape=Genotype,fill=Temperature)) + 
  geom_point(aes(size=Time),alpha=0.8,stroke=1)+
 scale_fill_manual(values=couleurs,labels=c("23°C","27°C","32°C","34°C"))+
  scale_size_manual(values=size,labels=c("Day 2", "Day 48"))+

  theme1 + 
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  labs(x="RDA1 (54.81%)",y="RDA2 (30.7%)")+ #contrained 
  geom_segment(data=df2subset, aes(x=0, xend=RDA1, y=0, yend=RDA2), 
               color="black", arrow=arrow(length=unit(0.01,"npc")),inherit.aes = FALSE)+
  annotate(geom="text",x=8,y=9,label=title_lab,size=5) +
  scale_shape_manual(values=shape,labels=c("Gambier","Marquesas"))+
  guides(fill=guide_legend(override.aes=list(shape=21,size=3)))
graph.rda.genet

