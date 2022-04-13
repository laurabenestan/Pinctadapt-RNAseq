## Pcadapt
ls()
rm(list=ls())
ls()

library(vcfR)
library(pcadapt)
setwd("06_outlier/")
#path_data<-"renamed.out.bed"
#filename <- read.pcadapt(path_data, type = "bed")

# Now import genotype 012
#choosing K
sel<-read.table("beagle.imputed.genotype.012",header=F)
dim(sel)
genotype <- sel[, 2:ncol(sel)]
dim(genotype)
pca_genotype <- read.pcadapt(t(genotype))
K <- 10
x <- pcadapt(pca_genotype, K = K)


#Add pop info
pop.map<-read.table("../01_info_files/pop.map.txt",header=F)
rownames(pop.map)<-pop.map$V1
ind.list.order.bed<-read.table("renamed.out.fam",header=F)
rownames(ind.list.order.bed)<-ind.list.order.bed$V1
target<-ind.list.order.bed$V1
merged.info.pop<-merge(ind.list.order.bed,pop.map,by=0)
#order list pop according tovector fam
library(dplyr)
merged.info.pop$names<-merged.info.pop$V1.x
ordered.merged<- merged.info.pop %>% arrange(factor(names, levels = target))

poplist.names <-ordered.merged$V2.y
print(poplist.names)


#stats
head(pca_genotype)
x<- pcadapt(pca_genotype, K = 2)
summary(x)
plot(x , option = "manhattan") + theme_bw()
plot(x, option = "qqplot")
plot(x, option = "scores",pop=poplist.names)+ theme_bw()
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
plot(x, option = "stat.distribution")

test<-prcomp(pca_genotype)

autoplot(x)

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
              legend.position="none")

df_scores<-x$scores
df_scores[1:10,]
ggplot(df_scores,aes(x=df_scores[,1]))
PCA<-plot(x, option = "scores",pop=poplist.names,pch=21,cex=6,col=c("tomato2","dodgerblue",))+theme1
PCA
#choosing a cutoff
library(qvalue)

qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.05
outliers <- which(qval < alpha)
length(outliers)
write.table(outliers,file="list_number_outliers_qval0.05.txt",col.names=F,row.names=F,quote=F)

padj <- p.adjust(x$pvalues,method="bonferroni")
alpha <- 0.1
outliers <- which(padj < alpha)
length(outliers)
x$zscores
head(outliers)
write.table(outliers,file="list_number_outliers.txt",col.names=F,row.names=F,quote=F)

head(x$scores)
