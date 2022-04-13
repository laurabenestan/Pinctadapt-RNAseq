          #####################################
          #####################################
          ## SNP analysis from RNAseq GAMMA  ##
          #####################################
          #####################################



######################################################
##  I.   clear environment and set working directory
######################################################

ls()
rm(list=ls())
ls()

#Set working directory to source file location


############################################################
##  II.   start analysis
############################################################



#Load all the packages needed:
library(grur)
library(radiator)
library(adegenet)
library(vcfR)
library(ggplot2)
library(plotly)
library(plyr)
library(hierfstat)
library(cluster)
library(zoo)
library(dartR)
library(vegan)
library(diveRsity)

setwd("~/06_outlier/")
#~50K variants à analyser
#vcf<-read.vcfR("combined.minDP10.MAF0.1.MISS0.9.Q30.biallelic.subset.recode.vcf")
vcf<-read.vcfR("beagle.imputed.vcf")
metadata<-read.table("../01_info_files/sample_metadata.txt",header=T)
row.names(metadata)<-metadata$Ind
metadata

#convert the VCF inxto a genlight object

data<-vcfR2genlight(vcf)
str(data)
data
data@ind.names

# Add population:
metadata<-metadata[row.names(metadata) %in% data$ind.names, ]
desiredOrder<-data$ind.names
metadata<-metadata[order(factor(metadata$Ind, levels=unique(desiredOrder))),]


data$pop<-as.factor(metadata$Genotype)
data$pop
data$other$Temperature<-as.factor(metadata$Temperature)
data$other$Days<-as.factor(metadata$Days)
data$other

strata<-metadata[,c(1,3)]
strata$STRATA<-strata$Genotype
strata$Genotype<-NULL
strata$INDIVIDUALS<-strata$Ind
strata$Ind<-NULL

## Vusualize missing data
#tidy.data<-data<-vcfR2tidy(vcf)

#library(grur)
#library(radiator)
#detect_genomic_format(tidy.data)
#missing_visualization(
 # data = tidy.data,
 # strata =NULL)

#ibm$ibm.plots$ibm.strata.POP_ID
########################################################################
## Differentiation and candidate SNPs with Fst and multivariate approach

######################################################
#Perform a principal component analysis (classic PCA):
# no supposition of any grouping
set.seed(1)
pca.1 <- glPca(data) 
3 # keep three axes
pca.1
# proportion of explained variance by first three axes
pca.1$eig[1]/sum(pca.1$eig) # proportion of variation explained by 1st axis
pca.1$eig[2]/sum(pca.1$eig) # proportion of variation explained by 2nd axis 
pca.1$eig[3]/sum(pca.1$eig) # proportion of variation explained by 3rd axis 

Genotype<-data@pop
Days<-data$other$Days
Temperature<-data$other$Temperature
gr1<-cbind(as.data.frame(pca.1$scores),Genotype,Days,Temperature)
gr1$Genotype<-revalue(gr1$Genotype, c("1"="Gambier", "2"="Marquesas"))

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
shape=c(21,24)
PCA_pop<-ggplot(gr1,aes(x=gr1$PC1,y=gr1$PC2,shape=Genotype))+ geom_point(size=4,stroke = 1,color="black",fill="grey90") + 
  scale_shape_manual(values=shape) +theme1 + xlab("PC1 (25.95%)")+ylab("PC2 (5.01%)")+
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  annotate(geom="text",x=-20,y=30,label="Gam. - Mar.\nFst = 0.21; P < 0.001",fontface =2,size=5,colour = "black")
PCA_pop
PCA_pop

# SNPR Relate -----
library(SNPRelate)
snpgdsVCF2GDS(vcf.fn="06_outlier/beagle.imputed.vcf",out.fn="test1.gds",method="biallelic.only")
snpgdsSummary("test1.gds")
f<-snpgdsOpen("test1.gds")
x2<-2-snpgdsGetGeno(f) 

data<-vcfR2genind(vcf)
pop.map<-read.table("../pcadapt_outflank2/pop.map.txt",header=F)
rownames(pop.map)<-pop.map$V1
ind.list.order.bed<-read.table("../pcadapt_outflank2/renamed.out.fam",header=F)
rownames(ind.list.order.bed)<-ind.list.order.bed$V1
target<-ind.list.order.bed$V1
merged.info.pop<-merge(ind.list.order.bed,pop.map,by=0)
#order list pop according tovector fam
library(dplyr)
merged.info.pop$names<-merged.info.pop$V1.x
ordered.merged<- merged.info.pop %>% arrange(factor(names, levels = target))
poplist.names <-ordered.merged$V2.y
print(poplist.names)
data@pop <- as.factor(poplist.names)

hierf<-genind2hierfstat(data,pop=as.factor(poplist.names))


#Using hierfstat for pop stats:
pairwise.neifst(hierf,diploid=TRUE)
#       1         2
#1      NA        0.2116
#2      0.2116    NA


## Fst andP values using Stammp
library(StAMPP)
library(vcfR)
library(adegenet)

vcf<-read.vcfR("beagle.imputed.vcf", verbose = FALSE)
data<-vcfR2genlight(vcf)
data@ind.names

# Add population:
metadata<-read.table("01_info_files/sample_metadata.txt",header=T)
row.names(metadata)<-metadata$Ind

# Add population:
metadata<-metadata[row.names(metadata) %in% data$ind.names, ]
desiredOrder<-data$ind.names
metadata<-metadata[order(factor(metadata$Ind, levels=unique(desiredOrder))),]
pop.names<-as.character(metadata$Genotype)


x2 <- as.matrix(data) #convert genlight object to matrix 
sample <- row.names(x2) #sample names 
ploidy <- ploidy(data) #extract ploidy info from genlight object 
x2 = x2 * (1/ploidy) #convert allele counts to frequency 
x2[is.na(x2)] = NaN 
format <- vector(length = length(sample))
#format id for the genotype data
format[1:length(format)] = "freq"  

x.stampp <- as.data.frame(cbind(sample, pop.names, ploidy, format, x2)) #convert to basic r data.frame suitable to stamppConvert 
geno <- stamppConvert(x.stampp, 'r') 
FST.stampp=stamppFst(geno, nboots = 1000, percent = 95, nclusters = 2)
FST <- FST.stampp$Fsts %>% as.data.frame(.)
FST[FST < 0 ] <- 0
print(FST)
Pval <- FST.stampp$Pvalues %>% as.data.frame(.)
print(Pval)

# nuleotide diversity ----
df_pi<-read.table("nucl.div.pop.txt",header=T)
couleurs_temp=c("tomato2","dodgerblue")
pi_plot<-ggplot(df_pi,aes(x=Population,y=PI,color=Population)) +
  geom_boxplot() + theme1 + ylab(expression("Nucleotide diversity"~pi))+xlab("")+
  scale_color_manual(values=couleurs_temp)+
  annotate(geom="text",x=1,y=0.6,label="Wilcox. P < 0.001",size=6,colour = "black")
pi_plot

library(Rmisc)
dftemp=summarySE(df_pi, measurevar="PI", groupvars="Population")
dftemp
#Population     N        PI        sd          se          ci
#1    Gambier 55546 0.3248005 0.1662315 0.000705321 0.001382434
#2  Marquesas 55546 0.3188667 0.1670771 0.000708909 0.001389466
