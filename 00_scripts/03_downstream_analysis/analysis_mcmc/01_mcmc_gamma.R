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
library("MCMCglmm")
library(flashClust)
library(adegenet)
#Global variables
setwd="analysis_mcmc/"



# Loading data
coldata<-read.table("../../01_info_files/design.txt", header=T,check.names=F)
coldata$Temperature_Time<-NULL

# keep only 27C and 34C at T0 + short term
coldata<-coldata[coldata$Temperature != 32,]
coldata<-coldata[coldata$Temperature != 23,]
coldata<-coldata[coldata$Time != 60,]
coldata<-coldata[coldata$Time != 0,]

rownames(coldata)<-coldata$Ind
coldata$Ind<-NULL
coldata <- coldata[order(rownames(coldata)),] 
coldata <- coldata[(rownames(coldata)) != 16,] 
cols <- c("Genotype", "Time", "Temperature")
coldata$group <- do.call(paste, c(coldata[cols], sep="_"))
# Load counts ----

cts<-read.table("join_global_genome_110621.txt", header=T,check.names=F) 
rownames(cts)<-cts$genes
cts$genes<-NULL

cts<-cts[,colnames(cts) %in% rownames(coldata) ]
cts <- cts[,order(colnames(cts))] 

colnames(cts)
rownames(coldata)

# Sequences major proportion
summary(colSums(cts, na.rm = FALSE, dims = 1))
colSums(cts, na.rm = FALSE, dims = 1)
mean(colSums(cts, na.rm = FALSE, dims = 1))

rownames(coldata)<-coldata$Name
colnames(cts)<-coldata$Name
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))

#Create dds
coldata$Temperature<-as.factor(coldata$Temperature)

coldata$Genotype<-as.factor(coldata$Genotype)
coldata$Time<-as.factor(coldata$Time)



# Built dds objects ----
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Temperature * Genotype)

name_filtered <-rownames(dds)


# filtering
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

#For test
idx <- rowSums(counts(dds,normalized=TRUE) >= 6 ) >= 4  # check use normalize
dds <- dds[idx,]
dim(dds)

# Data transformation -----
vsd.fast <- vst(dds, fitType='local',blind=FALSE)
write.table(assay(vsd.fast),file="vst_matrix_mcmc_110921.txt", quote=F)

rlogCOUNTS<-rlog(dds,blind=TRUE) #use blind=TRUE to not account for experimental design
dat=as.data.frame(assay(rlogCOUNTS)) 
head(dat)
boxplot(dat) #how do expression plots look overall?

# DAPC for genes expression -----
## extract from Kenkel and Matz
head(dat) 
degs10<-rownames(dat)

a.vsd<-dat[,grep("T27",colnames(dat))] #control ind
a.vsd.supp<-dat[,grep("T34",colnames(dat))] # challenged ind


#######################Some genes have insufficient variance for DFA in data subsets!...find those genes using loop below
dframe=(a.vsd[degs10,])
for(col in rownames(dframe))
{  min=min(dframe[col,])
max=max(dframe[col,])
if(min == max){print(col)}
}
col
#If print out above, copy gene names below and run line to remove genes from analysis. 
degs9=degs10[! degs10 %in% c("evm.TU.scaffold9size467357.9")]


# DAPC ------

pcp=prcomp(t(a.vsd[degs9,]), retx=TRUE, center=TRUE, scale.=TRUE) #scale.=TRUE
scores=pcp$x
screeplot(pcp,bstick=T) # only 1st PC is non-random when using diff exp genes only; higher PCs non-random when using all data...

# adegenet: finding clusters (even though we know what clusters we want) - choose 4 PCs and 2 groups
clus=find.clusters(t(a.vsd[degs9,]),max.n.clus=5) #[degs10,]

#Use clus$grp to rename to in2in and off2off -
clus$grp
clus$grp=c("Mar","Mar","Mar","Mar","Gam","Gam","Gam","Gam") #tell the DF which groups you want to cluster; in this case in2in and off2off


# now lets build a discriminant function for these two groups:
dp=dapc(t(a.vsd[degs9,]),clus$grp) #[degs10,]
# HOST: PCs: 6, functions: 1. For two groups only one discriminant function is possible.


quartz()
scatter(dp,bg="white",scree.da=FALSE,legend=TRUE,solid=.4)+ xlim(-10,10) #discriminant function for ORIGIN type expression

#Now, add in transplants and see where they fall along continuum

# On challenged individuals
pred.sup<-predict.dapc(dp,newdata=(t(a.vsd.supp[degs9,]))) #skip IO11C for host b/c outlier sample in WGCNA


names(pred.sup)
pred.sup$assign
names(a.vsd.supp)

#must create another dataframe structure in order to plot these predicted values
test<-dp
test$ind.coord<-pred.sup$ind.scores
test$posterior<-pred.sup$posterior
test$assign<-pred.sup$assign

test$grp
#HOST - for plotting distributions, must say which samples in which group
#test$grp<-as.factor(c(2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1)) #make sure origin is same num as before: IN=2, OFF=1



quartz()

scatter(test,bg="white",scree.da=FALSE,legend=TRUE,solid=.4,xlim = c(-10, 10))

#adjust axes to correspond to previous graph, then overlay plots in photoshop

##retain DFA values for additional calculations
dpc=data.frame(rbind(dp$ind.coord,pred.sup$ind.scores))

dpc


# plot density plot
data.merge<-merge(dpc,coldata,by=0)
data.merge$Temperature<-as.factor(data.merge$Temperature)
data.merge$Genotype<-as.factor(data.merge$Genotype)

data.merge

# plot density

scale_custom <- list(
    scale_color_manual(values = c("tomato2","tomato2","dodgerblue","dodgerblue")),
    scale_fill_manual(values = c("tomato2","grey90","dodgerblue","grey90"))
)

plotdens<-ggplot(data.merge, aes(x=LD1,group=group,fill=group)) + 
  geom_density(adjust=1, alpha=.6,aes(color=group))+
  scale_x_continuous(limits = c(-11, 11))+
  theme_bw()+ 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.y = element_text(colour="black",size=12,face="plain"), 
        axis.text.x = element_text(colour="black",size=12,face="plain"))+
  scale_custom

# compute mean
aggregate(. ~ data.merge$group, data.merge[2], mean)

# plot shift


theme<-theme(axis.text.x=element_text(colour="black",size=18),
             axis.text.y=element_text(colour="black",size=18),
             axis.title.x=element_text(colour="black",size=20,face="bold"),
             axis.title.y=element_text(colour="black",size=20,face="bold"),
             panel.background = element_blank(),
             panel.border=element_rect(colour = "black", fill=NA, size=1),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             plot.margin=unit(c(1,1,1,1),"line"),
             legend.title=element_blank(),
             aspect.ratio=1,
             legend.position = "none") 
# ad significance
#title_lab=expression(bold(paste("P MCMC = 0.01")))
a=0.01
title_lab= expression(bold(paste(italic('P')['MCMC']*' = 0.01')))

plotfin<- plotdens +theme +
  geom_segment(aes(x = 6.976, y = 0.6, xend = 4.1899, yend = 0.6),
               arrow = arrow(length = unit(0.4, "cm"),type="closed"),size=1,colour = "dodgerblue") +
  geom_segment(aes(x = -6.976, y = 0.95, xend = -1.43, yend = 0.95),
               arrow = arrow(length = unit(0.4, "cm"),type="closed"),size=1,colour = "tomato2") +
  geom_segment(aes(x = -11, y = 0, xend = 11, yend = 0),size=0.5,colour = "black") +
  annotate(geom="text",x=6,y=1.5,label=title_lab,size=6)+ labs(y = "Density")
 
plotfin 

###Testing significance of DFA differences - MCMCglmm ------
library(MCMCglmm)

# weak inverse wishart prior with parameter expansion for random effect of genotype (this is standard in MCMCglmm, the results are actually identical with the default uniform prior)
prior = list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V=1, nu=0.002,alpha.mu=0, alpha.V=1000)))

cd=MCMCglmm(LD1~Genotype+Genotype:Temperature,data=data.merge,nitt=75000, thin=25, burnin=5000)
summary(cd)
#Iterations = 5001:74976
#Thinning interval  = 25
#Sample size  = 2800 
#DIC: 51.71572 

#R-structure:  ~units

#post.mean l-95% CI u-95% CI eff.samp
#units     1.231    0.427    2.385     2800

#Location effects: LD1 ~ Genotype + Genotype:Temperature 

#post.mean l-95% CI u-95% CI eff.samp   pMCMC    
#(Intercept)                         6.994    5.949    8.173     2800 < 4e-04 ***
#  GenotypeMarquesas                 -13.972  -15.442  -12.295     2800 < 4e-04 ***
#  GenotypeGambier:Temperature34      -5.569   -7.164   -4.023     2800 < 4e-04 ***
#  GenotypeMarquesas:Temperature34     2.810    1.243    4.276     2840 0.00286 ** 

# calculating difference in magnitudes of orii:away and orio:away using sampled sets of parameters:
awayDelta=abs(cd$Sol[,"GenotypeGambier:Temperature34"])-abs(cd$Sol[,"GenotypeMarquesas:Temperature34"])
mean(abs(cd$Sol[,"GenotypeGambier:Temperature34"]))
mean(abs(cd$Sol[,"GenotypeMarquesas:Temperature34"]))
mean(abs(cd$Sol[,"GenotypeMarquesas"]))
# 95% credible interval:
HPDinterval(awayDelta)
#lower    upper
#var1 0.7784509 5.094481
#attr(,"Probability")
#[1] 0.95

#MCMC p-value:
if (is.na(table(awayDelta<0)[2])) {
  cat("p <",signif(1/length(awayDelta),1))
} else { cat("p =",signif(table(awayDelta<0)[2]/length(awayDelta),2)) } # 3.6e-4

