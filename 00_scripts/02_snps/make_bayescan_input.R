###############################
## R code for bayescan input ##
###############################
#CMO. Reisser Nov. 2018


setwd("")
#Get VCF into R
library(vcfR)
library(plyr)

#Read VCF
vcf<-read.vcfR("GAMMA_noComplex_Gambier.recode.vcf")
vcf
#get allelic frequencies
test2<-gt.to.popsum(vcf)

#Create a new variable to store number of alleles at each site for the population
test2$nAlleles<-rep("NA",length(rownames(test2)))
i<-0
for(i in 1:length(rownames(test2))){
  if(test2$Allele_counts[i]==as.character(test2$n[1]*2)){
    test2$nAlleles[i]=1}
  else if(test2$Allele_counts[i]==as.character(paste(0,",",test2$n[1]*2,sep=""))){
    test2$nAlleles[i]=1}
  else
    test2$nAlleles[i]=2
}
test2

#Create a new variable to transform the counts of the homozygous loci for reference allele
i<-0
test2$Alleles_counts2<-rep("NA",length(rownames(test2)))
for(i in 1:length(rownames(test2))){
  if(test2$Allele_counts[i]==as.character(test2$n[1]*2)){
    test2$Allele_counts2[i]=paste(test2$n[1]*2,",",0,sep="")}
  else
    test2$Allele_counts2[i]=as.character(test2$Allele_counts[i])
}

#Replace comma by space in allelic count column
test2$Allele_counts2<-sapply(test2$Allele_counts2, gsub, pattern=",", replacement=" ")

#Create object with all the data necessary
pop1<-cbind(rownames(test2),test2$n*2,test2$nAlleles,test2$Allele_counts2)
write.table(pop1,"Gambier_bayescan_input.txt",sep="\t",quote=F,row.names=F,col.names=F)


# Do the same for pop2:
vcf<-read.vcfR("GAMMA_noComplex_marquesas.recode.vcf")
vcf
test<-gt.to.popsum(vcf)
test
nn<-test$n[1]*2
test$nAlleles<-rep("NA",length(rownames(test)))
i<-0
for(i in 1:length(rownames(test))){
  if(test$Allele_counts[i]==as.character(test$n[1]*2)){
    test$nAlleles[i]=1}
  else if(test$Allele_counts[i]==as.character(paste(0,",",test$n[1]*2,sep=""))){
    test$nAlleles[i]=1}
  else
    test$nAlleles[i]=2
}

i<-0
test$Alleles_counts2<-rep("NA",length(rownames(test)))
for(i in 1:length(rownames(test))){
  if(test$Allele_counts[i]==as.character(test$n[1]*2)){
    test$Allele_counts2[i]=paste(test$n[1]*2,",",0,sep="")}
  else
    test$Allele_counts2[i]=as.character(test$Allele_counts[i])
}

test$Allele_counts2<-sapply(test$Allele_counts2, gsub, pattern=",", replacement=" ")

pop2<-cbind(rownames(test),test$n*2,test$nAlleles,test$Allele_counts2)
write.table(pop2,"Marquesas_bayescan_input.txt",sep="\t",quote=F,row.names=F,col.names=F)

