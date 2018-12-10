###############################
## R code for bayescan input ##
###############################
#CMO.Reisser

#Get VCF into R
library(vcfR)
library(plyr)

setwd("gamma/07_snps/05_bayescan")
vcf<-read.vcfR("gamma/07_snps/02_freebayes/GAMMA_noComplex_Gambier.recode.vcf")

#Get allelic summary:
test2<-gt.to.popsum(vcf)

#Create an additional colimn to get number of alleles
test2$nAlleles<-rep("NA",length(rownames(test2)))

#Add the number of alleles present in the population
i<-0
for(i in 1:length(rownames(test2))){
if(test2$Allele_counts[i]=="60"){
test2$nAlleles[i]=1}
else if(test2$Allele_counts[i]=="0,60"){
test2$nAlleles[i]=1}
else
test2$nAlleles[i]=2
}

#replace the comma separation with a space
test2$Allele_counts<-sapply(test2$Allele_counts, gsub, pattern=",", replacement=" ") 

#Add the missing 0 for reference homozygous variants
test2$Allele_counts<-revalue(test2$Allele_counts, c("60"=as.character("60 0")))

#Create the bayescan output format file for pop1:
pop1<-cbind(rownames(test2),rep(60,length(rownames(test2))),test2$nAlleles,test2$Allele_counts)
write.table(pop1,"Gambier_bayescan_input.txt",sep="\t",quote=F,row.names=F,col.names=F)


# Do the same for pop2:
vcf<-read.vcfR("gamma/07_snps/02_freebayes/GAMMA_noComplex_marquesas.recode.vcf")
vcf
test<-gt.to.popsum(vcf)
test$nAlleles<-rep("NA",length(rownames(test)))

i<-0
for(i in 1:length(rownames(test))){
  if(test$Allele_counts[i]=="72"){
    test$nAlleles[i]=1}
  else if(test$Allele_counts[i]=="0,72"){
    test$nAlleles[i]=1}
  else
    test$nAlleles[i]=2
}

test$Allele_counts<-sapply(test$Allele_counts, gsub, pattern=",", replacement=" ") 
test$Allele_counts<-revalue(test$Allele_counts, c("72"=as.character("72 0")))

pop2<-cbind(rownames(test),rep(72,length(rownames(test))),test$nAlleles,test$Allele_counts)
write.table(pop2,"Marquesas_bayescan_input.txt",sep="\t",quote=F,row.names=F,col.names=F)

