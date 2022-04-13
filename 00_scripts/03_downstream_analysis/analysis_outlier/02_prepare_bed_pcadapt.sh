#!/bin/bash

NAME=gamma.pcadapt.clean

vcftools --vcf ../SNP_exploration/beagle.imputed.vcf --plink-tped --out $NAME

plink --tfam "$NAME".tfam --tped "$NAME".tped --make-bed --out renamed.out

#prepare genotype matrix
vcftools --vcf ../SNP_exploration/beagle.imputed.vcf --012 --out beagle.imputed.genotype
