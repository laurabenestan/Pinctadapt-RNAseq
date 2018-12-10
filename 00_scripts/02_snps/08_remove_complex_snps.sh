#!/usr/bin/env bash
#PBS -q
#PBS -l walltime=24:00:00
#PBS -l mem=115g
#PBS -N vcftools

WORKDIR=gamma
INDIR=07_snps/02_freebayes

cd $WORKDIR/$INDIR

grep "^#" DP10_SNP_MAF0.1_miss0.1_GAMMA_genome_biallelic_IndFiltered.recode.vcf > header.txt
awk '$4=="A"' DP10_SNP_MAF0.1_miss0.1_GAMMA_genome_biallelic_IndFiltered.recode.vcf > A.txt
awk '$4=="C"' DP10_SNP_MAF0.1_miss0.1_GAMMA_genome_biallelic_IndFiltered.recode.vcf > C.txt
awk '$4=="G"' DP10_SNP_MAF0.1_miss0.1_GAMMA_genome_biallelic_IndFiltered.recode.vcf > G.txt
awk '$4=="T"' DP10_SNP_MAF0.1_miss0.1_GAMMA_genome_biallelic_IndFiltered.recode.vcf > T.txt

cat A.txt C.txt G.txt T.txt > temp_DP10_SNP_MAF0.1_miss0.1_GAMMA_genome_biallelic_IndFiltered.recode_noComplex.vcf
sort -k2 -n temp_DP10_SNP_MAF0.1_miss0.1_GAMMA_genome_biallelic_IndFiltered.recode_noComplex.vcf | sort -k1 > sorted_temp_DP10_SNP_MAF0.1_miss0.1_GAMMA_genome_biallelic_IndFiltered.recode_noComplex.vcf
cat header.txt sorted_temp_DP10_SNP_MAF0.1_miss0.1_GAMMA_genome_biallelic_IndFiltered.recode_noComplex.vcf > DP10_SNP_MAF0.1_miss0.1_GAMMA_genome_biallelic_IndFiltered.recode_noComplex.vcf

rm temp_DP10_SNP_MAF0.1_miss0.1_GAMMA_genome_biallelic_IndFiltered.recode_noComplex.vcf
rm sorted_temp_DP10_SNP_MAF0.1_miss0.1_GAMMA_genome_biallelic_IndFiltered.recode_noComplex.vcf
rm header.txt
rm A.txt
rm C.txt
rm G.txt
rm T.txt




