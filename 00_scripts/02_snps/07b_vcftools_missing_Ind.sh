#!/usr/bin/env bash
#PBS -q
#PBS -l walltime=24:00:00
#PBS -l mem=115g
#PBS -N vcftools

WORKDIR=gamma
INDIR=07_snps/01_freebayes
VCFTOOLSENV=


cd $WORKDIR/$INDIR

#Remove individuals taht had more than 15% missing data.
$VCFTOOLSENV/vcftools --remove-indv HI.4880.007.NEBNext_dual_i7_C1---NEBNext_dual_i5_C1.61 --remove-indv HI.4880.007.NEBNext_dual_i7_G1---NEBNext_dual_i5_G1.9 --remove-indv HI.4880.007.NEBNext_dual_i7_H12---NEBNext_dual_i5_H12.25 --remove-indv HI.4880.007.NEBNext_dual_i7_D12---NEBNext_dual_i5_D12.44 --vcf $INDIR/DP10_SNP_MAF0.1_miss0.1_GAMMA_genome_biallelic.vcf --recode --out DP10_SNP_MAF0.1_miss0.1_GAMMA_genome_biallelic_IndFiltered

