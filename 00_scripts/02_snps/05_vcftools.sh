#!/usr/bin/env bash
#PBS -q
#PBS -l walltime=24:00:00
#PBS -l mem=115g
#PBS -N vcftools

WORKDIR=gamma
INDIR=07_snps/02_freebayes
VCFTOOLSENV=

cd $WORKDIR/$INDIR

$VCFTOOLSENV/vcftools --maf 0.1 --max-missing 0.9 --vcf DP10_SNP_GAMMA_genome.vcf --recode --out DP10_SNP_MAF0.1_miss0.1_GAMMA_genome

