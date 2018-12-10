#!/usr/bin/env bash
#PBS -q
#PBS -l walltime=24:00:00
#PBS -l mem=115g
#PBS -N vcftools

WORKDIR=gamma
INDIR=07_snps/01_freebayes
VCFTOOLSENV=


cd $WORKDIR/$INDIR

#Get a file "out.miss" stating the percentage of missing data per individuals in the last column
$VCFTOOLSENV/vcftools --vcf DP10_SNP_MAF0.1_miss0.1_GAMMA_genome_biallelic.vcf --missing-indv
