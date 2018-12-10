#!/usr/bin/env bash
#PBS -q
#PBS -l walltime=24:00:00
#PBS -l mem=115g
#PBS -N BCFtools

WORKDIR=gamma
INDIR=07_snps/02_freebayes
BCFLIBENV=

$BCFLIBENV
cd $WORKDIR/$INDIR

# Goal: remove the multi allelic list
# Manual for bcftools : http://samtools.github.io/bcftools/bcftools.html#view
# Use -m2 -M2 -v snps to only view biallelic SNPs.

bcftools view -m2 -M2 -v snps DP10_SNP_MAF0.1_miss0.1_GAMMA_genome.recode.vcf -o DP10_SNP_MAF0.1_miss0.1_GAMMA_genome_biallelic.vcf

