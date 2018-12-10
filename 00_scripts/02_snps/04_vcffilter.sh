#!/usr/bin/env bash
#PBS -q
#PBS -l walltime=24:00:00
#PBS -l mem=115g
#PBS -N vcffilter

WORKDIR=gamma/07_snps
INDIR=02_freebayes
VCFLIBENV=


cd $WORKDIR/$INDIR

$VCFLIBENV/vcffilter -g "DP > 10" -f "TYPE = snp" GAMMA_genome_split_noRealignment_parallel.vcf > DP10_SNP_GAMMA_genome.vcf



