#!/usr/bin/env bash
#PBS -q
#PBS -l walltime=24:00:00
#PBS -l mem=115g
#PBS -N snpeff


SNPEFF=
DBNAME=

WORKDIR=gamma
OUTDIR=07_snps/03_SNPeff
VCF=07_snps/01_freebayes/DP10_SNP_MAF0.1_miss0.1_GAMMA_genome_biallelic_IndFiltered_noComplex_imputed.vcf
mkdir -p $WORKDIR/$OUTDIR

cd $WORKDIR/$OUTDIR

java -Xmx115G -jar $SNPEFF/snpEff.jar $DBNAME $WORKDIR/$VCF > DP10_SNP_MAF0.1_miss0.1_GAMMA_genome_biallelic_IndFiltered_noComplex_imputed_annotated.vcf ;

mv $SNPEFF/snpEff_summary.html DP10_SNP_MAF0.1_miss0.1_GAMMA_genome_biallelic_IndFiltered_imputed.html ;
mv $SNPEFF/snpEff_genes.txt DP10_SNP_MAF0.1_miss0.1_GAMMA_genome_biallelic_IndFiltered_imputed_genes_summary.txt ;






