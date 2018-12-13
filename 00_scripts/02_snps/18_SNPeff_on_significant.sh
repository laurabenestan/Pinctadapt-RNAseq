#!/usr/bin/env bash
#PBS -q mpi
#PBS -l walltime=24:00:00
#PBS -l mem=115g
#PBS -N snpeff


SNPEFF=
DBNAME=

WORKDIR="gamma"
OUTDIR=/07_snps/03_SNPeff
VCF=/07_snps/01_freebayes/Bayescan_significant_vcf.vcf
mkdir -p $OUTDIR
cd $SNPEFF

java -Xmx115G -jar $SNPEFF/snpEff.jar $DBNAME $WORKDIR/$VCF > $OUTDIR/Bayescan_significant_vcf.vcf ;

mv $SNPEFF/snpEff_summary.html $OUTDIR/Bayescan_significant_vcf.html ;
mv $SNPEFF/snpEff_genes.txt $OUTDIR/Bayescan_significant_vcf_genes_summary.txt ;
