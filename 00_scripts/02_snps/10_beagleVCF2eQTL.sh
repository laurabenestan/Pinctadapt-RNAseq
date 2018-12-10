#!/usr/bin/env bash
#PBS -q
#PBS -l walltime=24:00:00
#PBS -l mem=115g

WORKDIR=gamma
OUTDIR=07_snps/03_eQTL
INPUT=07_snps/02_freebayes/DP10_SNP_MAF0.1_miss0.1_GAMMA_genome_biallelic_IndFiltered_noComplex_imputed.vcf
GENO=07_snps/03_eQTL/genotype_full_samples.txt
SNP=07_snps/03_eQTL/snplocation_full_samples.txt


mkdir -p $WORKDIR/$OUTDIR

cd $WORKDIR/$OUTDIR
gunzip $WORKDIR/$INPUT.gz
grep -v "##" $WORKDIR/$INPUT | awk '{$1=$1"_"$2; print $0}' OFS="\t" | sed 's/0|0/0/g' | sed 's/0|1/1/g' | sed 's/1|0/1/g' | sed 's/1|1/2/g' | cut -f1,10- | sed 's/#CHROM_POS/id/g'> $WORKDIR/$GENO

grep -v "##" $WORKDIR/$INPUT | awk '{print $1"_"$2"\t"$1"\t"$2}' | sed 's/#CHROM_POS/snp/g' | sed 's/#CHROM/chr/g' | sed 's/POS/pos/g' > $WORKDIR/$SNP

