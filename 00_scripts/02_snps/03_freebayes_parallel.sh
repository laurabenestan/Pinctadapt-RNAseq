#!/bin/bash
#PBS -q 
#PBS -l walltime=100:00:00
#PBS -l ncpus=50
#PBS -lmem=115g
#PBS -N freebayes_parallel

WORKDIR=gamma
DATADIRECTORY=04_mapped
OUTDIR=07_snps/01_freebayes

FREEBAYESENV=
REF=$WORKDIR/01_info_files/***.fasta
INDEX=$WORKDIR/01_info_files/***.fasta.fai

#$FREEBAYESENV

mkdir -p $WORKDIR/$OUTDIR

#List all the BAM files contained in data folder
LS="ls $WORKDIR/$DATADIRECTORY/*split.bam"
$LS > $WORKDIR/00_scripts/02_snps/bam_list_freebayes.txt

#One individual with very little amount of reads mapped, need to be removed:
grep -v "A9.17" $WORKDIR/00_scripts/02_snps/bam_list_freebayes.txt > $WORKDIR/00_scripts/02_snps/bam_list_final_freebayes.txt
BAM=$WORKDIR/00_scripts/02_snps/bam_list_final_freebayes.txt

NCPU=50
nAlleles=4
minMapQ=30
minCOV=10

cd $WORKDIR/$OUTDIR

$FREEBAYESENV/freebayes-parallel <(fasta_generate_regions.py $INDEX 100000) "$NCPU" \
-f $REF --use-best-n-alleles $nAlleles -C 5 --min-mapping-quality $minMapQ --min-coverage $minCOV --genotype-qualities -L $BAM > GAMMA_genome_split_noRealignment_parallel.vcf


