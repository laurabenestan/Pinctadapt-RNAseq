#!/bin/bash
#PBS -q
#PBS -l walltime=20:00:00
#PBS -l mem=60g
#PBS -l ncpus=12
#PBS -N make_bayescan_input

WORKDIR=gamma
INDIR=07_snps/04_bayescan
RCODE=00_scripts/02_snps/Bayescan_results_analysis.R

mkdir -p $WORKDIR/$INDIR

source activate vcfR

Rscript --vanilla $WORKDIR/$RCODE &> $WORKDIR/00_scripts/make_bayescan_input.R.out ;

source deactivate vcfR
