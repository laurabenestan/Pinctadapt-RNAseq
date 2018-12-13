#!/usr/bin/env bash
#PBS -q
#PBS -l walltime=700:00:00
#PBS -l mem=115g
#PBS -l ncpus=56

WORKDIR=gamma
BAYESCAN=
INDIR=07_snps/04_bayescan
NCPU=56

#conda create -n bayescan python=2.7 anaconda
#source activate bayescan
#conda install -c bioconda bayescan
#source deactivate bayescan

mkdir -p $WORKDIR/$INDIR
cd $WORKDIR/$INDIR
source activate bayescan

bayescan2 Bayescan_input.txt -snp -threads $NCPU -out_freq -pr_odds 300 -od ./ -o GAMMA_bayescan_out

