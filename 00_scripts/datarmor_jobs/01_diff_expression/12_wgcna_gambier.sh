#!/bin/bash
#PBS -l walltime=20:00:00
#PBS -l mem=200gb
#PBS -N "wgcnai_gambier"
#PBS -o "log.wgcna.gambier.out"
#PBS -l ncpus=1
#PBS -r n

cd $PBS_O_WORKDIR

Rscript --vanilla wgcna_gambier.R
