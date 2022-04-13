#!/usr/bin/env bash
#PBS -N vcftools
#PBS -o 98_log_files/missing.vcftools.err
#PBS -l walltime=00:10:00
#PBS -l mem=5g
#PBS -r n



cd $PBS_O_WORKDIR

#load vcftools module
#vcftools0.1.16

OutDir="06_outlier"
IN="beagle.imputed.vcf"

## check missing per individual

vcftools --vcf "$OutDir"/"$IN" --missing-indv --out out.imputed
