#!/usr/bin/env bash
#PBS -N vcftools
#PBS -o 98_log_files/filtering.1st.step.vcftools.err
#PBS -l walltime=50:00:00
#PBS -l mem=5g
#PBS -r n



cd $PBS_O_WORKDIR

#load vcftools module
#vcftools0.1.16

InDir="04_mapped/gatk"
OutDir="06_outlier"

MAF="0.1"
MISS="0.9"
MINQ=30
MINDP=10

#prefilter
#vcftools --gzvcf "$InDir"/combined.raw.vcf.gz --remove-indels --min-alleles 2 --max-alleles 2 --max-missing "$MISS" --minQ "$MINQ" --minDP "$MINDP" --maf "$MAF" --recode --out 06_outlier/combined.minDP"$MINDP".MAF"$MAF".MISS"$MISS".biallelic 

#sed -i 's/HI.[0-9]*.[0-9]*.NEBNext_dual_i[0-9]*_[A-Z]*[0-9]*---NEBNext_dual_i[0-9]*_[A-Z]*[0-9]*./s/g' 06_outlier/combined.minDP"$MINDP".MAF"$MAF".MISS"$MISS".biallelic.recode.vcf

##Select only subset samples
IND="01_info_files/list_ind_outlier.txt"
vcftools --vcf "$OutDir"/raw.vcf --keep "$IND" --max-missing "$MISS" --minGQ "$MINQ" --minDP "$MINDP" --maf "$MAF" --recode --out 06_outlier/combined.minDP"$MINDP".MAF"$MAF".MISS"$MISS".Q"$MINQ".biallelic.subset

## check missing per individual

vcftools --vcf 06_outlier/combined.minDP"$MINDP".MAF"$MAF".MISS"$MISS".Q"$MINQ".biallelic.subset.recode.vcf --missing-indv
