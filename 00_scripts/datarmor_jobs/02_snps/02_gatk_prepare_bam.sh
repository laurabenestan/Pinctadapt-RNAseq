#!/bin/bash


# Clean past jobs

rm 00_scripts/datarmor_jobs/GATK_*sh


# launch scripts for Colosse
for file in $(ls scracth/gamma/04_mapped/genome/*.concordant_uniq|perl -pe 's/.concordant_uniq//'|sort -u)
do

base=$(basename "$file")

	toEval="cat 00_scripts/02_snps/02_gatk_prepare_bam.sh | sed 's/__BASE__/$base/g'"; eval $toEval > 00_scripts/datarmor_jobs/GATK_$base.sh
done

#Submit jobs
for i in $(ls 00_scripts/datarmor_jobs/GATK_*sh); do qsub $i; done


