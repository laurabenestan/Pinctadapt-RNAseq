#!/bin/bash


# Clean past jobs

rm 00_scripts/datarmor_jobs/GATK_*sh


# launch scripts for Colosse
for file in $(ls 04_mapped/transcriptome/subset/*sorted.bam|sed -e 's/.concordant_uniq.sorted.bam//'|sort -u)
do

base=$(basename "$file")

	toEval="cat 00_scripts/08_gatk.sh | sed 's/__BASE__/$base/g'"; eval $toEval > 00_scripts/datarmor_jobs/GATK_$base.sh
done


#Submit jobs
for i in $(ls 00_scripts/datarmor_jobs/GATK_*sh); do qsub $i; done


