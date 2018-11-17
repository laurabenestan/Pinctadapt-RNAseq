#!/bin/bash


# Clean past jobs

rm 00_scripts/datarmor_jobs/GENO_*sh


# launch scripts for Colosse
for file in $(ls 08_snps/*.bam|sed -e 's/.bam//'|sort -u)
do

base=$(basename "$file")

	toEval="cat 00_scripts/09_create.vcf | sed 's/__BASE__/$base/g'"; eval $toEval > 00_scripts/datarmor_jobs/GENO_$base.sh
done


#Submit jobs
for i in $(ls 00_scripts/datarmor_jobs/GENO_*sh); do qsub $i; done


