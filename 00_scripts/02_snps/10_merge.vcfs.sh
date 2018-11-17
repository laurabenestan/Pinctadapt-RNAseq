#!/bin/bash
#PBS -N vcf
#PBS -o 98_log_files/vcf.mergefull.err
#PBS -l walltime=23:00:00
#PBS -l mem=100g
#PBS -q omp
#PBS -l ncpus=16
#####PBS -m ea
#PBS -r n

cd $PBS_O_WORKDIR

Transcriptome="00_ressources/transcriptomes/P_margaritifera/Trinity.100aaorf.minexpr0.5.fa"
index="p_marg.dict"
FolderSNP="08_snps"
vcflist="vcf.test.list"

gatk CombineGVCFs -O "$FolderSNP"/raw.merged.full.vcf -V 08_snps/snps.noduplicates.HI.4287.003.Index_10.M5.vcf -V 08_snps/snps.noduplicates.HI.4287.003.Index_1.M1.vcf -V 08_snps/snps.noduplicates.HI.4287.003.Index_20.B2.vcf -V 08_snps/snps.noduplicates.HI.4287.003.Index_21.B5.vcf -V 08_snps/snps.noduplicates.HI.4287.003.Index_22.B3.vcf -V 08_snps/snps.noduplicates.HI.4287.003.Index_25.B4.vcf -V 08_snps/snps.noduplicates.HI.4287.003.Index_3.M2.vcf -V 08_snps/snps.noduplicates.HI.4287.003.Index_8.M4.vcf -V 08_snps/snps.noduplicates.HI.4287.003.Index_9.M3.vcf -R 00_ressources/transcriptomes/P_margaritifera/Trinity.100aaorf.minexpr0.5.fa
