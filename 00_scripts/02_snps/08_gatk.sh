#!/bin/bash
#PBS -N vcf.__BASE__
#PBS -o 98_log_files/vcf.__BASE__.err
#PBS -l walltime=23:00:00
#PBS -l mem=60g
#PBS -q omp
#PBS -l ncpus=16
#####PBS -m ea
#PBS -r n

cd $PBS_O_WORKDIR

Transcriptome="00_ressources/transcriptomes/P_margaritifera/Trinity.100aaorf.minexpr0.5.fa"
index="p_marg.dict"
FolderSNP="08_snps"
base=__BASE__
FolderIn="04_mapped/transcriptome/subset"



#clean sam

#samtools view -H "$FolderIn"/"$base".concordant_uniq.sorted.bam -o "$FolderIn"/"$base".concordant_uniq.sorted.sam

#gatk CleanSam --INPUT "$FolderIn"/"$base".concordant_uniq.sorted.sam --OUTPUT "$FolderIn"/"$base".clean.sam
#rm "$FolderIn"/"$base".concordant_uniq.sorted.sam

#samtools view -Sb "$FolderIn"/"$base".clean.sam > "$FolderIn"/"$base".bam

#rm "$FolderIn"/"$base".clean.sam

#samtools sort "$FolderIn"/"$base".bam -o "$FolderIn"/"$base".concordant_uniq.sorted2.bam

#samtools index "$FolderIn"/"$base".concordant_uniq.sorted2.bam

#gatk FilterSamReads --INPUT "$FolderIn"/"$base".concordant_uniq.sorted2.bam \
#		--OUTPUT "$FolderIn"/"$base".concordant_uniq.sorted3.bam \
#		--FILTER excludeReadList \
#		--READ_LIST_FILE list.exclude."$base"

#samtools view -h "$FolderIn"/"$base".concordant_uniq.sorted2.bam | grep -vf list.exclude."$base" | samtools view -bS -o "$FolderIn"/"$base".concordant_uniq.sorted3.bam -

# launck markdupli
gatk MarkDuplicates --INPUT "$FolderIn"/"$base".concordant_uniq.sorted.bam \
			--METRICS_FILE output.metrics."$base" \
			--OUTPUT "$FolderSNP"/noduplicates."$base".bam \
			--REMOVE_DUPLICATES true


