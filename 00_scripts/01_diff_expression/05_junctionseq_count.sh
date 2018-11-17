#!/usr/bin/env bash
#PBS -N jctseq.__BASE__
#PBS -o 98_log_files/log-jctseq.__BASE__.err
#PBS -l walltime=02:00:00
#PBS -m ea
#PBS -l ncpus=1
#PBS -l mem=30g
#PBS -r n

. samtools/latest/env.sh

# Move to present working dir
cd $PBS_O_WORKDIR


GTF="00_ressources/genomes/P_margaritifera/indexed_genome.gtf"
GENOME="00_ressources/genomes/P_margaritifera/Pmarg_trimmed.fasta"
base=__BASE__


# launch 
java -jar -Xmx25g softwares/QoRTs-STABLE.jar QC \
		--stranded \
		--nameSorted \
		--minMAPQ 30 \
		--maxReadLength 100 \
		--genomeFA "$GENOME" --keepMultiMapped \
		--runFunctions writeKnownSplices,writeNovelSplices,writeSpliceExon \
		04_mapped/"$base".concordant_mult.sorted.bam \
		$GTF \
		05_count_jctseq/concordant_mult/"$base"_jctseq/

