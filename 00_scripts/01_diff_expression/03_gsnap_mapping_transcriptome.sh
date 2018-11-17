#!/bin/bash
#PBS -N gsnap.trans.__BASE__
#PBS -o gsnap.trans.__BASE__.err
#PBS -l walltime=24:00:00
#PBS -l mem=30g
#####PBS -m ea
#PBS -l ncpus=12
#PBS -q omp
#PBS -r n


. samtools/1.4.1/env.sh

# Global variables
DATAOUTPUT="gamma/04_mapped"
DATAINPUT="gamma/03_trimmed"

# For transcriptome
GENOMEFOLDER="00_ressources/transcriptomes/P_margaritifera"
GENOME="gmap_0.5_pmargaritifera"

# For genome
#GENOMEFOLDER="00_ressources/genomes/P_margaritifera"
#GENOME="indexed_genome"
platform="Illumina"

#move to present working dir
cd $PBS_O_WORKDIR

base=__BASE__

    # Align reads
    echo "Aligning $base"

gsnap --gunzip -t 12 -A sam --min-coverage=0.5 \
	--dir="$GENOMEFOLDER" -d "$GENOME" \
	--split-output="$DATAOUTPUT"/"$base" \
	--max-mismatches=5 --novelsplicing=1 \
	--read-group-id="$base" \
	 --read-group-platform="$platform" \
	"$DATAINPUT"/"$base"_R1.paired.fastq.gz "$DATAINPUT"/"$base"_R2.paired.fastq.gz
    
# Create bam file
    echo "Creating bam for $base"
samtools view -b "$DATAOUTPUT"/"$base".concordant_uniq >"$DATAOUTPUT"/"$base".concordant_uniq.bam
samtools sort "$DATAOUTPUT"/"$base".concordant_uniq.bam -o "$DATAOUTPUT"/"$base".concordant_uniq.sorted.bam
samtools index "$DATAOUTPUT"/"$base".concordant_uniq.sorted.bam
    
# Clean up
    echo "Removing "$TMP"/"$base".sam"
    echo "Removing "$TMP"/"$base".bam"

   	rm "$DATAOUTPUT"/"$base".sam
    	rm "$DATAOUTPUT"/"$base".bam
