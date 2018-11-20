#!/bin/bash
#PBS -N gsnap.__BASE__
#PBS -o gsnap.__BASE__.err
#PBS -l walltime=23:00:00
#PBS -l mem=30g
#####PBS -m ea
#PBS -l ncpus=12
#PBS -q omp
#PBS -r n


. /appli/bioinfo/samtools/1.4.1/env.sh

# Global variables
DATAOUTPUT="/home1/scratch/jleluyer/gamma/04_mapped/genome"
DATAINPUT="/home1/scratch/jleluyer/gamma/03_trimmed"

# For genome
GENOMEFOLDER="/home1/datawork/jleluyer/01_projects/gamma/01_info_files/"
GENOME="gmap_genome"
platform="Illumina"

#move to present working dir
cd $PBS_O_WORKDIR

base=__BASE__

    # Align reads
    echo "Aligning $base"
 gsnap --gunzip -t "$NCPUS" -A sam --min-coverage=0.95 \
	--dir="$GENOMEFOLDER" -d "$GENOME" \
       	--max-mismatches=2 --novelsplicing=1 \
	--split-output="$DATAOUTPUT"/"$base" \
	--read-group-id="$base" \
	--read-group-platform="$platform" \
	"$DATAINPUT"/"$base"_R1.paired.fastq.gz "$DATAINPUT"/"$base"_R2.paired.fastq.gz

# concatenate sam
	samtools view -b "$DATAOUTPUT"/"$base".concordant_uniq >"$DATAOUTPUT"/"$base".concordant_uniq.bam
# name sorting bam
	echo "Creating sorted bam for $base"
	samtools sort -n "$DATAOUTPUT"/"$base".concordant_uniq.bam -o "$DATAOUTPUT"/"$base".sorted.bam
 	samtools index "$DATAOUTPUT"/"$base".sorted.bam  
# Clean up
    echo "Removing "$TMP"/"$base".sam"
    echo "Removing "$TMP"/"$base".bam"

