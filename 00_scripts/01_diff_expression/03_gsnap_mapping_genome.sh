#!/bin/bash
#PBS -N gsnap.__BASE__
#PBS -o gsnap.__BASE__.err
#PBS -l walltime=23:00:00
#PBS -l mem=15g
#####PBS -m ea
#PBS -l ncpus=12
#PBS -q omp
#PBS -r n

#Module load
#Samtools1.4.1

# Global variables
DATAOUTPUT="04_mapped"
DATAINPUT="03_trimmed"

# For genome
GENOMEFOLDER="path_to_genome_folder"
GENOME="gmap_genome_index"
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
	samtools sort "$DATAOUTPUT"/"$base".concordant_uniq.bam -o "$DATAOUTPUT"/"$base".sorted.bam
 	samtools index "$DATAOUTPUT"/"$base".sorted.bam  
# Clean up
    echo "Removing "$TMP"/"$base".sam"
    echo "Removing "$TMP"/"$base".bam"

