#!/usr/bin/env bash
#PBS -N htseq__BASE__
#PBS -o 98_log_files/htseq__BASE__.err
#PBS -l walltime=02:00:00
#PBS -m ea
#PBS -l ncpus=1
#PBS -l mem=20g
#PBS -r n


# Move to present working dir
cd $PBS_O_WORKDIR

# Module load htseq
#htseq0.6.1

#Global variables
DATAINPUT="04_mapped"
DATAOUTPUT="05_count"

GFF_FOLDER="path_to_genome_folder"
GFF_FILE="Pmarg.genome.gff3"
#launch script
base=__BASE__

# for gene expression
htseq-count -f "bam" -s "no" -r "pos" -t "gene" -i "ID" --mode "union" "$DATAINPUT"/"$base".sorted.bam "$GFF_FOLDER"/"$GFF_FILE" >>"$DATAOUTPUT"/htseq-count_"$base".txt
