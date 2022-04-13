#!/bin/bash
#PBS -N gmap_index
#PBS -o 98_log_index/gmap_index.out
#PBS -l walltime=24:00:00
#PBS -m ea 
#PBS -l ncpus=1
#PBS -l mem=30g
#PBS -r n

GENOMEFOLDER="path_to_genome_folder"
FASTA="path_to_genome_folder/genome.masked.fa"
GENOME="gmap_genome_index"

gmap_build --dir="$GENOMEFOLDER" "$FASTA" -d "$GENOME"
