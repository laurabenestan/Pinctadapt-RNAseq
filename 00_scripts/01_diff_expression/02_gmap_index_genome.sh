#!/bin/bash
#PBS -N gmap_index
#PBS -o 98_log_index/gmap_index.out
#PBS -l walltime=24:00:00
#PBS -m ea 
#PBS -l ncpus=1
#PBS -l mem=30g
#PBS -r n

GENOMEFOLDER="01_projects/gamma/01_info_files"
FASTA="01_projects/gamma/01_info_files/sspace.final.scaffolds.fasta"
GENOME="gmap_genome"

gmap_build --dir="$GENOMEFOLDER" "$FASTA" -d "$GENOME"
