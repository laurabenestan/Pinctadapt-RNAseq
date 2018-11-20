#!/usr/bin/env bash
#PBS -N htseq__BASE__
#PBS -o htseq__BASE__.err
#PBS -l walltime=02:00:00
#PBS -m ea
#PBS -l ncpus=1
#PBS -l mem=50g
#PBS -r n


# Move to present working dir
cd $PBS_O_WORKDIR

# install htseq
. /appli/bioinfo/htseq/0.6.1/env.sh

#Global variables
DATAINPUT="/home1/scratch/jleluyer/gamma/04_mapped/genome"
DATAOUTPUT="05_count/genome"
DATAOUTPUT_SPLICE=""

GFF_FOLDER="01_info_files"
GFF_FILE="genome_Pmarg_v2.gff3"
#launch script
base=__BASE__

# for gene expression
htseq-count -f="bam" -s="no" -r="pos" -t="gene" -i="ID" --mode="union" "$DATAINPUT"/"$base".sorted.bam "$GFF_FOLDER"/"$GFF_FILE" >>"$DATAOUTPUT"/htseq-count_"$base".txt

