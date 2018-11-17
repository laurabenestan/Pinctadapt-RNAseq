#!/bin/bash
#PBS -N dict.__BASE__
#PBS -o 98_log_files/dict.__BASE__.err
#PBS -l walltime=23:00:00
#PBS -l mem=30g
#####PBS -m ea
#PBS -r n


Transcriptome="00_ressources/transcriptomes/P_margaritifera/Trinity.100aaorf.minexpr0.5.fa"
index="p_marg.dict"


#launch
gatk CreateSequenceDictionary -R $Transcriptome -O $index

