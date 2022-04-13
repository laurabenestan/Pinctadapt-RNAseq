#!/usr/bin/env bash
#PBS -N beagle
#PBS -o 98_log_files/beagle.err
#PBS -l mem=110G
#PBS -q omp
#PBS -l ncpus=14
#PBS -l walltime=20:00:00
#PBS -j oe


#Module load
#beagle4.0_06Jun17

# Look at the manual for help : http://faculty.washington.edu/browning/beagle/beagle_4.0_06Jun17.pdf

cd $PBS_0_WORKDIR

INPUT="06_outlier/combined.minDP10.MAF0.1.MISS0.9.Q30.biallelic.subset.recode.vcf"
OUTPUT="06_outlier/beagle.imputed"
NCPUS=14
#BEAGLE="beagle.29May21.d6d.jar"
BEAGLE="beagle.27Jan18.7e1.jar"
java -Xmx100G -Djava.io.tmpdir=$HOME -jar $BEAGLE nthreads=$NCPUS gt="$INPUT" out="$OUTPUT"
 
