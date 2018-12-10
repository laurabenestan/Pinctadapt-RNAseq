#!/bin/bash
#PBS -q
#PBS -l walltime=20:00:00
#PBS -l mem=60g
#PBS -l ncpus=12
#PBS -N make_bayescan_input

WORKDIR=gamma
INDIR=07_snps/04_bayescan
RCODE=00_scripts/02_snps/make_bayescan_input.R


source activate vcfR

Rscript --vanilla $WORKDIR/$RCODE &> $WORKDIR/00_scripts/make_bayescan_input.R.out ;

source deactivate vcfR

cd $WORKDIR/$INDIR
echo '[populations]=2' >> head.txt ;
sed -i '1 i\[loci]=27394' head.txt ;

sed -i '1 i\[pop]=1' Gambier_bayescan_input.txt ;
sed -i '1 i\[pop]=2' Marquesas_bayescan_input.txt ;

cat head.txt Gambier_bayescan_input.txt Marquesas_bayescan_input.txt > Bayescan_input.txt ;

rm head.txt
rm Gambier_bayescan_input.txt
rm Marquesas_bayescan_input.txt



