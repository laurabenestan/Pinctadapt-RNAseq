#!/usr/bin/env bash
#PBS -q
#PBS -l walltime=24:00:00
#PBS -l mem=115g
#PBS -N beagle

WORKDIR=gamma
INDIR=07_snps/02_freebayes
BEAGLE=
TEMP_DIR="Path to your temporary output folder (large volume needed)"

cd $WORKDIR/$INDIR

cat DP10_SNP_MAF0.1_miss0.1_GAMMA_genome_biallelic_IndFiltered.recode_noComplex.vcf | sed 's#.:.:.:.:.:.:.:.:.#./.:.:.:.:.:.:.:.:.#g' > vcf_input_beagle.vcf ;


java -Xmx115G -Djava.io.tmpdir=$TEMP_DIR -jar $BEAGLE gt=vcf_input_beagle.vcf out=DP10_SNP_MAF0.1_miss0.1_GAMMA_genome_biallelic_IndFiltered_noComplex_imputed ;

rm vcf_input_beagle.vcf ;

