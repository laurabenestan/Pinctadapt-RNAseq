#!/usr/bin/env bash
#PBS -q
#PBS -l walltime=24:00:00
##PBS -l ncpus=4
#PBS -l mem=40g
#PBS -N vcftools_stat

WORKDIR=gamma
INDIR=07_snps/01_freebayes
VCFTOOLSENV=
VCF1=GAMMA_noComplex_Gambier.recode.vcf
VCF2=GAMMA_noComplex_marquesas.recode.vcf

cd $WORKDIR/$INDIR
mkdir -p stats_popVCF


$VCFTOOLSENV/vcftools --site-pi --vcf $VCF1
mv out.sites.pi stats_popVCF/Gambier_sites_pi.txt
$VCFTOOLSENV/vcftools --site-pi --vcf $VCF2
mv out.sites.pi stats_popVCF/Marquesas_sites_pi.txt
$VCFTOOLSENV/vcftools --het --vcf $VCF1
mv out.het stats_popVCF/Gambier_het.txt
$VCFTOOLSENV/vcftools --het --vcf $VCF2
mv out.het stats_popVCF/Marquesas_het.txt

