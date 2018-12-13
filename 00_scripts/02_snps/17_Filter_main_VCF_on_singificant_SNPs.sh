#!/bin/bash
#PBS -q 
#PBS -l walltime=20:00:00
#PBS -l mem=60g
#PBS -l ncpus=12
#PBS -N Sig_filtering_VCF

WORKDIR="gamma"
VCF=named_DP10_SNP_MAF0.1_miss0.1_GAMMA_genome_biallelic_IndFiltered_noComplex_imputed.vcf
SIG="07_snps/04_bayescan/locus_id.txt"
INDIR="07_snps/01_freebayes"

#---------------------------------------------
#Step 1: creation of all the necessary inputs

# We need to separate the header from the rest of the vcf file in order for the code to work:

cd $INDIR

grep -v "^#" $VCF > no_header_vcf.vcf
grep "^#" $VCF > header.txt

#---------------------------------------------------------------------------
#Step 2: creating a "single_ID" filed in the 1st column of no_header_vcf.vcf

awk '{OFS = "_" ; print $1,$2}' no_header_vcf.vcf | paste - no_header_vcf.vcf > single_ID_no_header_vcf.vcf

#Do the same for the locus_id file (columns 2 and 3) and place tham in a new file:
awk '{OFS = "_" ; print $2,$3}' $WORKDIR/$SIG > locus_id.txt


#-----------------------------------------------
#Step 3: keep only the SNPs that are significant:

#grep -f $SIG single_ID_no_header_vcf.vcf > single_ID_no_header_vcf_significant.vcf
awk 'FNR==NR {a[$1]=$0; next}; $1 in a {print a[$1]}' single_ID_no_header_vcf.vcf locus_id.txt > single_ID_no_header_vcf_significant.vcf

#----------------------------------------
#Step 4: reformat it so that it is a VCF:

#Remove first column (and the tabulation separating $1 and $2) and reattach the header
awk '{ $1=""; print $0 }' single_ID_no_header_vcf_significant.vcf | sed 's/^[ \t]+//g' - | sed 's#^[ ]##g' | sed 's# #\t#g' - | cat header.txt - > Bayescan_significant_vcf.vcf

#-----------------------------
#Step 5: clean your workspace:

rm single_ID_no_header_vcf_significant.vcf
rm single_ID_no_header_vcf_significant.vcf
rm single_ID_no_header_vcf.vcf
rm no_header_vcf.vcf
rm locus_id.txt
rm header.txt
