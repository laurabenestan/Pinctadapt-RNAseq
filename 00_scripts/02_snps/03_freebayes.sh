#!/usr/bin/env bash
#PBS -q omp
#PBS -l walltime=36:00:00
#PBS -l ncpus=2
#PBS -l mem=20g


cd $PBS_O_WORKDIR

#Working directories, input and output files
ls -1 scracth/gamma/04_mapped/genome/gatk/*split.bam|grep -v "A9.17" >01_info_files/list.individuals.freebayes.txt
data="01_info_files/list.individuals.freebayes.txt"
outdir=07_snps
ref="01_info_files/sspace.final.scaffolds.fasta"
tag="snp_freebayes_v1"
tmp=scracth


#Freebayes parameters 
nAlleles="0"
minMapQ="30"
minCOV=10
Ploidy=2


. module/freebayes/latest/env.sh 

#Calling SNP with Freebayes on trimmed bam file
echo "Running Freebayes on ${data} samples..."

time freebayes -f $ref \
                --use-best-n-alleles $nAlleles \
                --min-mapping-quality $minMapQ \
                --no-indels \
   		--no-complex \
		--min-coverage $minCOV \
                --genotype-qualities \
                --bam-list ${data} \
                --vcf ${outdir}/${tag}.vcf 
echo "Running Freebayes on ${data} samples done."

