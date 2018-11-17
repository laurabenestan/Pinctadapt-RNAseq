#!/usr/bin/env bash
#PBS -q omp
#PBS -l walltime=36:00:00
#PBS -l ncpus=2
#PBS -l mem=20g


cd $PBS_O_WORKDIR

#Working directories, input and output files
data=list.input.freebayes.txt
outdir=08_snps
ref="00_ressources/transcriptomes/P_margaritifera/Trinity.100aaorf.minexpr0.5.fa"
tag="snp_subset_freebayes_v1"
tmp=/home1/scratch/jleluyer


#Freebayes parameters 
nAlleles="0"
minMapQ="30"
minCOV=10
Ploidy=10


. freebayes/latest/env.sh 

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
                --vcf ${outdir}/diploid.${tag}.vcf 
echo "Running Freebayes on ${data} samples done."


#    --pooled-discrete --ploidy $Ploidy \

