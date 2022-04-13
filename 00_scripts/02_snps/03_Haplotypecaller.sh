#!/usr/bin/env bash
#PBS -N haplotypecaller
#PBS -o 98_log_files/haplo.__BASE__.err
#PBS -l walltime=50:00:00
#PBS -l mem=15g
#PBS -r n


cd $PBS_O_WORKDIR
Reference="path_to_genome_folder/Pmarg.genome.masked.fa"
index="path_to_genome_folder/Pmarg.genome.masked.dict"
InDir="/home1/scratch/jleluyer/gamma/04_mapped"
OutDir="/home1/scratch/jleluyer/gamma/04_mapped/gatk"

GATK_jar="gatk-4.0.3.0/gatk-package-4.0.3.0-local.jar"
TmpDir="gamma/gatk"
base=__BASE__

java -Xmx10G -Djava.io.tmpdir="$TmpDir" -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar "$GATK_jar" HaplotypeCaller \
        -I "$OutDir"/"$base".split.bam \
        -R "$Reference" \
	-O "$OutDir"/"$base".g.vcf.gz \
	-ERC GVCF
