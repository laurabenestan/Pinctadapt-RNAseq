#!/bin/bash
#PBS -N dict.__BASE__
#PBS -o 98_log_files/gatk.prepare.ref.err
#PBS -l walltime=23:00:00
#PBS -l mem=70g
#PBS -r n

cd $PBS_O_WORKDIR

REFERENCE="work/01_projects/gamma/01_info_files/sspace.final.scaffolds.fasta"
index="work/01_projects/gamma/01_info_files/sspace.final.scaffolds.dict"

TmpDir="scracth/gamma/gatk/"
GATK_jar="home/softwares/gatk-4.0.3.0/gatk-package-4.0.3.0-local.jar"


#launch
java -Xmx40G -Djava.io.tmpdir="$TmpDir" -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar "$GATK_jar" CreateSequenceDictionary \
	-R $REFERENCE -O $index

# Samtools index fasta
samtools faidx $REFERENCE
