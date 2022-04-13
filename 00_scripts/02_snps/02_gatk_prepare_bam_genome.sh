#!/bin/bash
#PBS -N gatk.__BASE__
#PBS -o 98_log_files/gatk.__BASE__.err
#PBS -l walltime=23:00:00
#PBS -l mem=40g
#####PBS -m ea
#PBS -r n


cd $PBS_O_WORKDIR

REFERENCE="path_to_genome_folder/Pmarg.genome.masked.fa"
index="path_to_genome_folder/Pmarg.genome.masked.dict"
InDir="/home1/scratch/jleluyer/gamma/04_mapped"
OutDir="/home1/scratch/jleluyer/gamma/04_mapped/gatk"

GATK_jar="gatk-4.0.3.0/gatk-package-4.0.3.0-local.jar"
TmpDir="gamma/gatk"

base=__BASE__

# clean sam
java -Xmx35G -Djava.io.tmpdir="$TmpDir" -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar "$GATK_jar" CleanSam \
	--INPUT "$InDir"/"$base".concordant_uniq \
	--OUTPUT "$OutDir"/"$base".clean.sam 2>&1 | tee 98_log_files/error."$base".log

# remove reads problematic
grep "error" 98_log_files/gatk."$base".err |awk '{print $8}'|sed 's/,//g' >REMOVE_"$base".log

grep -v -f REMOVE_"$base".log "$OutDir"/"$base".clean.sam >"$OutDir"/tmp."$base"

mv "$OutDir"/tmp."$base" "$OutDir"/"$base".clean.sam

# Sort clean Sam file
java -Xmx35G -Djava.io.tmpdir="$TmpDir" -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar "$GATK_jar" SortSam \
	--INPUT "$OutDir"/"$base".clean.sam \
	--OUTPUT "$OutDir"/"$base".sorted.bam \
	--SORT_ORDER coordinate  

# dedup
java -Xmx35G -Djava.io.tmpdir="$TmpDir" -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar "$GATK_jar" MarkDuplicates \
        --INPUT "$OutDir"/"$base".sorted.bam \
        --OUTPUT "$OutDir"/"$base".dedup.bam \
	--REMOVE_DUPLICATES true \
        -M "$OutDir"/"$base".metrics.txt

#bam index
java -Xmx35G -Djava.io.tmpdir="$TmpDir" -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar "$GATK_jar" BuildBamIndex \
        --INPUT "$OutDir"/"$base".dedup.bam

# Split Cigar strings reads
java -Xmx35G -Djava.io.tmpdir="$TmpDir" -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar "$GATK_jar" SplitNCigarReads \
	-I "$OutDir"/"$base".dedup.bam \
	-O "$OutDir"/"$base".split.bam \
	--refactor-cigar-string true \
	-RF CigarContainsNoNOperator -RF GoodCigarReadFilter -R "$REFERENCE" 

## Base recalibration
### Need to have vcf file with known variation first !
#java -Xmx35G -Djava.io.tmpdir="$TmpDir" -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar "$GATK_jar" BaseRecalibrator \
 #  -I "$OutDir"/"$base".split.bam \
  # -R "$REFERENCE" \
  # -O "$OutDir"/"$base".recal_data.table
