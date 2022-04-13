#!/usr/bin/env bash
#PBS -N combineGvcf
#PBS -o 98_log_files/combine.err
#PBS -l walltime=50:00:00
#PBS -l mem=15g
#PBS -r n



cd $PBS_O_WORKDIR
Reference="path_to_genome_folder/Pmarg.genome.masked.fa"
index="path_to_genome_folder/Pmarg.genome.masked.dict"
OutDir="04_mapped/gatk"

GATK_jar="gatk-4.0.3.0/gatk-package-4.0.3.0-local.jar"
TmpDir="gamma/gatk"
base=__BASE__

#Combine gVCF
#java -Xmx10G -Djava.io.tmpdir="$TmpDir" -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar "$GATK_jar" CombineGVCFs \
 #  -R "$Reference" \
 #  -O "$OutDir"/combined.raw.gvcf.gz \
#	--variant 04_mapped/gatk/HI.4880.006.NEBNext_dual_i7_A10---NEBNext_dual_i5_A10.1.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.006.NEBNext_dual_i7_A11---NEBNext_dual_i5_A11.21.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.006.NEBNext_dual_i7_A9---NEBNext_dual_i5_A9.17.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.006.NEBNext_dual_i7_B10---NEBNext_dual_i5_B10.50.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.006.NEBNext_dual_i7_B11---NEBNext_dual_i5_B11.51.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.006.NEBNext_dual_i7_B9---NEBNext_dual_i5_B9.49.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.006.NEBNext_dual_i7_C10---NEBNext_dual_i5_C10.56.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.006.NEBNext_dual_i7_C11---NEBNext_dual_i5_C11.37.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.006.NEBNext_dual_i7_C9---NEBNext_dual_i5_C9.55.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.006.NEBNext_dual_i7_D10---NEBNext_dual_i5_D10.42.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.006.NEBNext_dual_i7_D11---NEBNext_dual_i5_D11.43.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.006.NEBNext_dual_i7_D9---NEBNext_dual_i5_D9.41.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.006.NEBNext_dual_i7_E10---NEBNext_dual_i5_E10.48.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.006.NEBNext_dual_i7_E11---NEBNext_dual_i5_E11.113.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.006.NEBNext_dual_i7_E9---NEBNext_dual_i5_E9.47.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.006.NEBNext_dual_i7_F10---NEBNext_dual_i5_F10.118.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.006.NEBNext_dual_i7_F11---NEBNext_dual_i5_F11.119.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.006.NEBNext_dual_i7_F9---NEBNext_dual_i5_F9.117.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.006.NEBNext_dual_i7_G10---NEBNext_dual_i5_G10.121.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.006.NEBNext_dual_i7_G11---NEBNext_dual_i5_G11.106.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.006.NEBNext_dual_i7_G9---NEBNext_dual_i5_G9.124.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.006.NEBNext_dual_i7_H10---NEBNext_dual_i5_H10.111.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.006.NEBNext_dual_i7_H11---NEBNext_dual_i5_H11.112.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.006.NEBNext_dual_i7_H9---NEBNext_dual_i5_H9.110.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.007.NEBNext_dual_i7_A12---NEBNext_dual_i5_A12.5.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.007.NEBNext_dual_i7_A1---NEBNext_dual_i5_A1.33.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.007.NEBNext_dual_i7_A2---NEBNext_dual_i5_A2.40.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.007.NEBNext_dual_i7_B12---NEBNext_dual_i5_B12.52.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.007.NEBNext_dual_i7_B1---NEBNext_dual_i5_B1.39.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.007.NEBNext_dual_i7_B2---NEBNext_dual_i5_B2.62.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.007.NEBNext_dual_i7_C12---NEBNext_dual_i5_C12.38.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.007.NEBNext_dual_i7_C1---NEBNext_dual_i5_C1.61.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.007.NEBNext_dual_i7_C2---NEBNext_dual_i5_C2.116.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.007.NEBNext_dual_i7_D12---NEBNext_dual_i5_D12.44.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.007.NEBNext_dual_i7_D2---NEBNext_dual_i5_D2.102.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.007.NEBNext_dual_i7_E12---NEBNext_dual_i5_E12.114.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.007.NEBNext_dual_i7_E1---NEBNext_dual_i5_E1.101.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.007.NEBNext_dual_i7_E2---NEBNext_dual_i5_E2.105.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.007.NEBNext_dual_i7_F12---NEBNext_dual_i5_F12.120.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.007.NEBNext_dual_i7_F1---NEBNext_dual_i5_F1.108.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.007.NEBNext_dual_i7_F2---NEBNext_dual_i5_F2.29.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.007.NEBNext_dual_i7_G12---NEBNext_dual_i5_G12.107.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.007.NEBNext_dual_i7_G1---NEBNext_dual_i5_G1.9.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.007.NEBNext_dual_i7_G2---NEBNext_dual_i5_G2.35.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.007.NEBNext_dual_i7_H12---NEBNext_dual_i5_H12.25.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.007.NEBNext_dual_i7_H1---NEBNext_dual_i5_H1.34.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.007.NEBNext_dual_i7_H2---NEBNext_dual_i5_H2.57.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.008.NEBNext_dual_i7_A3---NEBNext_dual_i5_A3.63.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.008.NEBNext_dual_i7_A4---NEBNext_dual_i5_A4.98.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.008.NEBNext_dual_i7_A5---NEBNext_dual_i5_A5.122.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.008.NEBNext_dual_i7_B3---NEBNext_dual_i5_B3.97.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.008.NEBNext_dual_i7_B4---NEBNext_dual_i5_B4.104.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.008.NEBNext_dual_i7_B5---NEBNext_dual_i5_B5.128.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.008.NEBNext_dual_i7_C3---NEBNext_dual_i5_C3.103.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.008.NEBNext_dual_i7_C4---NEBNext_dual_i5_C4.126.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.008.NEBNext_dual_i7_C5---NEBNext_dual_i5_C5.16.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.008.NEBNext_dual_i7_D12---NEBNext_dual_i5_D12.109.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.008.NEBNext_dual_i7_D3---NEBNext_dual_i5_D3.125.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.008.NEBNext_dual_i7_D4---NEBNext_dual_i5_D4.13.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.008.NEBNext_dual_i7_D5---NEBNext_dual_i5_D5.54.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.008.NEBNext_dual_i7_E3---NEBNext_dual_i5_E3.32.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.008.NEBNext_dual_i7_E4---NEBNext_dual_i5_E4.53.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.008.NEBNext_dual_i7_E5---NEBNext_dual_i5_E5.60.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.008.NEBNext_dual_i7_F3---NEBNext_dual_i5_F3.36.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.008.NEBNext_dual_i7_F4---NEBNext_dual_i5_F4.59.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.008.NEBNext_dual_i7_F5---NEBNext_dual_i5_F5.46.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.008.NEBNext_dual_i7_G4---NEBNext_dual_i5_G4.45.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.008.NEBNext_dual_i7_G5---NEBNext_dual_i5_G5.100.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.008.NEBNext_dual_i7_H3---NEBNext_dual_i5_H3.64.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.008.NEBNext_dual_i7_H4---NEBNext_dual_i5_H4.99.g.vcf.gz \
#	--variant 04_mapped/gatk/HI.4880.008.NEBNext_dual_i7_H5---NEBNext_dual_i5_H5.123.g.vcf.gz

## Create vcf output
java -Xmx10G -Djava.io.tmpdir="$TmpDir" -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar "$GATK_jar" GenotypeGVCFs \
   -R "$Reference" \
   -V "$OutDir"/combined.raw.gvcf.gz \
   -O "$OutDir"/combined.raw.vcf.gz
