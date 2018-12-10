#!/usr/bin/env bash
#PBS -q
#PBS -l walltime=24:00:00
#PBS -l mem=115g
#PBS -N snpeff_database_building

SNPEFF=
WORKDIR=gamma
GFF=01_info_files/***.gff3
FASTA=01_info_files/***.fasta
DBNAME=***

cd $SNPEFF

#create necessary directories for analysis: a directory data in wich two other directories are created
mkdir -p data
mkdir -p data/$DBNAME
mkdir -p data/genomes


# Copy files where they need to be. GFF in teh $DBNAME directory, fasta in the genome directory
cp $WORKDIR/$GFF ./data/$DBNAME/genes.gff
cp $WORKDIR/$FASTA ./data/genomes/***.fa


#Modify the config file in order to add the new genome
echo "# Genome of Pinctada margaritifera, ***.fasta" >> snpEff.config
echo "$DBNAME.genome : $DBNAME" >> snpEff.config


java -Xmx115G -jar $SNPEFF/snpEff.jar build -gff3 -v $DBNAME


