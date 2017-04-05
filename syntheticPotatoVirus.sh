#!/bin/bash

#$ -cwd
#$ -j yes
#$ -pe smp 4
THREADS=4
DBNAME=/mnt/shared/scratch/ae42909/201609_BBSRC_Diagnostics/kraken/customdataset/CD_build
export PATH=/home/ae42909/Programs/Kraken/kraken-0.10.5-beta/builddir:/mnt/apps/jellyfish/1.1.11/bin:$PATH

#To classify a set of sequences (reads), use the kraken command. To obtain optimum speeds, Kraken's database should be loaded into RAM first - use the --preload switch.

kraken --preload --threads $THREADS --db $DBNAME --fastq-input --paired /home/ae42909/synthetic_RNAseq/mappingRNAseq/concatenated_fastaFiles/Potato_withViruses_1.fastq /home/ae42909/synthetic_RNAseq/mappingRNAseq/concatenated_fastaFiles/Potato_withViruses_2.fastq > PotatoViruses_CD_results

 kraken-translate --mpa-format --db $DBNAME PotatoViruses_CD_results > Potato_withViruses_classified.labels
