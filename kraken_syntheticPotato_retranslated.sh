#!/bin/bash

#$ -cwd
#$ -j yes
#$ -pe smp 4
THREADS=4
DBNAME=/home/ae42909/Scratch/kraken/kraken_analysis/retranslateDatabase
FILE1=/home/ae42909/Scratch/kraken/kraken_analysis/retranslateDatabase/retranslated_Potato_withViruses_1
FILE2=/home/ae42909/Scratch/kraken/kraken_analysis/retranslateDatabase/retranslated_Potato_withViruses_2
export PATH=/home/ae42909/Programs/Kraken/kraken-0.10.5-beta/builddir:/mnt/apps/jellyfish/1.1.11/bin:$PATH

#To classify a set of sequences (reads), use the kraken command. To obtain optimum speeds, Kraken's database should be loaded into RAM first - use the --preload switch.

kraken --preload --threads $THREADS --db $DBNAME --fasta-input --paired $FILE1 $FILE2 > PotatoViruses_CD_results

 kraken-translate --mpa-format --db $DBNAME PotatoViruses_RD_results > Potato_withViruses__RD_classified.labels
