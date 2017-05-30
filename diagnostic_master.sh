#!/bin/bash

#$ -cwd
#$ -j yes
#$ -pe smp 4
set -e

export THREADS=4
export FILE1=/home/ae42909/synthetic_RNAseq/mappingRNAseq/concatenated_fastaFiles/Potato_withViruses_1.fastq
export FILE2=/home/ae42909/synthetic_RNAseq/mappingRNAseq/concatenated_fastaFiles/Potato_withViruses_2.fastq
export FORMAT=fastq
export OUTDIR=/home/ae42909/Scratch/synthPotato_pipeline/
export TRIMLEN=50

# Test data format
python /home/ae42909/scripts_inProcess/testing_format.py

# QC and trim data
/home/ae42909/scripts_inProcess/fastqc.sh $TRIMLEN $FILE1 $FILE2

# Order and replace names
/home/ae42909/scripts_inProcess/rename_sequences.sh

# Kraken nucleotide analysis
/home/ae42909/scripts_inProcess/kraken_nucleotides.sh

# Pull out sequences for re-analysis
python /home/ae42909/scripts_inProcess/sequences_reanalysis.py

# Kaiju classification
/home/ae42909/scripts_inProcess/kaiju_nucleotide.sh

# Emboss translation
/home/ae42909/scripts_inProcess/translate_nucleotide.sh

# Backtranslate
python /home/ae42909/scripts_inProcess/backtranslate.py

# Kraken analysis of backtranslate
/home/ae42909/scripts_inProcess/kraken_aminoacids.sh

# Result analysis
