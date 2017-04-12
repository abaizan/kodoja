#!/bin/bash
#$ -cwd
#$ -j yes

export PATH=/home/ae42909/Programs/Taxonomer/taxonomer_0.5/taxonomer:$PATH

### Building nucleotide database from our data for taxonomer

## 1. Create k-mer count files using kanalze of fasta files of reference nucleotide sequences:

 java -jar /home/ae42909/Programs/Kanalyze/kanalyze-2.0.0/kanalyze.jar count -k21 -d 3 -l 3 -rcanonical -o syntheticPotatoInfected -m hex /home/ae42909/synthetic_RNAseq/mappingRNAseq/concatenated_fastaFiles/Potato_withViruses_1.fastq

 
