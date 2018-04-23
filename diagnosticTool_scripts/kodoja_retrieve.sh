#!/bin/bash

#$ -cwd
#$ -j yes
#$ -pe smp 1

/home/ae42909/kodoja/viral_diagnostics/diagnosticTool_scripts/kodoja_retrieve.py -r1 '/home/ae42909/kodoja/Stuart_data/Dee-B5_S1_R1.fastq' -o '/home/ae42909/kodoja/diagnostic_results/Stuarts_data/Dee-B5_S1/'  -r2 '/home/ae42909/kodoja/Stuart_data/Dee-B5_S1_R2.fastq' -f 'fastq' -t 12275
