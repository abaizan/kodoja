#!/bin/bash

#$ -cwd
#$ -j yes
#$ -pe smp 1

python /home/ae42909/viral_diagnostics/diagnosticTool_scripts/diagnostic_seqRetrive.py -r1 '/home/ae42909/Scratch/smallTest_data/100_Potato_withViruses_1.fastq' -o '/home/ae42909/Scratch/100Seq_krakenDB_viral/' -r2 '/home/ae42909/Scratch/smallTest_data/100_Potato_withViruses_2.fastq' -f 'fastq' -t 137758
