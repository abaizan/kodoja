#!/bin/bash

#$ -cwd
#$ -j yes
#$ -pe smp 1


/home/ae42909/viral_diagnostics/diagnosticTool_scripts/diagnostic_master.py -r1 /home/ae42909/Scratch/smallTest_data/100_Potato_withViruses_1.fastq -o /home/ae42909/Scratch/100Seq_PE/ -d1 /home/ae42909/Scratch/kodoja_db/krakenDB_test/ -d2 /home/ae42909/Scratch/kodoja_db/kaijuDB_test/ -r2 /home/ae42909/Scratch/smallTest_data/100_Potato_withViruses_2.fastq
