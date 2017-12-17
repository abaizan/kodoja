#!/bin/bash

#$ -cwd
#$ -j yes
#$ -pe smp 10


/home/ae42909/viral_diagnostics/diagnosticTool_scripts/diagnostic_master.py -r1 /home/ae42909/data_forTesting/Apple_data/AppleVirus_SRR1089477.fastq -o /home/ae42909/Scratch/100Seq_SE/ -d1 /home/ae42909/Scratch/diagnostic_databases/krakenDB_k31_plnVir/ -d2 /home/ae42909/Scratch/diagnostic_databases/kaijuDB_plnVir/ -m 30 -t 10
