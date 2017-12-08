#!/bin/bash

#$ -cwd
#$ -j yes
#$ -pe smp 4

export PATH=/home/ae42909/Programs/Kraken/kraken-1.0/builddir:/mnt/apps/jellyfish/1.1.11/bin:/home/ae42909/Programs/Kaiju/kaiju-v1.5.0-linux-x86_64-static/bin:/mnt/apps/trimmomatic/0.32:$PATH

python /home/ae42909/viral_diagnostics/diagnosticTool_scripts/diagnostic_master.py -r1 /home/ae42909/Scratch/smallTest_data/100_Potato_withViruses_1.fastq -r2 /home/ae42909/Scratch/smallTest_data/100_Potato_withViruses_2.fastq -o /home/ae42909/Scratch/100Seq_PE/ -d1 /home/ae42909/Scratch/diagnostic_databases/krakenDB_k31_plnVir/ -d2 /home/ae42909/Scratch/diagnostic_databases/kaijuDB_plnVir/ -t 4
