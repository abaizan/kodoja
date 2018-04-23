#!/bin/bash

#$ -cwd
#$ -j yes
#$ -pe smp 4


/home/ae42909/kodoja/viral_diagnostics/diagnosticTool_scripts/kodoja_search.py -r1 /home/ae42909/data_forTesting/Algae_data/AlgaeVirus60_SRR924352.fastq -o /home/ae42909/Scratch/Algae_60/ -d1 /home/ae42909/Scratch/kodoja_db/krakenDB_algae_554065/ -d2 /home/ae42909/Scratch/kodoja_db/kaijuDB_algae_554065/ -m 30 -t 4
