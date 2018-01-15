#!/bin/bash

#$ -cwd
#$ -j yes
#$ -pe smp 10



# Download test viruses and add them to mini db
# /home/ae42909/viral_diagnostics/diagnosticTool_scripts/database_master.py --output_dir /home/ae42909/Scratch/kodoja_db/ --threads 10 --db_tag 'test' --test --kraken_tax /home/ae42909/Scratch/kraken_taxonomy/

# Download all virus genomes but only add plant viruses to databases (no host)
# /home/ae42909/viral_diagnostics/diagnosticTool_scripts/database_master.py --output_dir /home/ae42909/Scratch/kodoja_db/ --threads 10 --db_tag 'plnVir'

# All vir genomes downloaded already in out_dir, add only plant viruses (no host)
# /home/ae42909/viral_diagnostics/diagnosticTool_scripts/database_master.py --output_dir /home/ae42909/Scratch/kodoja_db/ --threads 10 --db_tag 'allVir' --no_download --all_viruses


