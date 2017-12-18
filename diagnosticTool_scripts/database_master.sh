#!/bin/bash

#$ -cwd
#$ -j yes
#$ -pe smp 10



# Download test viruses and add them to mini db

# Pre-create output folder with symlink to already downloaded taxonomy
# (to avoid kraken downloading it all again)
rm -rf /home/ae42909/Scratch/kodoja_db/
mkdir -p /home/ae42909/Scratch/kodoja_db/krakenDB_test
ln -s /home/ae42909/viral_diagnostics/test/taxonomy /home/ae42909/Scratch/kodoja_db/krakenDB_test/taxonomy
/home/ae42909/viral_diagnostics/diagnosticTool_scripts/database_master.py --output_dir /home/ae42909/Scratch/kodoja_db/ --threads 10 --db_tag 'test' --test

# Download all viruses  and add to db plant viruses
# /home/ae42909/viral_diagnostics/diagnosticTool_scripts/database_master.py --output_dir /home/ae42909/Scratch/kodoja_db/ --threads 10 --db_tag 'plnVir'

# All vir genomes downloaded already in out_dir, add only plant viruses
# /home/ae42909/viral_diagnostics/diagnosticTool_scripts/database_master.py --output_dir /home/ae42909/Scratch/kodoja_db/ --threads 10 --db_tag 'allVir' --no_download --all_viruses
