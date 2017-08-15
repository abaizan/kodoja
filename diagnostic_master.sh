#!/bin/bash

#$ -cwd
#$ -j yes
#$ -pe smp 4
export PATH=/home/ae42909/Programs/Kraken/kraken-0.10.5-beta/builddir:/mnt/apps/jellyfish/1.1.11/bin:$PATH
export PERL5LIB=/mnt/apps/BioPerl-1.6.901

python /home/ae42909/scripts_inProcess/diagnostic_master.py
