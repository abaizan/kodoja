#!/bin/bash
#$ -cwd
#$ -j yes
export PERL5LIB=/mnt/apps/BioPerl-1.6.901
perl /home/ae42909/finished_scripts/kraken_ftp_viral.pl
