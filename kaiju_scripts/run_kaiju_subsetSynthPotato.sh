#!/bin/bash

#$ -cwd
#$ -j yes
#$ -pe smp 4
THREADS=4
FASTQ1=/home/ae42909/Scratch/kraken/kraken_analysis/customDatabase/subset_Potato_withViruses_1.fastq
FASTQ2=/home/ae42909/Scratch/kraken/kraken_analysis/customDatabase/subset_Potato_withViruses_2.fastq
OUTPUTNAME=subsetSynthPotato_ouput

#cd /mnt/shared/scratch/ae42909/201609_BBSRC_Diagnostics/kaiju

/home/ae42909/Programs/Kaiju/kaiju-v1.5.0-linux-x86_64-static/bin/kaiju -z $THREADS -t kaijudb/nodes.dmp -f kaijudb/kaiju_db.fmi -i $FASTQ1 and -j $FASTQ2 -o $OUTPUTNAME -x -v -a greedy -e 5
