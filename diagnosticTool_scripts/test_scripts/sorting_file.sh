#!/bin/bash

#$ -cwd
#$ -j yes



cat /mnt/shared/projects/virology/201702_Stuart_local_raspvirus/static/Ample-M3_S3_R1.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > /home/ae42909/Scratch/berry/sorted1.fastq

cat /mnt/shared/projects/virology/201702_Stuart_local_raspvirus/static/Ample-M3_S3_R2.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > /home/ae42909/Scratch/berry/sorted2.fastq
