#!/bin/bash

KAIJU_nodes=/mnt/shared/scratch/ae42909/201609_BBSRC_Diagnostics/kaiju/kaijudb/nodes.dmp
KAIJU_fmi=/mnt/shared/scratch/ae42909/201609_BBSRC_Diagnostics/kaiju/kaijudb/kaiju_db.fmi
OUTPUTNAME=kaiju_ouput

#cd /mnt/shared/scratch/ae42909/201609_BBSRC_Diagnostics/kaiju

/home/ae42909/Programs/Kaiju/kaiju-v1.5.0-linux-x86_64-static/bin/kaiju -z $THREADS -t $KAIJU_nodes -f $KAIJU_fmi -i PE1_subset and -j PE2_subset -o $OUTPUTNAME -x -v -a greedy -e 5
