#!/bin/bash

#$ -cwd
#$ -j yes
#$ -pe smp 4
THREADS=4
FASTQ1=R1.fastq
FASTQ2=R2.fastq
OUTPUTNAME=kaiju_ouput

# cd /mnt/shared/scratch/ae42909/201609_BBSRC_Diagnostics/kaiju

# If you chose options -n or -e in makeDB.sh, then use -f kaiju_db_nr.fmi or -f kaiju_db_nr_euk.fmi

/home/ae42909/Programs/Kaiju/kaiju-v1.5.0-linux-x86_64-static/bin/kaiju -z $THREADS -t kaijudb/nodes.dmp -f kaijudb/kaiju_db.fmi -i $FASTQ1 and -j $FASTQ2 -o $OUTPUTNAME -x -v -a greedy -e 5

# The default run mode is MEM, which only considers exact matches. For using the Greedy mode, which allows mismatches, set the mode via the option -a and the number of allowed substitutions using option -e e.g. -a greedy -e 5
# The cutoffs for minimum required match length and match score can be changed using the options -m (default: 11) and -s (default: 65)
# Option -x can be used to enable filtering of query sequences containing low-complexity regions by using the SEG algorithm from the blast+ package. Enabling this option is always recommended in order to avoid false positive matches caused by spurious matches due to simple repeat patterns or other sequencing noise. -v option is to enable verbose results.

# The default run mode is MEM, which only considers exact matches. For using the Greedy mode, which allows mismatches, set the mode via the option -a and the number of allowed substitutions using option -e
