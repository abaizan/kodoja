#!/bin/bash

#$ -cwd
#$ -j yes
#$ -pe smp 10
THREADS=10
DBNAME=/mnt/shared/scratch/ae42909/201609_BBSRC_Diagnostics/kraken/customdataset/ProteinCD_build
export PATH=/home/ae42909/Programs/Kraken/kraken-0.10.5-beta/builddir:/mnt/apps/jellyfish/1.1.11/bin:$PATH

# Download taxonomy from NCBI using kraken
#kraken-build --download-taxonomy --threads $THREADS --db $DBNAME

#kraken-build --add-to-library /home/ae42909/Scratch/emboss_transeq/concatinatedRefseqViral_transeq.fna --db $DBNAME

# Build the actual database. Supply a smaller hash size to Jellyfish using kraken-build's --jellyfish-hash-size switch. Each space in the hash table uses approximately 6.9 bytes, so using "--jellyfish-hash-size 6400M" will use a hash table size of 6.4 billion spaces and require 44.3 GB of RAM. 11400 = 80 GB. After jellyfish finished there wasn't enough memory to continue so used --max-db-size so that it reduced the database to max 70 GB.
kraken-build --build --threads $THREADS --jellyfish-hash-size 11400M --max-db-size 70 --db $DBNAME
