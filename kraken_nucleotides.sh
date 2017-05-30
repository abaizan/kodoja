
#!/bin/bash

DBNAME=/mnt/shared/scratch/ae42909/201609_BBSRC_Diagnostics/kraken/customdataset/CD_build
export PATH=/home/ae42909/Programs/Kraken/kraken-0.10.5-beta/builddir:/mnt/apps/jellyfish/1.1.11/bin:$PATH

#To classify a set of sequences (reads), use the kraken command. To obtain optimum speeds, Kraken's database should be loaded into RAM first - use the --preload switch.

kraken --preload --threads $THREADS --db $DBNAME --fastq-input --paired PE1_renamed PE2_renamed > kraken_nt_results

kraken-translate --mpa-format --db $DBNAME kraken_nt_results > kraken_nt_labels
