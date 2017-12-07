#!/bin/bash

#$ -cwd
#$ -j yes
#$ -pe smp 1

python diagnosticTool_scripts/diagnostic_master.py -r1 data/testData_1.fastq -r2 data/testData_2.fastq -o 
diff data/example_virus_table.txt 
