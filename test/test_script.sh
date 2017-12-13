#!/bin/bash

#$ -cwd
#$ -j yes
#$ -pe smp 2

# Create test database
# diagnosticTool_scripts/database_master.py -o test/example_db/ -t 1 -q -a 'test' -t 2

# Run diagnostic tool
diagnosticTool_scripts/diagnostic_master.py -r1 test/data/testData_1.fastq -o test/PE_test/ -d1 test/databases/krakenDB_test/ -d2 test/databases/kaijuDB_test/ -r2 ./test/data/testData_1.fastq -t 2
diff test/PE_test/virus_table.txt test/data/example_virus_table.txt
rm -r test/PE_test/

diagnosticTool_scripts/diagnostic_master.py -r1 test/data/testData_1.fastq -o test/SE_test/ -d1 test/databases/krakenDB_test/ -d2 test/databases/kaijuDB_test/ -t 2
diff test/SE_test/virus_table.txt test/data/example_virus_table.txt
rm -r test/SE_test/

diagnosticTool_scripts/diagnostic_master.py -r1 test/data/testData_1.fasta -o test/fasta_test/ -d1 test/databases/krakenDB_test/ -d2 test/databases/kaijuDB_test/ -r2 ./test/data/testData_1.fasta -f fasta
diff test/fasta_test/virus_table.txt test/data/fastaEx_virus_table.txt
rm -r test/fasta_test/


