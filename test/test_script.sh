#!/bin/bash

#$ -cwd
#$ -j yes
#$ -pe smp 2

# Enable strict bash mode - we want to abort on the first error, see also
# http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail

# Create test database - ought work but we run out of disk space on TravisCI
# while downloading and unzipping the taxonomy data.
#
diagnosticTool_scripts/database_master.py -o test/building_db/ -t 1 --db_tag 'test' --test
ls test/building_db
diff test/example_db/krakenDB_test/database.idx test/building_db/krakenDB_test/database.idx
# TODO: Compare more of the output files?

# Run diagnostic tool
diagnosticTool_scripts/diagnostic_master.py -r1 test/data/testData_1.fastq -o test/PE_test/ -d1 test/example_db/krakenDB_test/ -d2 test/example_db/kaijuDB_test/ -r2 ./test/data/testData_1.fastq -t 2
diff test/PE_test/virus_table.txt test/data/example_virus_table.txt
rm -r test/PE_test/

diagnosticTool_scripts/diagnostic_master.py -r1 test/data/testData_1.fastq -o test/SE_test/ -d1 test/example_db/krakenDB_test/ -d2 test/example_db/kaijuDB_test/ -t 2
diff test/SE_test/virus_table.txt test/data/example_virus_table.txt
rm -r test/SE_test/

diagnosticTool_scripts/diagnostic_master.py -r1 test/data/testData_1.fasta -o test/fasta_test/ -d1 test/example_db/krakenDB_test/ -d2 test/example_db/kaijuDB_test/ -r2 ./test/data/testData_1.fasta -f fasta
diff test/fasta_test/virus_table.txt test/data/fastaEx_virus_table.txt
rm -r test/fasta_test/


