#!/bin/bash

#$ -cwd
#$ -j yes
#$ -pe smp 2

# Enable strict bash mode - we want to abort on the first error, see also
# http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail

if [ -f "test/test_script.sh" ]
then
    echo "Good, in the expected directory. Starting tests..."
else
    echo "ERROR. Run this from the GitHub repository root directory."
    echo "Please don't run from the test directory itself."
    exit 1
fi

# Create test database - ought work but we run out of disk space on TravisCI
# while downloading and unzipping the taxonomy data.
#
diagnosticTool_scripts/kodoja_build.py -o test/building_db/ -t 1 --db_tag 'test' --test
ls test/building_db
diff test/example_db/krakenDB_test/database.idx test/building_db/krakenDB_test/database.idx
# TODO: Compare more of the output files?

# Run diagnostic tool
# Paired end FASTQ
diagnosticTool_scripts/kodoja_search.py -r1 test/data/testData_1.fastq -o test/PE_test/ -d1 test/example_db/krakenDB_test/ -d2 test/example_db/kaijuDB_test/ -r2 ./test/data/testData_2.fastq -t 2
diff test/PE_test/virus_table.txt test/data/virus_table_PE_fastq.txt
rm -r test/PE_test/

# Single end FASTQ
diagnosticTool_scripts/kodoja_search.py -r1 test/data/testData_1.fastq -o test/SE_test/ -d1 test/example_db/krakenDB_test/ -d2 test/example_db/kaijuDB_test/ -t 2
diff test/SE_test/virus_table.txt test/data/virus_table_SE_fastq.txt
rm -r test/SE_test/

# Paired end FASTA
diagnosticTool_scripts/kodoja_search.py -r1 test/data/testData_1.fasta -o test/fasta_test/ -d1 test/example_db/krakenDB_test/ -d2 test/example_db/kaijuDB_test/ -r2 ./test/data/testData_2.fasta -f fasta
diff test/fasta_test/virus_table.txt test/data/virus_table_PE_fasta.txt
rm -r test/fasta_test/
