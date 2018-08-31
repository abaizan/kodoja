#!/bin/bash

#$ -cwd
#$ -j yes
#$ -pe smp 2

# Enable strict bash mode - we want to abort on the first error, see also
# http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail

if [ ! -f "test/test_script.sh" ]
then
    echo "ERROR. Run this from the GitHub repository root directory."
    echo "Please don't run from the test directory itself."
    exit 1
fi

# Confirm taxonomy files present
if [ ! -f "test/taxonomy/nodes.dmp" ] || [ ! -f "test/taxonomy/names.dmp" ]
then
    echo "ERROR. Please cache the NCBI taxonomy files before running tests."
    echo "Please consult instuctions in file test/taxonomy/README.txt"
    exit 1
fi

echo "Which version of Python do we have?"
python --version

# Confirm binary dependencies present (strict bash mode will abort on failure)
echo "Do we have ncbi-genome-download?"
ncbi-genome-download --version

echo "Do we have kraken?"
kraken-build --version

echo "Do we have kaiju?"
# Hiding kaiju call in a subprocess because it will give non-zero return code
# First grep is to get just one line of output; second is to fail if missing
echo `kaiju 2>&1 | grep Kaiju` | grep Kaiju

echo "Do we have trimmomatic?"
trimmomatic -version

echo "Do we have fastqc?"
fastqc --version

echo "Do we have numpy?"
python -c "from numpy import __version__; print(__version__)"

echo "Do we have pandas?"
python -c "from pandas import __version__; print(__version__)"

echo "Do we have biopython?"
python -c "from Bio import __version__; print(__version__)"

# Catches any syntax errors, or broken version option:
echo "What version of kodoja build?"
./diagnosticTool_scripts/kodoja_build.py --version
echo "What version of kodoja search?"
./diagnosticTool_scripts/kodoja_search.py --version
echo "What version of kodoja retrieve?"
./diagnosticTool_scripts/kodoja_retrieve.py --version

# Generating test data using modern Illumina paired naming without /1 and /2
# Also deliberately using .fq rather than .fastq to ensure that works too
if [ ! -f "test/data/test_names_R1.fq" ]
then
    echo "Generating test_names_R1.fq"
    sed -e 's#/1$# 1:N:0:CTCCGC#g' test/data/testData_1.fastq > test/data/test_names_R1.fq
fi
if [ ! -f "test/data/test_names_R2.fq" ]
then
    echo "Generating test_names_R2.fq"
    sed -e 's#/2$# 2:N:0:CTCCGC#g' test/data/testData_2.fastq > test/data/test_names_R2.fq
fi
if [ ! -f "test/data/filtered_PE_fastq_names_stringent_R1.fq" ]
then
    echo "Generating filtered_PE_fastq_names_stringent_R1.fq"
    sed -e 's#/1$# 1:N:0:CTCCGC#g' test/data/filtered_PE_fastq_stringent_R1.fastq > test/data/filtered_PE_fastq_names_stringent_R1.fq
fi
if [ ! -f "test/data/filtered_PE_fastq_names_stringent_R2.fq" ]
then
    echo "Generating filtered_PE_fastq_names_stringent_R2.fq"
    sed -e 's#/2$# 2:N:0:CTCCGC#g' test/data/filtered_PE_fastq_stringent_R2.fastq > test/data/filtered_PE_fastq_names_stringent_R2.fq
fi

echo "Beginning tests..."

echo "=============================================================="
echo "Testing kodoja_build.py"
echo "=============================================================="
# Create test database - doing this from a clean slate ought to work,
# but we run out of disk space on TravisCI while downloading and
# unzipping the taxonomy data. Instead we provide a prepopulated
# minimal test/building_db/krakenDB_test/taxonomy/ folder.

# Download three specific viruses for a mini test database:
export TMP=${TMP:-/tmp}

if [ ! -f "137758.fna.gz" ]
then
    ncbi-genome-download --verbose -o $TMP/ -F fasta -t 137758 viral
    ln -f -s $TMP/refseq/viral/GCF_000884835.1/GCF_000884835.1_ViralProj38085_genomic.fna.gz 137758.fna.gz
fi
if [ ! -f "137758.faa.gz" ]
then
    ncbi-genome-download --verbose -o $TMP/ -F protein-fasta -t 137758 viral
    ln -f -s $TMP/refseq/viral/GCF_000884835.1/GCF_000884835.1_ViralProj38085_protein.faa.gz 137758.faa.gz
fi
if [ ! -f "946046.fna.gz" ]
then
    ncbi-genome-download --verbose -o $TMP/ -F fasta -t 946046 viral
    ln -f -s $TMP/refseq/viral/GCF_000888855.1/GCF_000888855.1_ViralProj61097_genomic.fna.gz 946046.fna.gz
fi
if [ ! -f "946046.faa.gz" ]
then
    ncbi-genome-download --verbose -o $TMP/ -F protein-fasta -t 946046 viral
    ln -f -s $TMP/refseq/viral/GCF_000888855.1/GCF_000888855.1_ViralProj61097_protein.faa.gz 946046.faa.gz
fi
if [ ! -f "12227.fna.gz" ]
then
    ncbi-genome-download --verbose -o $TMP/ -F fasta -t 12227 viral
    ln -f -s $TMP/refseq/viral/GCF_000861345.1/GCF_000861345.1_ViralProj15325_genomic.fna.gz 12227.fna.gz
fi
if [ ! -f "12227.faa.gz" ]
then
    ncbi-genome-download --verbose -o $TMP/ -F protein-fasta -t 12227 viral
    ln -f -s $TMP/refseq/viral/GCF_000861345.1/GCF_000861345.1_ViralProj15325_protein.faa.gz 12227.faa.gz
fi

echo "Running kodoja_build.py with three viruses as input"
diagnosticTool_scripts/kodoja_build.py -o test/building_db/ -t 1 --db_tag 'test' -v -n -e 137758.fna.gz 946046.fna.gz 12227.fna.gz 137758.faa.gz 946046.faa.gz 12227.faa.gz -x 137758 946046 12227 137758 946046 12227 -k 18 -m 5
ls test/building_db
diff test/example_db/krakenDB_test/database.idx test/building_db/krakenDB_test/database.idx
diff test/example_db/krakenDB_test/database.kdb test/building_db/krakenDB_test/database.kdb
#Binary files differ slightly - although same size
#diff test/example_db/kaijuDB_test/kaiju_library.fmi test/building_db/kaijuDB_test/kaiju_library.fmi
if [ ! -s test/building_db/kaijuDB_test/kaiju_library.fmi ]; then echo "Missing kaiju_library.fmi" && false; fi

if [ -d test/building_db/extra/ ]
then
    echo "Removing test/building_db/extra/ from a previous test run"
    rm -rf test/building_db/extra/
fi

echo "=============================================================="
echo "Testing kodoja_search.py with paired end FASTQ with /1, /2"
echo "=============================================================="
diagnosticTool_scripts/kodoja_search.py -r1 test/data/testData_1.fastq -o test/PE_test/ -d1 test/example_db/krakenDB_test/ -d2 test/example_db/kaijuDB_test/ -r2 ./test/data/testData_2.fastq -t 2
diff test/PE_test/virus_table.txt test/data/virus_table_PE_fastq.txt
diff test/PE_test/kodoja_VRL.txt test/data/virus_reads_PE_fastq.txt

echo "=============================================================="
echo "Testing kodoja_retrieve.py with paired end FASTQ with /1, /2"
echo "=============================================================="
# Using -s for stringent, defaults to fastq
diagnosticTool_scripts/kodoja_retrieve.py -o test/PE_test/ -r1 test/data/testData_1.fastq -r2 test/data/testData_2.fastq -s
ls test/PE_test/subset_files/
diff test/PE_test/subset_files/virus_all_sequences1.fastq test/data/filtered_PE_fastq_stringent_R1.fastq
diff test/PE_test/subset_files/virus_all_sequences2.fastq test/data/filtered_PE_fastq_stringent_R2.fastq
rm -r test/PE_test/

echo "=============================================================="
echo "Testing kodoja_search.py with paired end FASTQ without /1, /2"
echo "=============================================================="
diagnosticTool_scripts/kodoja_search.py -r1 test/data/test_names_R1.fq -r2 test/data/test_names_R2.fq -o test/PE_test/ -d1 test/example_db/krakenDB_test/ -d2 test/example_db/kaijuDB_test -t 2
diff test/PE_test/virus_table.txt test/data/virus_table_PE_fastq.txt
diff test/PE_test/kodoja_VRL.txt test/data/virus_reads_PE_names_fastq.txt

echo "=============================================================="
echo "Testing kodoja_retrieve.py on paired end FASTQ without /1, /2"
echo "=============================================================="
# Using -s for stringent, defaults to fastq format and extension
diagnosticTool_scripts/kodoja_retrieve.py -o test/PE_test/ -r1 test/data/test_names_R1.fq -r2 test/data/test_names_R2.fq -s
ls test/PE_test/subset_files/
diff test/PE_test/subset_files/virus_all_sequences1.fastq test/data/filtered_PE_fastq_names_stringent_R1.fq
diff test/PE_test/subset_files/virus_all_sequences2.fastq test/data/filtered_PE_fastq_names_stringent_R2.fq
rm -r test/PE_test/

echo "=============================================================="
echo "Testing kodoja_search.py with single end FASTQ"
echo "=============================================================="
diagnosticTool_scripts/kodoja_search.py -r1 test/data/testData_1.fastq -o test/SE_test/ -d1 test/example_db/krakenDB_test/ -d2 test/example_db/kaijuDB_test/ -t 2
diff test/SE_test/virus_table.txt test/data/virus_table_SE_fastq.txt
diff test/SE_test/kodoja_VRL.txt test/data/virus_reads_SE_fastq.txt

echo "=============================================================="
echo "Testing kodoja_retrieve.py with single end FASTQ"
echo "=============================================================="
diagnosticTool_scripts/kodoja_retrieve.py -o test/SE_test/ -r1 test/data/testData_1.fastq -g
ls test/SE_test/subset_files/
diff test/SE_test/subset_files/virus_all_sequences1.fastq test/data/filtered_SE_fastq_genus.fastq
# Now again but with taxid option
rm -rf test/SE_test/subset_files/
diagnosticTool_scripts/kodoja_retrieve.py -o test/SE_test/ -r1 test/data/testData_1.fastq -t 12227
diff test/SE_test/subset_files/virus_12227_sequences1.fastq test/data/filtered_SE_fastq_12227.fastq
rm -r test/SE_test/

echo "=============================================================="
echo "Testing kodoja_search.py with paired end FASTA"
echo "=============================================================="
diagnosticTool_scripts/kodoja_search.py -r1 test/data/testData_1.fasta -o test/fasta_test/ -d1 test/example_db/krakenDB_test/ -d2 test/example_db/kaijuDB_test/ -r2 ./test/data/testData_2.fasta -f fasta
diff test/fasta_test/virus_table.txt test/data/virus_table_PE_fasta.txt
diff test/fasta_test/kodoja_VRL.txt test/data/virus_reads_PE_fasta.txt

echo "=============================================================="
echo "Testing kodoja_retreive.py with paired end FASTA"
echo "=============================================================="
# Using -g for genus level, explicitly saying using fasta
diagnosticTool_scripts/kodoja_retrieve.py -o test/fasta_test/ -r1 test/data/testData_1.fasta -r2 test/data/testData_2.fasta -f fasta -g
ls test/fasta_test/subset_files/
diff test/fasta_test/subset_files/virus_all_sequences1.fasta test/data/filtered_PE_fasta_genus_R1.fasta
diff test/fasta_test/subset_files/virus_all_sequences2.fasta test/data/filtered_PE_fasta_genus_R2.fasta
rm -r test/fasta_test/

echo "=============================================================="
echo "Testing finished."
