#!/usr/bin/env python
import os
import time
from diagnostic_modules import check_path
from diagnostic_modules import test_format
from diagnostic_modules import check_file
from diagnostic_modules import fastqc_trim
from diagnostic_modules import kraken_classify
from diagnostic_modules import seq_reanalysis
from diagnostic_modules import kaiju_classify
from diagnostic_modules import result_analysis
import argparse

parser = argparse.ArgumentParser(description='Kodoja')
parser.add_argument('-r1', '--read1', type=str, required=True,
                    help='Read 1 file path, required')
parser.add_argument('-o', '--output_dir', type=str, required=True,
                    help='Output directory path, required')
parser.add_argument('-d1', '--kraken_db', type=str, required=True,
                    help='Kraken database path, required')
parser.add_argument('-d2', '--kaiju_db', type=str, required=True,
                    help='Kaiju database path, required')
parser.add_argument('-r2', '--read2', type=str, default=False,
                    help='Read 2 file path')
parser.add_argument('-f', '--data_format', type=str, default='fastq',
                    help='Sequence data format')
parser.add_argument('-t', '--threads', type=int, default=1,
                    help='Number of threads')
parser.add_argument('-s', '--host_subset', action='store_true',
                    help='Subset host sequences before Kaiju')
parser.add_argument('-m', '--trim_minlen', type=int, default=50,
                    help='Trimmomatic minimum length')
parser.add_argument('-a', '--trim_adapt', type=str, default=False,
                    help='Illumina adapter sequence file')
parser.add_argument('-q', '--kraken_quick', type=int, default=False,
                    help='Number of minium hits by Kraken')
parser.add_argument('-p', '--kraken_preload', action='store_true',
                    help='Kraken preload database')
parser.add_argument('-c', '--kaiju_score', type=int, default=85,
                    help='Kaju alignment score')
parser.add_argument('-l', '--kaiju_minlen', type=int, default=15,
                    help='Kaju minimum length')
parser.add_argument('-i', '--kaiju_mismatch', type=int, default=85,
                    help='Kaju allowed mismatches')
args = parser.parse_args()

# Check that dirs have "/" at the end
args.output_dir += check_path(args.output_dir)
args.kraken_db += check_path(args.kraken_db)
args.kaiju_db += check_path(args.kaiju_db)

# Check args.output_dir exits else make dir
if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

# os.chdir(args.output_dir)

# Write a log_file:
with open(args.output_dir + "log_file.txt", "w") as log_file:
    log_file.write("General parameters:\n" + "file1 = " + args.read1 + "\n" +
                   "file2 = " + str(args.read2) + "\n" +
                   "output directory = " +  args.output_dir + "\n" +
                   "threads =" + str(args.threads) + "\n" +
                   "host_subset =" + str(args.host_subset) + "\n" +
                   "Trimmomatic parameters:\n" + "trim_minlen = " + str(args.trim_minlen) + "\n" +
                   "Kraken parameters:\n" + "kraken database = " + args.kraken_db + "\n" +
                   "quick_minhits = " + str(args.kraken_quick) + "\n" + "preload =" +
                   str(args.kraken_preload) + "\n" +
                   "Kaiju parameters:\n" + "args.kaiju_db =" + args.kaiju_db + "\n" +
                   "kaiju_minlen = " + str(args.kaiju_minlen) + "\n" +
                   "kaiju_score = " + str(args.kaiju_score) + "\n" +
                   "kaiju_mismatch = " + str(args.kaiju_mismatch) + "\n")


t0 = time.time()
# Test format, change seqIDs and check paired files are correct
test_format(args.read1, args.data_format)
check_file(args.read1, args.output_dir, args.data_format, args.read2)
t1 = time.time()

# Set all variables
initial_file1 = args.output_dir + 'renamed_file_1.' + args.data_format
kraken_file1 = kaiju_file1 = args.output_dir + "trimmed_read1"

if args.read2:
    # Set tool files
    kraken_file2 = kaiju_file2 = args.output_dir + "trimmed_read2"
    initial_file2 = args.output_dir + 'renamed_file_2.' + args.data_format

else:
    kraken_file2 = kaiju_file2 = False
    initial_file2 = False

if args.data_format == "fastq":
    # fasta files cannot be QC'd - only for fastq files
    fastqc_trim(args.output_dir, initial_file1, args.trim_minlen, args.threads, args.trim_adapt, initial_file2)
else:
    kraken_file1 = kaiju_file1 = initial_file1
    kraken_file2 = kaiju_file2 = initial_file2


if args.host_subset:
    kaiju_file1 = args.output_dir + "subset_file1." + args.data_format
    if args.read2:
        kaiju_file2 = args.output_dir + "subset_file2." + args.data_format
t2 = time.time()

# Kraken classification
kraken_classify(args.output_dir, kraken_file1, args.threads, args.data_format, args.kraken_db, kraken_file2,
                quick_minhits = args.kraken_quick, preload = args.kraken_preload)
t3 = time.time()

# Format kraken data and subset unclassified and non-host sequences
seq_reanalysis("kraken_table.txt", "kraken_labels.txt", args.output_dir,
               args.data_format, kraken_file1, args.host_subset, kraken_file2)
t4 = time.time()

# Kaiju classification of all sequences or subset sequences
kaiju_classify(kaiju_file1, args.threads, args.output_dir, args.kaiju_db, args.kaiju_minlen, args.kraken_db,
               kaiju_file2, kaiju_mismatch = args.kaiju_mismatch,
               kaiju_score = args.kaiju_score)
t5 = time.time()     

# Merege results
result_analysis(args.output_dir, "kraken_VRL.txt", "kaiju_table.txt", "kaiju_labels.txt")
t6 = time.time()


# Create log file
if args.host_subset:
    print_statment = "subset sequences = " + str((t4-t3)/60) + " min\n"
else:
    print_statment = "formatting kraken data = " + str((t4-t3)/60) + " min\n"

with open(args.output_dir + "log_file.txt", "a") as log_file:
    log_file.write("Script timer:\n" + "testing format/replace seqID = " + str((t1-t0)) + " s\n" +
                   "fastq and trim = " + str((t2-t1)/60) + " min\n" +
                   "kraken classification = " + str((t3-t2)/3600) + " h\n" +
                   print_statment + "kaiju classification = " + str((t5-t6)/3600) +
                   " h\n" + "Results = " + str((t6-t5)/3600) + " h\n" + "total = " +
                   str((t6-t0)/3600) + " h\n")
