import os
import time
import sys
sys.path.insert(0, '/home/ae42909/viral_diagnostics/diagnosticTool_scripts/')
from diagnostic_modules import *

# General parameters
file1 = "/home/ae42909/Scratch/100_Potato_withViruses_1.fastq"
file2 = "/home/ae42909/Scratch/100_Potato_withViruses_2.fastq"
user_format = "fastq"
out_dir = "/home/ae42909/Scratch/100Seq_krakenDB_viral/"
threads = 4
ncbi_file = '/home/ae42909/Scratch/kraken/kraken_analysis/customDatabase/NCBI_taxonomy.csv'
subset = True

# Trimmomatic paramters
trim_minlen = 50
adapter_file = "/mnt/apps/trimmomatic/IlluminaAdapters.fasta"

# Kraken parameters
kraken_db = "/home/ae42909/Scratch/krakenDB_k31_m15/"
quick_minhits = False
preload = False

# Kaiju parameters
kaiju_db = "/home/ae42909/Scratch/kaijuDB_viral_2/"
kaiju_minlen = 15
kaiju_score = False
kaiju_mismatch = False

# Check that dirs have "/" at the end
out_dir += check_path(out_dir)
kraken_db += check_path(kraken_db)
kaiju_db += check_path(kaiju_db)

# Check out_dir exits else make dir
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

os.chdir(out_dir)

t0 = time.time()

test_format(file1, user_format)
t1 = time.time()

if file2:
    # Check paired reads are paired correctly
    paired_test(file1, file2, user_format, out_dir)

    # Set tool files
    kraken_file1 = kaiju_file1 = "PE_trimmed_data_1P"
    kraken_file2 = kaiju_file2 = "PE_trimmed_data_2P"

else:
    kraken_file1 = kaiju_file1 = "SE_trimmed_data"
    kraken_file2 = kaiju_file2 = False

if user_format == "fastq":
    # fasta files cannot be QC'd - only for fastq files
    fastqc_trim(out_dir, file1, trim_minlen, threads, adapter_file, file2)
else:
    kraken_file1 = kaiju_file1 = file1
    kraken_file2 = kaiju_file2 = file2

if subset:
    kaiju_file1 = "subset_file1." + user_format
    if file2:
        kaiju_file2 = "subset_file2." + user_format


# Kraken classification
kraken_classify(kraken_file1, threads, user_format, kraken_db, kraken_file2,
                quick_minhits = quick_minhits, preload = preload)
t4 = time.time()

# Format kraken data and subset viral and unclassified sequences
seq_reanalysis("kraken_table.txt", "kraken_labels.txt", ncbi_file, out_dir,
               user_format, kraken_file1, subset, kraken_file2)
t5 = time.time()

# Kaiju classification of all sequences or subset sequences
kaiju_classify(kaiju_file1, threads, out_dir, kaiju_db, kaiju_minlen, kraken_db,
               kaiju_file2, kaiju_mismatch = kaiju_mismatch,
               kaiju_score = kaiju_score)
t6 = time.time()     

# Merege results
result_analysis(out_dir, "kraken_VRL.txt", "kaiju_table.txt", "kaiju_labels.txt",
                ncbi_file)
t7 = time.time()


# Create log file
if subset:
	print_statment = "subset sequences = " + str((t5-t4)/60) + " min\n"
else:
	print_statment = "formatting kraken data = " + str((t5-t4)/60) + " min\n"

log_file = open(out_dir + "log_file.txt", "w")
log_file.write("General parameters:\n" + "file1 = " + file1 + "\n" +
               "file2 = " + str(file2) + "\n" + "output directory = " +
               out_dir + "\n" + "threads =" + str(threads) + "\n" +
               "subset =" + str(subset) + "\n")
log_file.write("Trimmomatic parameters:\n" + "trim_minlen = " + str(trim_minlen) + "\n")
log_file.write("Kraken parameters:\n" + "kraken database = " + kraken_db + "\n" +
               "quick_minhits = " + str(quick_minhits) + "\n" + "preload =" +
               str(preload) + "\n")
log_file.write("Kaiju parameters:\n" + "kaiju_db =" + kaiju_db + "\n" +
               "kaiju_minlen = " + str(kaiju_minlen) + "\n" +
               "kaiju_score =" + str(kaiju_score) + "\n" +
               "kaiju_mismatch =" + str(kaiju_mismatch) + "\n")
log_file.write("Script timer:\n" + "testing format = " + str((t1-t0)) + " s\n" +
               "fastq and trim = " + str((t2-t1)/60) + " min\n" +
               "kraken classification = " + str((t4-t2)/3600) + " h\n" +
               print_statment + "kaiju classification = " + str((t6-t5)/3600) +
               " h\n" + "Results = " + str((t7-t6)/3600) + " h\n" + "total = " +
               str((t7-t0)/3600) + " h\n")
log_file.close()


