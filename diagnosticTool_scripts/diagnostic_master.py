import os
import time
from diagnostic_modules import *

# General parameters
file2 = False
# file1 = "/home/ae42909/synthetic_RNAseq/mappingRNAseq/concatenated_fastaFiles/Potato_withViruses_1.fastq"
# file2 = "/home/ae42909/synthetic_RNAseq/mappingRNAseq/concatenated_fastaFiles/Potato_withViruses_2.fastq"
# file1 = "/mnt/shared/projects/virology/201702_Stuart_local_raspvirus/static/Dee-B5_S1_R1.fastq"
# file2 = "/mnt/shared/projects/virology/201702_Stuart_local_raspvirus/static/Dee-B5_S1_R2.fastq"
file1 = '/home/ae42909/data_forTesting/Barerro_data/PB64-S7_clean.fq.gz_trim.fastq'
user_format = "fastq"
rna_type = "sRNA"
out_dir = "/home/ae42909/Scratch/Barrero_data_results/PB64-S7_1/"
threads = 4
subset = False

# Trimmomatic paramters
if rna_type == "RNA":
    trim_minlen = 50
elif rna_type == "sRNA":
    trim_minlen = 18

adapter_file = "/mnt/apps/trimmomatic/IlluminaAdapters.fasta"

# Kraken parameters
kraken_db = "/home/ae42909/Scratch/krakenDB_k18_m5/"
quick_minhits = False
preload = False

# Kaiju parameters
kaiju_db = "/home/ae42909/Scratch/kaijuDB_viral_2/"
kaiju_minlen = 3
kaiju_score = 45
kaiju_mismatch = 1

# Check that dirs have "/" at the end
out_dir += check_path(out_dir)
kraken_db += check_path(kraken_db)
kaiju_db += check_path(kaiju_db)

# Check out_dir exits else make dir
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

os.chdir(out_dir)

# Check if data is compressed - if so uncompress


t0 = time.time()

# Test data format - fasta or fastq
test_format(file1, user_format)
t1 = time.time()


if file2:    
    # strip metadata from ids (if any), assert the paired sequences have the same number of sequences and they are synchronised
    paired_test(file1, file2, user_format, out_dir)
    
    # fasta files cannot be QC'd - only for fastq files
    if user_format == "fastq":
        # QC and trim data
        fastqc_trim(out_dir, "renamed_1", trim_minlen, threads, adapter_file, "renamed_2")
        kraken_file1 = "PE_trimmed_data_1P"
        kraken_file2 = "PE_trimmed_data_2P"
        t2 = time.time()
    else:
        kraken_file1 = file1
        kraken_file2 = file2
        t2 = t1

    # Kraken classification
    kraken_classify(kraken_file1, threads, user_format, kraken_db, kraken_file2)
    t4 = time.time()

    # Format kraken data and subset viral and unclassified sequences
    seq_reanalysis("kraken_table.txt", "kraken_labels.txt", out_dir, user_format, "PE_trimmed_data_1P", subset, "PE_trimmed_data_2P")
    t5 = time.time()

    # Set data for kaiju analysis
    if subset:
        kaiju_file1 = "subset_file1." + user_format
        kaiju_file2 = "subset_file2." + user_format
    else:
        
        kaiju_file1 = "PE_trimmed_data_1P"
        kaiju_file2 = "PE_trimmed_data_2P"

    # Kaiju classification of all sequences or subset sequences
    kaiju_classify(kaiju_file1, threads, kaiju_db, kaiju_minlen, kraken_db, kaiju_file2, kaiju_mismatch, kaiju_score)
    t6 = time.time()     

    # Merege results
    result_analysis(out_dir, "kraken_VRL.txt", "kaiju_table.txt", "kaiju_labels.txt", "ID1.txt", "ID2.txt")
    t7 = time.time()
else:
    # fasta files cannot be QC'd - only for fastq files
    if user_format == "fastq":
        # QC and trim data
        fastqc_trim(out_dir, file1, trim_minlen, threads, adapter_file)
        kraken_file1 = "SE_trimmed_data"
        t2 = time.time()
    else:
        kraken_file1 = file1
        t2 = t1

    # Kraken classification
    kraken_classify(kraken_file1, threads, user_format, kraken_db,)
    t4 = time.time()

    # Subset viral and unclassified sequences
    seq_reanalysis("kraken_table.txt", "kraken_labels.txt", out_dir, user_format, "SE_trimmed_data", subset)
    t5 = time.time()

        # Set data for kaiju analysis
    if subset:
        kaiju_file1 = "subset_file1." + user_format
    else:
        
        kaiju_file1 = "SE_trimmed_data"

    # Kaiju classification of subset sequences
    kaiju_classify(kaiju_file1, threads, kaiju_db, kaiju_minlen, kraken_db, kaiju_file2 = False, kaiju_mismatch = kaiju_mismatch, kaiju_score = kaiju_mismatch)
    t6 = time.time()

    # Merege results
    result_analysis(out_dir, "kraken_VRL.txt", "kaiju_table.txt", "kaiju_labels.txt", "ID.txt")
    t7 = time.time()

# Creating a log file	
if subset:
	print_statment = "subset sequences = " + str((t5-t4)/60) + " min\n"
else:
	print_statment = "formatting kraken data = " + str((t5-t4)/60) + " min\n"

log_file = open(out_dir + "log_file.txt", "w")
log_file.write("General parameters:\n" + "file1 = " + file1 + "\n" + "file2 = " + str(file2) + "\n" + "output directory = " + out_dir + "\n" + "threads =" + str(threads) + "\n" + "subset =" + str(subset) + "\n")
log_file.write("Trimmomatic parameters:\n" + "trim_minlen = " + str(trim_minlen) + "\n")
log_file.write("Kraken parameters:\n" + "kraken database = " + kraken_db + "\n" + "quick_minhits = " + str(quick_minhits) + "\n" + "preload =" + str(preload) + "\n")
log_file.write("Kaiju parameters:\n" + "kaiju_db =" + kaiju_db + "\n" + "kaiju_minlen = " + str(kaiju_minlen) + "\n" + "kaiju_score =" + str(kaiju_score) + "\n" + "kaiju_mismatch =" + str(kaiju_mismatch) + "\n")
log_file.write("Script timer:\n" + "testing format = " + str((t1-t0)) + " s\n" + "fastq and trim = " + str((t2-t1)/60) + " min\n" + "kraken classification = " + str((t4-t2)/3600) + " h\n" + print_statment + "kaiju classification = " + str((t6-t5)/3600) + " h\n" + "Results = " + str((t7-t6)/3600) + " h\n" + "total = " + str((t7-t0)/3600) + " h\n")
log_file.close()


