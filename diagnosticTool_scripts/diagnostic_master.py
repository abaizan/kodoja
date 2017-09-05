import os
import time
from diagnostic_modules import *

# General parameters
# file1 = "/home/ae42909/synthetic_RNAseq/mappingRNAseq/concatenated_fastaFiles/Potato_withViruses_1.fastq"
# file2 = "/home/ae42909/synthetic_RNAseq/mappingRNAseq/concatenated_fastaFiles/Potato_withViruses_2.fastq"
# file1 = "/home/ae42909/Scratch/100_Potato_withViruses_1.fastq"
# file2 = "/home/ae42909/Scratch/100_Potato_withViruses_2.fastq"
# file1 = "/home/ae42909/Scratch/100_Potato_withViruses_1.fasta"
# file2 = "/home/ae42909/Scratch/100_Potato_withViruses_2.fasta"
file1 = "/mnt/shared/projects/virology/201702_Stuart_local_raspvirus/static/Dee-B5_S1_R1.fastq"
file2 = "/mnt/shared/projects/virology/201702_Stuart_local_raspvirus/static/Dee-B5_S1_R2.fastq"
user_format = "fastq"
# out_dir = "/home/ae42909/Scratch/synthPotato_krakenDB_viral"
out_dir = "/home/ae42909/Scratch/berry/Dee-B5_S1"
user_format = "fastq"
threads = 4
# Trimmomatic paramters
trim_minlen = 50
adapter_file = "/mnt/apps/trimmomatic/IlluminaAdapters.fasta"
# Kraken parameters
kraken_db = "/home/ae42909/Scratch/krakenDB_viral/"
# quick_minhits = 2
# Kaiju parameters
kaiju_db = "/home/ae42909/Scratch/kaijuDB_viral/"
# kaiju_nodes = "/mnt/shared/scratch/ae42909/201609_BBSRC_Diagnostics/kaiju/kaijudb/nodes.dmp"
# kaiju_names = "/mnt/shared/scratch/ae42909/201609_BBSRC_Diagnostics/kaiju/kaijudb/names.dmp"
# kaiju_fmi = "/mnt/shared/scratch/ae42909/201609_BBSRC_Diagnostics/kaiju/kaijudb/kaiju_db.fmi"
kaiju_minlen = 11
kaiju_score = 65
kaiju_mismatch = 5

# Check that dirs have "/" at the end
out_dir += check_path(out_dir)
kraken_db += check_path(kraken_db)
kaiju_db += check_path(kaiju_db)

# Check out_dir exits else make dir
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

os.chdir(out_dir)

t0 = time.time()

# Test data format - fasta or fastq
test_format(file1, user_format)
t1 = time.time()


if file2:
    # fasta files cannot be QC'd - only for fastq files
    if user_format == "fastq":
        # QC and trim data
        fastqc_trim(out_dir, file1, trim_minlen, threads, adapter_file, file2)
        t2 = time.time()
        # Order and replace names from trimmed sequences
        rename_seq("PE_trimmed_data_1P", out_dir, user_format, paired = "1")
        rename_seq("PE_trimmed_data_2P", out_dir, user_format, paired = "2")
        t3 = time.time()
        
    elif user_format == "fasta":
        rename_seq(file1, out_dir, user_format, paired = "1")
        rename_seq(file2, out_dir, user_format, paired = "2")

    # Kraken classification
    kraken_classify("renamed_file1.", threads, user_format, kraken_db, "renamed_file2.")
    t4 = time.time()

    # Subset viral and unclassified sequences
    seq_reanalysis("kraken_table.txt", "kraken_labels.txt", out_dir, user_format, "renamed_file1.", "renamed_file2.")
    t5 = time.time()

    # Kaiju classification of subset sequences
    kaiju_classify("subset_file1.", user_format, threads, kaiju_db, kaiju_mismatch, kraken_db, "subset_file2.")
    t6 = time.time()     

    # Merege results
    result_analysis(out_dir, "kraken_FormattedTable.txt", "kaiju_table.txt", "kaiju_labels.txt", "ID1.txt", "ID2.txt")
    t7 = time.time()
else:
    # fasta files cannot be QC'd - only for fastq files
    if user_format == "fastq":
        # QC and trim data
        fastqc_trim(out_dir, file1, trimlen, threads)
        t2 = time.time()

        # Order and replace names from trimmed sequences
        rename_seq("SE_trimmed_data", out_dir, user_format, paired = False)
        t3 = time.time()
        
    elif user_format == "fasta":
        rename_seq(file1, out_dir, user_format, paired = False)
        t3 = time.time()

    # Kraken classification
    kraken_classify("renamed_file.", threads, user_format, kraken_db,)
    t4 = time.time()

    # Subset viral and unclassified sequences
    seq_reanalysis("kraken_table.txt", "kraken_labels.txt", out_dir, user_format, "renamed_file.")
    t5 = time.time()

    # Kaiju classification of subset sequences
    kaiju_classify("subset_file1.", user_format, threads, kaiju_db, kaiju_mismatch, kraken_db)
    t6 = time.time()

    # Merege results
    result_analysis(out_dir, "kraken_FormattedTable.txt", "kaiju_table.txt", "kaiju_labels.txt", "ID.txt")
    t7 = time.time()

print "testing format = " + str((t1-t0)/3600)
print "fastq and trim = " + str((t2-t1)/3600)
print "rename sequences = " + str((t3-t2)/3600)
print "kraken classification = " + str((t4-t3)/3600)
print "subset vrl and unclas sequences = " + str((t5-t4)/3600)
print "kaiju classification = " + str((t6-t5)/3600)
print "Results = " + str((t7-t6)/3600)
print "total = " + str((t7-t0)/3600)
