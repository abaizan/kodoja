import os
import time
import sys
from kraken_parameterValues import krakenDB_parameters

sys.path.insert(0, '/home/ae42909/viral_diagnostics/diagnosticTool_scripts/')
from database_modules import *

genome_download_dir = "/home/ae42909/Scratch/refDB_files/"
# ncbi_download_parallel = 4
threads = 4
# core_ram = 8 
#jellyfish_hash_size = str(int((threads * core_ram * 1000)/6.9)) + "M"

for kmers in krakenDB_parameters["kraken_kmer"]:
    for minimizers in krakenDB_parameters["kraken_minimizer"]:
        if kmers > minimizers:
            kraken_db_dir = "/home/ae42909/Scratch/parameter_test/kraken/databases/krakenDB_k" + str(kmers) + "_m" + str(minimizers)
            kraken_kmer = kmers
            kraken_minimizer = minimizers
            t0 = time.time()
            krakenDB_build(genome_download_dir, kraken_db_dir, threads, kraken_kmer, kraken_minimizer)
            t1 = time.time()
            with open(kraken_db_dir + '/' + 'log_file.txt', 'w') as out_file:
                text = "k-mer = " + str(kmers) + "\n"
                text += "minimizer = " + str(minimizers) + "\n"
                text += "threads = " + str(threads) + "\n"
                out_file.write(text)
