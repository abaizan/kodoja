import os
import time
from database_modules import *
from kraken_parameterValues import krakenDB_parameters


genome_download_dir = "/home/ae42909/Scratch/ncbi_download/"
ncbi_download_parallel = 4
threads = 4
core_ram = 8 
#jellyfish_hash_size = str(int((threads * core_ram * 1000)/6.9)) + "M"

for kmers in krakenDB_parameters["kraken_kmer"]:
    for minimizers in krakenDB_parameters["kraken_minimizer"]:
        if kmers > minimizers:
            kraken_db_dir = "/home/ae42909/Scratch/parameter_test/kraken/krakenDB_k" + str(kmers) + "_m" + str(minimizers)
            kraken_kmer = kmers
            kraken_minimizer = minimizers
            t0 = time.time()
            krakenDB_build(genome_download_dir, kraken_db_dir, threads, kraken_kmer, kraken_minimizer)
            t1 = time.time()
            print "krakenDB_k" + str(kmers) + "_m" + str(minimizers) + " = " + str(t1-t0)
