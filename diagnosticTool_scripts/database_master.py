import os
import time
from database_modules import *

genome_download_dir = "/home/ae42909/Scratch/refDB_files/"
genome_download_dir += check_path(genome_download_dir)
threads = 4

tool = "kraken"
if tool == "kraken":
    kraken_kmer = 5
    kraken_minimizer = 2
    kraken_db_dir =  "/home/ae42909/Scratch/parameter_test/kraken/krakenDB_k5_m2"
    kraken_db_dir += check_path(kraken_db_dir)
    # core_ram = 8
    # jellyfish_hash_size = str(int((threads * core_ram * 1000)/6.9)) + "M"
elif tool == "kaiju":
    kaiju_db_dir =  "/home/ae42909/Scratch/kaijuDB_viral_2"
    kaiju_db_dir += check_path(kaiju_db_dir)

download_files = False
if download_files:
    ncbi_download_parallel = 4
    
extra_files = False
if extra_files:
    extra_files = ["Rubus_occidentalis_v1.0.a1.scaffolds.fna.gz",]
    extra_taxid = [75079,]


if download_files:
    if extra_files:
        # Check there are the same number of files as taxids
        assert len(extra_files) == len(extra_taxid), "Each extra file provided needs to have a correspondinf ncbi taxid"

        # Check extra files provided are compressed and have the right file extension
        for files in extra_files:
            assert files.endswith(".fna.gz" or ".faa.gz"),"File extensions need to be either '.fna' for genomic data or '.faa' for protein data, and be compressed ('.gz')"

    # Download NCBI data
    ncbi_download(tool, genome_download_dir, ncbi_download_parallel)

    # Rename downloaded genomic files for Kraken or protein file for Kaiju
    ncbi_rename_customDB(tool, genome_download_dir, extra_files, extra_taxid)

if tool == "kraken":
    # Make Kraken database
    krakenDB_build(genome_download_dir, kraken_db_dir, threads, kraken_kmer, kraken_minimizer)
    with open(kraken_db_dir + "log_file.txt", "w") as out_file:
        text = 'genome_download_dir = ' + genome_download_dir + '\n'
        text += 'kraken_kmer = ' + str(kraken_kmer) + '\n'
        text += 'kraken_minimizer = ' + str(kraken_minimizer) + '\n'
        out_file.write(text)
else:
    # Make Kaiju database
    kaijuDB_build(genome_download_dir, kaiju_db_dir)
    with open(kaiju_db_dir + "log_file.txt", "w") as out_file:
        text = 'genome_download_dir = ' + genome_download_dir + '\n'
        out_file.write(text)





