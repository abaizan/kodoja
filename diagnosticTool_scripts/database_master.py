import os
import time
from database_modules import *

#tool = "kraken"
tool = "kaiju"
genome_download_dir = "/home/ae42909/Scratch/ncbi_download/"
ncbi_download_parallel = "4"
threads = "10"
kraken_kmer = "31"
kraken_minimizer = "15"
# kraken_max_dbSize = 
kraken_db_dir =  "/home/ae42909/Scratch/krakenDB_viral"
jellyfish_hash_size = "11400M"
kaiju_db_dir =  "/home/ae42909/Scratch/kaijuDB_viral"


# Check directory paths have "/" at the end
genome_download_dir += check_path(genome_download_dir)
kraken_db_dir += check_path(kraken_db_dir)
kaiju_db_dir += check_path(kaiju_db_dir)

# Download NCBI data
# ncbi_download(tool, genome_download_dir, ncbi_download_parallel)

# Rename downloaded genomic files for Kraken or protein file for Kaiju
ncbi_rename_customDB(tool, genome_download_dir)

if tool == "kraken":
    # Make Kraken database
    krakenDB_build(genome_download_dir, kraken_db_dir, threads, kraken_kmer, kraken_minimizer, jellyfish_hash_size)
else:
    # Make Kaiju database
    kaijuDB_build(genome_download_dir, kaiju_db_dir)




