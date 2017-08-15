import os
import time
from diagnostic_module import *

threads = "10"
genome_download_dir = "/home/ae42909/Scratch/ncbi_downloads"
kraken_kmer = "31"
kraken_minimizer = "15"
# kraken_max_dbSize = 
kraken_db_dir =  "/home/ae42909/Scratch/krakenDB_viral"

# Download NCBI data
# krakenDB_download(genome_download_dir, viral = True, plant = True, bacteria = False)

# Make Kraken db
krakenDB_build(genome_download_dir, kraken_db_dir, threads, kraken_kmer, kraken_minimizer)

# krakenDB_shrink(kraken_db_dir, "1000", "/home/ae42909/Scratch/minikraken_1000")





