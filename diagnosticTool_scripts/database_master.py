import os
import time
from database_modules import *

threads = "10"
genome_download_dir = "/home/ae42909/Scratch/ncbi_downloads/"
kraken_kmer = "31"
kraken_minimizer = "15"
# kraken_max_dbSize = 
# kraken_db_dir =  "/home/ae42909/Scratch/krakenDB_viral"
kraken_db_dir =  "/home/ae42909/Scratch/new_krakenDB_viral"
jellyfish_hash_size = "11400M"


# Check that genome_download_dir and kraken_db_dir have "/" at the end
if genome_download_dir[-1] != "/":
    genome_download_dir += "/"

if kraken_db_dir[-1] != "/":
    kraken_db_dir += "/"


# Download NCBI data
# For RefSeq genomes file_format = "fasta", for RefSeq protein file_format = "protein-fasta".
# "parallel" referes to how many genomes to download in parallel (e.g. parallel = "4")
# ncbi_download("fasta", genome_download_dir, "4")
# ncbi_download("protein-fasta", genome_download_dir, "4")

# Rename downloaded genomic files for Kraken
# kraken_fna_rename(genome_download_dir)

# Make Kraken db
krakenDB_build(genome_download_dir, kraken_db_dir, threads, kraken_kmer, kraken_minimizer, jellyfish_hash_size)

# krakenDB_shrink(kraken_db_dir, "1000", "/home/ae42909/Scratch/minikraken_1000")





