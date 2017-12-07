import urllib
import os
import pandas as pd
# import argparse
from diagnostic_modules import check_path
from database_modules import ncbi_download 
from database_modules import ncbi_rename_customDB
from database_modules import krakenDB_build
from database_modules import kaijuDB_build

# parser = argparse.ArgumentParser(description='Kodoja database construction')
# group = parser.add_mutually_exclusive_group()
# group.add_argument('kraken', action='store_true',
#                     help='')
# group.add_argument("kaiju", action="store_true")
# parser.add_argument('output_dir', type=str,
#                     help='Output directory path')
# parser.add_argument('-t', '--threads', type=int, default=1,
#                     help='Number of threads')
# parser.add_argument('-q', '--test', action='store_true',
#                     help='Make test database')

# args = parser.parse_args()

genome_download_dir = "/home/ae42909/Scratch/refDB_files/"
genome_download_dir += check_path(genome_download_dir)
threads = 4
tool = "kaiju"
download_files = False
extra_files = False
# extra_files = ["Rubus_occidentalis_v1.0.a1.scaffolds.fna.gz",]
extra_taxid = [75079,]
vir_host = 'plants'
# either specific taxIDs or 'plants' (plant viruses) or False for all viruses regardless of host
test = False

if test:
    vir_host = [137758, 946046, 12227]
    download_files = True

if tool == "kraken":
    kraken_kmer = 31
    kraken_minimizer = 15
    kraken_db_dir = "/home/ae42909/Scratch/diagnostic_databases/krakenDB_k31_plnVir/"
    kraken_db_dir += check_path(kraken_db_dir)
    # core_ram = 8
    # jellyfish_hash_size = str(int((threads * core_ram * 1000)/6.9)) + "M"
elif tool == "kaiju":
    kaiju_db_dir = "/home/ae42909/Scratch/diagnostic_databases/kaijuDB_plnVir"
    kaiju_db_dir += check_path(kaiju_db_dir)

if extra_files:
    # Check there are the same number of files as taxids
    assert len(extra_files) == len(extra_taxid), \
        "Each extra file provided needs to have a correspondinf ncbi taxid"

    # Check extra files provided are compressed and have the right file extension
    for files in extra_files:
        assert files.endswith(".fna.gz" or ".faa.gz"), \
            "File extensions need to be either '.fna' for genomic data" + \
            " or '.faa' for protein data, and be compressed ('.gz')"

# Download virus/host table and assembly summary for refseq
if not os.path.exists(genome_download_dir + "refseq/viral_assembly_summary.txt"):
    os.chdir(genome_download_dir)
    urllib.urlretrieve('ftp://ftp.ncbi.nih.gov/genomes/refseq/viral/assembly_summary.txt', 'refseq/viral_assembly_summary.txt')
path_assembly_summary = genome_download_dir + "refseq/viral_assembly_summary.txt"
vir_assembly = pd.read_table(path_assembly_summary, sep='\t', skiprows=1, header = 0)
vir_assembly = vir_assembly.rename(columns = {'# assembly_accession':'assembly_accession'})

if vir_host == 'plants':
    if not os.path.exists(genome_download_dir + "virushostdb.tsv"):
        os.chdir(genome_download_dir)
        urllib.urlretrieve('ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv', 'virushostdb.tsv')
    host_list = pd.read_csv(genome_download_dir + "virushostdb.tsv", sep="\t").fillna('')
    plnVir = host_list[host_list['host lineage'].str.contains("Viridiplantae")]
    vir_host = list(plnVir['virus tax id'])
    plnVir_assembly = list(vir_assembly.assembly_accession[vir_assembly['taxid'].isin(vir_host)])
elif not vir_host:
    plnVir_assembly = False

if download_files:
    ncbi_download_parallel = 4
    # Download NCBI data
    ncbi_download(tool, genome_download_dir, ncbi_download_parallel)

    # Rename downloaded genomic files for Kraken or protein file for Kaiju
    ncbi_rename_customDB(tool, genome_download_dir, extra_files, extra_taxid)

if tool == "kraken":
    # Make Kraken database
    krakenDB_build(genome_download_dir, kraken_db_dir, threads, kraken_kmer, kraken_minimizer, plnVir_assembly)
    with open(kraken_db_dir + "log_file.txt", "w") as out_file:
        text = 'genome_download_dir = ' + genome_download_dir + '\n'
        text += 'kraken_kmer = ' + str(kraken_kmer) + '\n'
        text += 'kraken_minimizer = ' + str(kraken_minimizer) + '\n'
        out_file.write(text)
elif tool == "kaiju":
    # Make Kaiju database
    kaijuDB_build(genome_download_dir, kaiju_db_dir, plnVir_assembly)
    with open(kaiju_db_dir + "log_file.txt", "w") as out_file:
        text = 'genome_download_dir = ' + genome_download_dir + '\n'
        out_file.write(text)
