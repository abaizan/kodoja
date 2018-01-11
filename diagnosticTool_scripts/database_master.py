#!/usr/bin/env python
"""Script for running Kodoja database construction modules."""
import urllib
import os
import pandas as pd
import argparse
from diagnostic_modules import check_path
from database_modules import ncbi_download
from database_modules import ncbi_rename_customDB
from database_modules import krakenDB_build
from database_modules import kaijuDB_build

parser = argparse.ArgumentParser(description='Kodoja database construction')
parser.add_argument('-o', '--output_dir', type=str, required=True,
                    help='Output directory path, required')
parser.add_argument('-t', '--threads', type=int, default=1,
                    help='Number of threads, default=1')
parser.add_argument('-q', '--test', action='store_true',
                    help='Make test database')
parser.add_argument('-p', '--host', nargs='+', default=False,
                    help='Host tax ID')
parser.add_argument('-d', '--download_parallel', type=int, default=4,
                    help='Parallel genome download, default=4')
parser.add_argument('-n', '--no_download', action='store_false',
                    help='Genomes have already been downloaded')
parser.add_argument('-e', '--extra_files', type=str, nargs='*',
                    help='List of extra files added to "extra" dir')
parser.add_argument('-x', '--extra_taxids', type=str, nargs='*',
                    help='List of taxID of extra files')
parser.add_argument('-v', '--all_viruses', action='store_true',
                    help='Build databases with all viruses (not plant specific)')
parser.add_argument('-b', '--kraken_tax', type=str, default=False,
                    help='Path to taxonomy directory')
parser.add_argument('-k', '--kraken_kmer', type=int, default=31,
                    help='Kraken kmer size, default=31')
parser.add_argument('-m', '--kraken_minimizer', type=int, default=15,
                    help='Kraken minimizer size, default=15')
parser.add_argument('-a', '--db_tag', type=str, default=False,
                    help='Suffix for databases')
args = parser.parse_args()

# extra_files = ["Rubus_occidentalis_v1.0.a1.scaffolds.fna.gz",]
# extra_taxid = [75079,]

tool_list = ['kraken', 'kaiju']

args.output_dir += check_path(args.output_dir)
kraken_db_dir = args.output_dir + "krakenDB"
kaiju_db_dir = args.output_dir + "kaijuDB"

# Check args.output_dir exits, else make dir
if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

# Name databases with tag
if args.db_tag:
    kraken_db_dir += '_' + args.db_tag
    kaiju_db_dir += '_' + args.db_tag

kraken_db_dir += check_path(kraken_db_dir)
kaiju_db_dir += check_path(kaiju_db_dir)

# Ensure extra files have taxIDs and are in the right format
if args.extra_files:
    # Check there are the same number of files as taxids
    assert len(args.extra_files) == len(args.extra_taxids), \
        "Each extra file provided needs to have a corresponding ncbi taxid"

    # Check extra files provided are compressed and have the right file extension
    for files in args.extra_files:
        assert files.endswith(".fna.gz" or ".faa.gz"), \
            "File extensions need to be either '.fna' for genomic data" + \
            " or '.faa' for protein data, and be compressed ('.gz')"

# Download virus assembly summary for refseq
if not os.path.exists(args.output_dir + "viral_assembly_summary.txt"):
    urllib.urlretrieve('ftp://ftp.ncbi.nih.gov/genomes/refseq/viral/assembly_summary.txt',
                       args.output_dir + 'viral_assembly_summary.txt')
path_assembly_summary = args.output_dir + "viral_assembly_summary.txt"
vir_assembly = pd.read_table(path_assembly_summary, sep='\t', skiprows=1, header=0)
vir_assembly = vir_assembly.rename(columns={'# assembly_accession': 'assembly_accession'})

# Set subset_vir_assembly and vir_host
# subset_vir_assembly - list of virus accession names which will be added to databases (used in krakenDB_build and kaijuDB_build)
# vir_host - list of viral taxIDs for plant viruses. A subset of these are in refseq, be added to datbase by setting 'subset_vir_assembly'.
if args.all_viruses:
    subset_vir_assembly = False
    vir_host = False
else:
    if args.test:
        vir_host = args.test = [137758, 946046, 12227]
        args.kraken_kmer = 18
        args.kraken_minimizer = 5
    else:
        if not os.path.exists(args.output_dir + "virushostdb.tsv"):
            os.chdir(args.output_dir)
            urllib.urlretrieve('ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv',
                               'virushostdb.tsv')
        host_list = pd.read_csv(args.output_dir + "virushostdb.tsv", sep="\t").fillna('')
        plnVir = host_list[host_list['host lineage'].str.contains("Viridiplantae")]
        vir_host = list(plnVir['virus tax id'])

    subset_vir_assembly = list(vir_assembly.assembly_accession[vir_assembly['taxid'].isin(vir_host)])


for tool in tool_list:
    if args.no_download:
        # Download NCBI data
        ncbi_download(tool, args.output_dir, args.download_parallel, args.host, args.test)
        print "DONE with downloading"
    # Rename downloaded genomic files for Kraken or protein file for Kaiju
    ncbi_rename_customDB(tool, args.output_dir, args.extra_files, args.extra_taxids)
    print "DONE with renaming"
    # Make Kraken database
    if tool == "kraken":
        krakenDB_build(args.output_dir, kraken_db_dir, args.threads, args.kraken_kmer,
                       args.kraken_minimizer, subset_vir_assembly, args.kraken_tax)
        print "DONE with kraken db"
        with open(kraken_db_dir + "log_file.txt", "w") as out_file:
            text = 'output_dir = ' + args.output_dir + '\n'
            text += 'kraken_kmer = ' + str(args.kraken_kmer) + '\n'
            text += 'kraken_minimizer = ' + str(args.kraken_minimizer) + '\n'
            text += 'Genomes added to db = ' + str(subset_vir_assembly) + '\n'
            out_file.write(text)
    elif tool == "kaiju":
        # Make Kaiju database
        kaijuDB_build(args.output_dir, kaiju_db_dir, subset_vir_assembly)
        with open(kaiju_db_dir + "log_file.txt", "w") as out_file:
            text = 'output_dir = ' + args.output_dir + '\n'
            text += 'Genomes added to db = ' + str(subset_vir_assembly) + '\n'
            out_file.write(text)
