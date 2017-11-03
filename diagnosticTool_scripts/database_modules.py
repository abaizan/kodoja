# python 2.7.13
import subprocess
import pandas as pd
import numpy as np
import os
import urllib
import re

# Check directory paths have "/" at end
def check_path(dirs):
    if dirs[-1] != "/":
        return "/"
    else:
        return ""

# Download refseq files from ncbi ftp site - use ncbi-genome-download
def ncbi_download(tool, genome_download_dir, parallel=False):
    assert (tool == "kraken") | (tool == "kaiju"), "Argument 'tool' must be either 'kraken' or 'kaiju'."
    if tool == "kraken":
        file_format = "fasta"
    else:
        file_format = "protein-fasta"

    # Check directory exists
    if not os.path.exists(genome_download_dir):
        os.makedirs(genome_download_dir)
        
    ngd_command = "ncbi-genome-download -F " + file_format +  " -o " + genome_download_dir + " viral "
    if parallel:
        ngd_command += "--parallel " + str(parallel)
    subprocess.call(ngd_command, shell = True)
    
# Rename ncbi data files for custom databases
def ncbi_rename_customDB(tool, genome_download_dir, extra_files = False, extra_taxid = False):
    assert (tool == "kraken") | (tool == "kaiju"), "Argument 'tool' must be either 'kraken' or 'kaiju'."
    if tool == "kraken":
        file_extension = ".fna.gz"
    else:
        file_extension = ".faa.gz"

    # Download assembly summary file from ncbi ftp site and load as pandas dataframe
    if not os.path.exists(genome_download_dir + "refseq/viral_assembly_summary.txt"):
        os.chdir(genome_download_dir)
        urllib.urlretrieve('ftp://ftp.ncbi.nih.gov/genomes/refseq/viral/assembly_summary.txt', 'refseq/viral_assembly_summary.txt')
    path_assembly_summary = genome_download_dir + "refseq/viral_assembly_summary.txt"
    assembly_summary = pd.read_table(path_assembly_summary, sep='\t', skiprows=1, header = 0)
    assembly_summary.rename(columns={'# assembly_accession':'assembly_accession'}, inplace=True)  # rename column to exclude "#"

    # Count for protein sequences
    kaiju_count = 1
    for root, subdirs, files in os.walk(genome_download_dir):
        for filename in files:
            if filename.endswith(file_extension) and not filename.endswith(tool + file_extension):
                zip_filename = os.path.join(root, filename)
                # Uncompress ".gz" file
                subprocess.call("gunzip " + zip_filename, shell =True)
                unzip_filename = zip_filename[:-3]

                if root.endswith("extra"):
                    id_loc = [i for i,x in enumerate(extra_files) if x == filename][0]
                    assert 'id_loc' in locals() or 'id_loc' in globals(), "Error: problem with name of the extra files provided"
                    taxid = extra_taxid[id_loc]
                    assert 'taxid' in locals() or 'taxid' in globals(), "Error: problem with the taxid of the extra files provided"
                else:
                    # Retrieve assembly accession number for file path
                    assembly_accession = re.findall(r'/viral/([^(]*)/', unzip_filename)
                    assert 'assembly_accession' in locals() or 'assembly_accession' in globals(), "Can't locate assemble accession"
                    # retrieve taxid for file
                    taxid_list = list(assembly_summary.loc[assembly_summary['assembly_accession'] == assembly_accession[0]]["taxid"])
                    assert (len(taxid_list) == 1), "Taxid has " + len(taxid) + "vales. Should only have 1 value"
                    taxid = taxid_list[0]

                # Create new genomic file with rename sequence identifier to comply with tool requirements for custom database
                renamed_file = unzip_filename[:-4] + "." + tool + unzip_filename[-4:]
                with open(renamed_file, 'w') as out_file, open(unzip_filename, 'r') as in_file:
                    for line in in_file:
                        if line[0] == ">":
                            if tool == "kraken":
                                out_file.write(line[:line.index(" ")] + "|kraken:taxid|" + str(taxid[0]) + line[line.index(" "):])
                            else:
                                out_file.write(">" + str(taxid[0]) + "\n")
                        else:
                            out_file.write(line)
                # Delete original file
                subprocess.call("rm " + unzip_filename,shell=True)
                # Compress modified file
                subprocess.call("gzip " + renamed_file, shell =True)

# Count number of dirs and files ### DELETE ####
# count_dir = 0
# count_file = 0
# for root, dirs, filenames in os.walk(genome_download_dir + "refseq/viral/"):
#     for directory in dirs:
#         count_dir+= 1
#     for files in filenames:
#         if files[-3:] == ".gz":
#             if files[-7:] == ".faa.gz":
#                 count_file += 1
 

# Build the database with the renamed genome files from ncbi
def krakenDB_build(genome_download_dir, kraken_db_dir, threads, kraken_kmer, kraken_minimizer, jellyfish_hash_size = False, kraken_max_dbSize = False):
    # Make a kraken database directory
    if not os.path.exists(kraken_db_dir):
        os.makedirs(kraken_db_dir)
    
    # Download taxonomy for Kraken database
    subprocess.call("kraken-build --download-taxonomy --threads " + str(threads) + " --db " + kraken_db_dir, shell = True)

    # Add files downloaded and ready for kraken ("<file>.tax.fna") through krakenDB_download() to kraken library
    for root, subdirs, files in os.walk(genome_download_dir):
        for filename in files:
            if filename.endswith("kraken.fna.gz"):
                zip_filename = os.path.join(root,filename)
                subprocess.call("gunzip " + zip_filename,shell = True)
                unzip_filename = zip_filename[:-3]
                subprocess.call("kraken-build --add-to-library "  + unzip_filename + " --db " + kraken_db_dir, shell = True)
                subprocess.call("gzip " + unzip_filename, shell = True)

    # Build the command to run kraken-build based on parameters specifed
    kraken_command = "kraken-build --build --threads " + str(threads) + " --db " + kraken_db_dir + " --kmer-len " + str(kraken_kmer) + " --minimizer-len " + str(kraken_minimizer)
    if kraken_max_dbSize:
        kraken_command += " --max-db-size " + str(kraken_max_dbSize)

    if jellyfish_hash_size:
        kraken_command += " --jellyfish-hash-size " + jellyfish_hash_size

    subprocess.call(kraken_command, shell = True)
    
    # Clear unnecessary files from kraken database directory
    subprocess.call("kraken-build --clean --db " + kraken_db_dir, shell = True)
#    subprocess.call("tar -czvf kraken_db.tar.gz " + kraken_db_dir, shell = True)

def kaijuDB_build(genome_download_dir, kaiju_db_dir):
    # Make a kaiju database directory
    if not os.path.exists(kaiju_db_dir):
        os.makedirs(kaiju_db_dir)

    # Add all .faa files to one fasta file
    kaijuDB_fasta = kaiju_db_dir + "kaiju_library.faa"
    count = 0
    with open(kaijuDB_fasta, "w") as out_file:
        for root, subdirs, files in os.walk(genome_download_dir + "refseq/"):
            for filename in files:
                if filename.endswith("kaiju.faa.gz"):
                    zip_filename = os.path.join(root,filename)
                    subprocess.call("gunzip " + zip_filename,shell = True)
                    unzip_filename = zip_filename[:-3]
                    with open(unzip_filename, 'r') as in_file:
                        for line in in_file:
                            if line[0] == ">":
                                out_file.write(line[:1] + str(count) + "_" + line[1:])
                                count += 1
                            else:
                                 out_file.write(line)

                    subprocess.call("gzip " + unzip_filename, shell = True)

    os.chdir(kaiju_db_dir)
    # Fetch "nodes.dmp" and "names.dmp"
    if not os.path.exists(kaiju_db_dir + "names.dmp"):
        urllib.urlretrieve('ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz', 'taxdump.tar.gz')
        subprocess.call("tar -xzvf taxdump.tar.gz", shell = True)
        subprocess.call("rm taxdump.tar.gz citations.dmp delnodes.dmp division.dmp gc.prt gencode.dmp merged.dmp readme.txt", shell = True)
        
    # Build Kaiju database
    subprocess.call("mkbwt -n 5 -a ACDEFGHIKLMNPQRSTVWY -o kaiju_library kaiju_library.faa", shell = True)
    subprocess.call("mkfmi kaiju_library", shell = True)
    subprocess.call("rm kaiju_library.faa kaiju_library.bwt kaiju_library.sa", shell = True)

