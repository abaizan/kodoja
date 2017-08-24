# python 2.7.13
import subprocess
import pandas as pd
import numpy as np
import os
import urllib
import re
import StringIO
# from Bio import SeqIO


# Download refseq files from ncbi ftp site - use ncbi-genome-download - and the assembly_summary.txt" file for taxid information
def ncbi_download(file_format, genome_download_dir, parallel=False):
    ngd_command = "ncbi-genome-download -F " + file_format +  " -o " + genome_download_dir + " viral "
    if parallel:
        ngd_command += "--parallel " + parallel
    subprocess.call(ngd_command, shell = True)
    
# Rename sequence identifiers for ".fna" genomic files downloaded using "ncbi_download"
def kraken_fna_rename(genome_download_dir):
    # Download assembly summary file from ncbi ftp site and load as pandas dataframe
    genome_download_dir += "refseq/"
    os.chdir(genome_download_dir)  
    urllib.urlretrieve('ftp://ftp.ncbi.nih.gov/genomes/refseq/viral/assembly_summary.txt', 'viral_assembly_summary.txt')
    path_assembly_summary = genome_download_dir + "viral_assembly_summary.txt"
    assembly_summary = pd.read_table(path_assembly_summary, sep='\t', skiprows=1, header = 0)
    assembly_summary.rename(columns={'# assembly_accession':'assembly_accession'}, inplace=True)  # rename column to exclude "#"

    # Rename all .fna files so that sequence identifiers comply with Kraken requirments
    for root, subdirs, files in os.walk(genome_download_dir + "viral/"):
        for filename in files:
            if filename[-7:] == ".fna.gz":
                zip_file = os.path.join(root, filename)
                # Uncompress ".gz" file
                subprocess.call("gunzip " + zip_file, shell =True)
                unzip_file = zip_file[:-3]

                # Retrieve assembly accession number for file
                assembly_accession = re.findall(r'/viral/([^(]*)/', unzip_file)
                # retrieve taxid for file
                taxid = list(assembly_summary.loc[assembly_summary['assembly_accession'] == assembly_accession[0]]["taxid"])[0]

                # Create new genomic file with rename sequence identifier to comply with Kraken requirements
                renamed_file = unzip_file[:-4] + ".kraken" + unzip_file[-4:]
                with open(renamed_file, 'w') as out_file, open(unzip_file, 'r') as in_file:
                    for line in in_file:
                        if line[0] == ">":
                            out_file.write(line[:line.index(" ")] + "|kraken:taxid|" + str(taxid) + line[line.index(" "):])
                        else:
                            out_file.write(line)
                # Delete original file
                subprocess.call("rm " + unzip_file,shell=True)
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
def krakenDB_build(genome_download_dir, kraken_db_dir, threads, kraken_kmer, kraken_minimizer, jellyfish_hash_size, kraken_max_dbSize = False):
    # Make a kraken database directory
    os.makedirs(kraken_db_dir)
    
    # Download taxonomy for Kraken database
#    subprocess.call("kraken-build --download-taxonomy --threads " + threads + " --db " + kraken_db_dir, shell = True)

    # Add files downloaded and ready for kraken ("<file>.tax.fna") through krakenDB_download() to kraken library
for path, subdirs, files in os.walk(genome_download_dir + "refseq/"):
    for filename in files:
        if filename.endswith("kraken.fna.gz"):
            subprocess.call("gunzip " + names,shell = True)
            unzip_names = names[:-3]
            subprocess.call("kraken-build --add-to-library "  + os.path.join(path,unzip_names) + " --db " + kraken_db_dir, shell = True)
            subprocess.call("gzip " + unzip_names, shell = True)

    # Build the command to run kraken-build based on parameters specifed
    kraken_command = "kraken-build --build --threads " + threads + " --db " + kraken_db_dir + " --kmer-len " + kraken_kmer + " --minimizer-len " + kraken_minimizer
    if kraken_max_dbSize:
        kraken_command += " --max-db-size " + kraken_max_dbSize

    #kraken_command += " --jellyfish-hash-size 11400M"
    kraken_command += " --jellyfish-hash-size " + jellyfish_hash_size

    subprocess.call(kraken_command, shell = True)
    
    # Clear unnecessary files from kraken database directory
    subprocess.call("kraken-build --clean --db " + kraken_db_dir, shell = True)
#    subprocess.call("tar -czvf kraken_db.tar.gz " + kraken_db_dir, shell = True)

# def krakenDB_shrink(kraken_db_dir, num_kmer, newDB_name):
#     # mkdir newDB_name
#     subprocess.call("kraken-build --shrink " + num_kmer + " --db " + kraken_db_dir + " --new-db " + newDB_name, shell = True)

# def kaijuDB_faa_rename(genome_download_dir):
    genome_download_dir += "refseq/"
    # Download assembly summary file from ncbi ftp site and load as pandas dataframe
    if not os.path.isfile(genome_download_dir + "viral_assembly_summary.txt"):
        os.chdir(genome_download_dir)
        urllib.urlretrieve('ftp://ftp.ncbi.nih.gov/genomes/refseq/viral/assembly_summary.txt', 'viral_assembly_summary.txt')
    path_assembly_summary = genome_download_dir + "viral_assembly_summary.txt"
    assembly_summary = pd.read_table(path_assembly_summary, sep='\t', skiprows=1, header = 0)
    assembly_summary.rename(columns={'# assembly_accession':'assembly_accession'}, inplace=True)  # rename column to exclude "#"

    # Rename all .fna files so that sequence identifiers comply with Kraken requirments
    for root, dirs, filenames in os.walk(genome_download_dir + "viral/"):
        for f in filenames:
            if f[-7:] == ".faa.gz":
                zip_file = os.path.join(root, f)
                # Uncompress ".gz" file
                subprocess.call("gunzip " + zip_file, shell =True)
                unzip_file = zip_file[:-3]

                # Retrieve assembly accession number for file
                assembly_accession = re.findall(r'/viral/([^(]*)/', unzip_file)
                # retrieve taxid for file
                taxid = list(assembly_summary.loc[assembly_summary['assembly_accession'] == assembly_accession[0]]["taxid"])[0]

                # Create new genomic file with rename sequence identifier to comply with Kraken requirements
                renamed_file = unzip_file[:-4] + ".kraken" + unzip_file[-4:]
                with open(renamed_file, 'w') as out_file, open(unzip_file, 'r') as in_file:
                    for line in in_file:
                        if line[0] == ">":
                            out_file.write(line[:line.index(" ")] + "|kraken:taxid|" + str(taxid) + line[line.index(" "):])
                        else:
                            out_file.write(line)
                # Delete original file
                subprocess.call("rm " + unzip_file,shell=True)
                # Compress modified file
                subprocess.call("gzip " + renamed_file, shell =True)

