"""Database construction modules."""
import subprocess
import pandas as pd
import os
import re
import sys
import time


try:
    # Python 3
    from urllib.request import urlretrieve, urlcleanup
    from urllib.error import URLError
except ImportError:
    # Python 2
    from urllib import urlretrieve, urlcleanup
    URLError = IOError


def download_with_retries(url, destination, retries=5):
    """Download file using urlretrieve, with automatic retries.

    If the n-th attempt fails, will wait for n-seconds before
    trying again.

    If the final retry fails, will abort the script with message
    to stderr.
    """
    for attempt in range(1, retries + 1):
        urlcleanup()  # Seems important with FTP on Python 2.7
        try:
            return urlretrieve(url, destination)
        except URLError as err:
            if attempt < retries:
                time.sleep(attempt)
                sys.stderr.write("Will retry downloading %r attempt %i\n"
                                 % (url, attempt + 1))
            else:
                sys.stderr.write("ERROR: Failed to download %r\n" % url)
                raise err


# Download refseq files from ncbi ftp site - use ncbi-genome-download
def ncbi_download(tool, genome_download_dir, parallel, host_taxid):
    """Download genomic or protein data from NCBI ftp site using ncbi-genome-download."""
    assert (tool == "kraken") | (tool == "kaiju"),\
        "Argument 'tool' must be either 'kraken' or 'kaiju'."
    if tool == "kraken":
        file_format = "fasta"
    else:
        file_format = "protein-fasta"

    # Check directory exists
    if not os.path.exists(genome_download_dir):
        os.makedirs(genome_download_dir)

    ngd_command = "ncbi-genome-download -F " + file_format + " -o " + genome_download_dir

    if host_taxid:
        # Single host ID, so no need to set the parallel option
        taxid_ngd_command = ngd_command + " --species-taxid " + str(host_taxid) + " plant"
        subprocess.check_call(taxid_ngd_command, shell=True)

    ngd_command += " --parallel " + str(parallel) + " viral"
    subprocess.check_call(ngd_command, shell=True)


def ncbi_rename_customDB(tool, genome_download_dir, host_taxid, extra_files=False, extra_taxid=False):
    """Rename ncbi data files for custom databases.

    To add NCBI files to database Kraken and Kaiju require the files to have formatted
    identifiers. This script modifies identifiers of files ending in .fna to kraken
    format, and files ending in .faa to kaiju format. Once renamed, original files are
    deleted.

    """
    assert (tool == "kraken") | (tool == "kaiju"), "Argument 'tool' must be either 'kraken' or 'kaiju'."
    if tool == "kraken":
        file_extension = ".fna.gz"
    else:
        file_extension = ".faa.gz"

    path_assembly_summary = os.path.join(genome_download_dir, "viral_assembly_summary.txt")
    assembly_summary = pd.read_table(path_assembly_summary, sep='\t',
                                     skiprows=1, header=0)
    assembly_summary.rename(columns={'# assembly_accession': 'assembly_accession'},
                            inplace=True)  # rename column to exclude "#"

    # If extra files, create dictionary containing file names and taxIDs
    if extra_files:
        new_extra = []
        for extraFiles in extra_files:
            new_extra.append(extraFiles.split('/')[-1])
        extra_dict = dict(zip(new_extra, extra_taxid))

    kaiju_count = 1  # Count for protein sequences
    for root, subdirs, files in os.walk(genome_download_dir):
        for filename in files:
            if filename.endswith(file_extension) and not filename.endswith(tool + file_extension):
                zip_filename = os.path.join(root, filename)
                subprocess.check_call("gunzip " + zip_filename, shell=True)  # Uncompress ".gz" file
                unzip_filename = zip_filename[:-3]

                if root.endswith("extra"):
                    taxid = extra_dict[filename]
                    assert 'taxid' in locals() or 'taxid' in globals(),\
                        "Error: no taxid assigned for extra files provided"
                elif root.split('/')[-2] == 'plant':
                    taxid = host_taxid
                else:
                    # Retrieve assembly accession number for file path
                    assembly_accession = re.findall(r'/viral/([^(]*)/', unzip_filename)
                    assert 'assembly_accession' in locals() or 'assembly_accession' in globals(),\
                        "Can't locate assemble accession"
                    # retrieve taxid for file
                    taxid_list = list(assembly_summary.loc[assembly_summary['assembly_accession'] == assembly_accession[0]]["taxid"])
                    assert (len(taxid_list) == 1),\
                        "Taxid has " + len(taxid) + "values. Should only have 1 value"
                    taxid = taxid_list[0]

                # Create new genomic file with rename sequence identifier to comply with tool
                #  requirements for custom database
                renamed_file = unzip_filename[:-4] + "." + tool + unzip_filename[-4:]
                with open(renamed_file, 'w') as out_file, open(unzip_filename, 'r') as in_file:
                    for line in in_file:
                        if line[0] == ">":
                            if tool == "kraken":
                                if " " in line:
                                    insert = line.index(" ")
                                else:
                                    insert = len(line) - 1
                                out_file.write(line[:insert] + "|kraken:taxid|" + str(taxid) + line[insert:])
                            else:
                                out_file.write(">" + str(kaiju_count) + "_" + str(taxid) + "\n")
                                kaiju_count += 1
                        else:
                            out_file.write(line)
                # Delete original file
                os.remove(unzip_filename)
                # Compress modified file
                subprocess.check_call("gzip " + renamed_file, shell=True)


def krakenDB_build(genome_download_dir, kraken_db_dir, threads, kraken_kmer, kraken_minimizer,
                   subset_vir_assembly, taxonomy, jellyfish_hash_size=False, kraken_max_dbSize=False):
    """Build kraken database with the renamed .fna files from ncbi."""
    # Make a kraken database directory
    if not os.path.exists(kraken_db_dir):
        os.makedirs(kraken_db_dir)

    # Download or create symlink of taxonomy for Kraken database
    if taxonomy:
        os.symlink(taxonomy, os.path.join(kraken_db_dir, "taxonomy"))
    else:
        subprocess.check_call("kraken-build --download-taxonomy --threads " +
                              str(threads) + " --db " + kraken_db_dir, shell=True)

    file_list = []
    # Add files downloaded and ready for kraken ("<file>.kraken.fna") to kraken library
    for root, subdirs, files in os.walk(genome_download_dir):
        for filename in files:
            if subset_vir_assembly:
                if root.split('/')[-1] in subset_vir_assembly and filename.endswith("kraken.fna.gz"):
                    file_list.append(os.path.join(root, filename))
                elif root.split('/')[-2] == 'plant' and filename.endswith("kraken.fna.gz"):
                    file_list.append(os.path.join(root, filename))
                elif root.endswith('extra') and filename.endswith("kraken.fna.gz"):
                    file_list.append(os.path.join(root, filename))
            else:
                if filename.endswith("kraken.fna.gz"):
                    file_list.append(os.path.join(root, filename))

    for genome_file in file_list:
        zip_filename = genome_file
        subprocess.check_call("gunzip " + zip_filename, shell=True)
        unzip_filename = zip_filename[:-3]
        subprocess.check_call("kraken-build --add-to-library " + unzip_filename +
                              " --db " + kraken_db_dir, shell=True)
        subprocess.check_call("gzip " + unzip_filename, shell=True)

    kraken_command = "kraken-build --build --threads " + str(threads) + " --db " + \
                     kraken_db_dir + " --kmer-len " + str(kraken_kmer) + \
                     " --minimizer-len " + str(kraken_minimizer)
    if kraken_max_dbSize:
        kraken_command += " --max-db-size " + str(kraken_max_dbSize)

    if jellyfish_hash_size:
        kraken_command += " --jellyfish-hash-size " + jellyfish_hash_size

    subprocess.check_call(kraken_command, shell=True)

    # Clear unnecessary files from kraken database directory
    # subprocess.check_call("kraken-build --clean --db " + kraken_db_dir, shell=True)


def kaijuDB_build(genome_download_dir, kaiju_db_dir, subset_vir_assembly):
    """Build kraken database with the renamed .faa files from ncbi."""
    # Make a kaiju database directory
    if not os.path.exists(kaiju_db_dir):
        os.makedirs(kaiju_db_dir)

    # Add files downloaded and ready for kaiju ("<file>.kaiju.faa") to one fasta file
    kaijuDB_fasta = os.path.join(kaiju_db_dir, "kaiju_library.faa")
    count = 0
    file_list = []
    for root, subdirs, files in os.walk(genome_download_dir):
        for filename in files:
            if subset_vir_assembly:
                if root.split('/')[-1] in subset_vir_assembly and filename.endswith("kaiju.faa.gz"):
                    file_list.append(os.path.join(root, filename))
                elif root.split('/')[-2] == 'plant' and filename.endswith("kaiju.faa.gz"):
                    file_list.append(os.path.join(root, filename))
                elif root.endswith('extra') and filename.endswith("kaiju.faa.gz"):
                    file_list.append(os.path.join(root, filename))
            else:
                if filename.endswith("kaiju.faa.gz"):
                    file_list.append(os.path.join(root, filename))

    with open(kaijuDB_fasta, "w") as out_file:
        for protein_file in file_list:
            zip_filename = protein_file
            subprocess.check_call("gunzip " + zip_filename, shell=True)
            unzip_filename = zip_filename[:-3]
            with open(unzip_filename, 'r') as in_file:
                for line in in_file:
                    if line[0] == ">":
                        out_file.write(line[:1] + str(count) + "_" + line[1:])
                        count += 1
                    else:
                        out_file.write(line)

            subprocess.check_call("gzip " + unzip_filename, shell=True)

    try:
        # Assume kaiju v1.7.0 onwards:
        subprocess.check_output(["kaiju-mkbwt", "-help"])
        prefix = "kaiju-"
    except OSError:  # Expect FileNotFoundError on Python 3.3+
        # kaiju prior to v1.7.0
        prefix = ""

    # Build Kaiju database
    subprocess.check_call(prefix + "mkbwt -n 5 -a ACDEFGHIKLMNPQRSTVWY -o " +
                          os.path.join(kaiju_db_dir, "kaiju_library") + " " +
                          os.path.join(kaiju_db_dir, "kaiju_library.faa"),
                          shell=True)
    subprocess.check_call(prefix + "mkfmi " + os.path.join(kaiju_db_dir, "kaiju_library"),
                          shell=True)
    os.remove(os.path.join(kaiju_db_dir, "kaiju_library.faa"))
    os.remove(os.path.join(kaiju_db_dir, "kaiju_library.bwt"))
    os.remove(os.path.join(kaiju_db_dir, "kaiju_library.sa"))
