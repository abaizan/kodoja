#!/usr/bin/env python
"""Script for running Kodoja database construction modules."""
from __future__ import print_function

import os
import pandas as pd
import argparse
import shutil
import sys

from diagnostic_modules import version
from diagnostic_modules import check_path
from database_modules import download_with_retries
from database_modules import ncbi_download
from database_modules import ncbi_rename_customDB
from database_modules import krakenDB_build
from database_modules import kaijuDB_build

help_text = """Kodoja Build creates a database for use with Kodoja Search."""

help_epilog = """
See also https://github.com/abaizan/kodoja/wiki/Kodoja-Manual
"""


def main():
    """Run kodoja build."""
    parser = argparse.ArgumentParser(description=help_text,
                                     epilog=help_epilog,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--version',
                        action='version',
                        version='Kodoja v' + version)
    parser.add_argument('-o', '--output_dir', type=str, required=True,
                        help='Output directory path, required')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='Number of threads, default=1')
    parser.add_argument('-p', '--host_taxid', type=int, default=False,
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
    # args.extra_files = ['/home/ae42909/Scratch/GCF_000147415.1_v_1.0_genomic.fna.gz']
    # args.extra_taxids = [554065]

    tool_list = ['kraken', 'kaiju']

    args.output_dir += check_path(args.output_dir)
    kraken_db_dir = os.path.join(args.output_dir, "krakenDB")
    kaiju_db_dir = os.path.join(args.output_dir, "kaijuDB")

    # Check args.output_dir exits, else make dir
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Name databases with tag
    if args.db_tag:
        kraken_db_dir += '_' + args.db_tag
        kaiju_db_dir += '_' + args.db_tag

    kraken_db_dir += check_path(kraken_db_dir)
    kaiju_db_dir += check_path(kaiju_db_dir)

    # Ensure extra files have taxIDs/are in the right format and create symlinks
    if args.extra_files:
        # Check there are the same number of files as taxids
        assert len(args.extra_files) == len(args.extra_taxids), \
            "Each extra file provided needs to have a corresponding ncbi taxid"
        # Check extra files provided are compressed and have the right file extension
        for f in args.extra_files:
            if not f.endswith((".fna.gz", ".faa.gz")):
                sys.exit("File extensions need to be either compressed '.fna.gz' "
                         "for genomic data, or '.faa.gz' for protein data. "
                         "Got %r" % f)
        # Make a copy of each file in extra_files into 'extra' directory
        os.makedirs(os.path.join(args.output_dir, "extra/"))
        for extraFile in args.extra_files:
            shutil.copy(extraFile, os.path.join(args.output_dir, 'extra/'))

    # Download virus assembly summary for refseq
    if not os.path.exists(os.path.join(args.output_dir, "viral_assembly_summary.txt")):
        download_with_retries('https://ftp.ncbi.nih.gov/genomes/refseq/viral/assembly_summary.txt',
                              os.path.join(args.output_dir, 'viral_assembly_summary.txt'))
    path_assembly_summary = os.path.join(args.output_dir, "viral_assembly_summary.txt")
    vir_assembly = pd.read_table(path_assembly_summary, sep='\t', skiprows=1, header=0)
    vir_assembly = vir_assembly.rename(columns={'# assembly_accession': 'assembly_accession'})

    # Set subset_vir_assembly and vir_host
    # subset_vir_assembly - list of virus accession names which will be added to
    #                       databases (used in krakenDB_build and kaijuDB_build)
    # vir_host - list of viral taxIDs for plant viruses. A subset of the genomes
    #            in refseq, be added to datbase by setting 'subset_vir_assembly'.
    if args.all_viruses:
        subset_vir_assembly = False
        vir_host = False
    else:
        # After downloading, will filter for plant-host virus
        if not os.path.exists(os.path.join(args.output_dir, "virushostdb.tsv")):
            # os.chdir(args.output_dir)
            download_with_retries('ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv',
                                  os.path.join(args.output_dir, 'virushostdb.tsv'))
        virHost_table = pd.read_csv(os.path.join(args.output_dir, "virushostdb.tsv"), sep="\t").fillna('')
        plnVir = virHost_table[virHost_table['host lineage'].str.contains("Viridiplantae")]
        vir_host = list(plnVir['virus tax id'])

        subset_vir_assembly = list(vir_assembly.assembly_accession[vir_assembly['taxid'].isin(vir_host)])

    for tool in tool_list:
        if args.no_download:
            # Download NCBI genomes
            ncbi_download(tool, args.output_dir, args.download_parallel, args.host_taxid)
            print("DONE with downloading")
        # Rename downloaded genomic files for Kraken or protein file for Kaiju
        ncbi_rename_customDB(tool, args.output_dir, args.host_taxid, args.extra_files, args.extra_taxids)
        print("DONE with renaming")
        # Make Kraken database
        if tool == "kraken":
            krakenDB_build(args.output_dir, kraken_db_dir, args.threads, args.kraken_kmer,
                           args.kraken_minimizer, subset_vir_assembly, args.kraken_tax)
            print("DONE with kraken db")
            if subset_vir_assembly:
                vir_genomes_text = 'plant viruses'
            else:
                vir_genomes_text = 'all viruses'

            other_genomes_text = ''
            if args.host_taxid:
                other_genomes_text += str(args.host_taxid) + ', '
            if args.extra_taxids:
                other_genomes_text += str(args.extra_taxids)
            with open(os.path.join(kraken_db_dir, "log_file.txt"), "w") as out_file:
                text = 'output_dir = ' + args.output_dir + '\n'
                text += 'kraken_kmer = ' + str(args.kraken_kmer) + '\n'
                text += 'kraken_minimizer = ' + str(args.kraken_minimizer) + '\n'
                text += 'Viral genomes added to db = ' + vir_genomes_text + '\n'
                text += 'Other genome added to db = ' + other_genomes_text + '\n'
                out_file.write(text)
        elif tool == "kaiju":
            # Make Kaiju database
            kaijuDB_build(args.output_dir, kaiju_db_dir, subset_vir_assembly)
            with open(os.path.join(kaiju_db_dir, "log_file.txt"), "w") as out_file:
                text = 'output_dir = ' + args.output_dir + '\n'
                text += 'Viral genomes added to db = ' + vir_genomes_text + '\n'
                text += 'Other genome added to db = ' + other_genomes_text + '\n'
                out_file.write(text)


try:
    main()
except KeyboardInterrupt:
    msg = "Kodoja was interupted by the user.\n"
    sys.exit(msg)
except Exception:
    import traceback
    msg = ("Kodoja failed unexpectedly with the following:\n"
           "\n"
           "%s\n"
           "\n"
           "If this happens consistently, please check you are using the\n"
           "latest version, then check the issue tracker to see if this is\n"
           "a known issue, and if not please report the problem:\n"
           "\n"
           "https://github.com/abaizan/kodoja/issues\n"
           "\n"
           "Kodoja aborted.\n" % traceback.format_exc())
    sys.exit(msg)
