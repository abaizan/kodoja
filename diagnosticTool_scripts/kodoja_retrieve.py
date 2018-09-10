#!/usr/bin/env python
"""Retrieve viral sequences of interest."""
from __future__ import print_function

import argparse
import os
import pickle
import sys

import pandas as pd
from diagnostic_modules import version
from diagnostic_modules import sequence_subset
from diagnostic_modules import check_path


help_text = """Kodoja Retrieve is used with the output of Kodoja Search to
give subsets of your input sequencing reads matching viruses."""

help_epilog = """
The main output of ``kodoja_search.py`` is a file called ``virus_table.txt``
(a table summarising the potential viruses found), but the specified output
directory will also contain ``kodoja_VRL.txt`` (a table listing the read
identifiers). This second file is used as input to ``kodoja_retrieve.py``
along with the original sequencing read files.

A sub-directory ``subreads/`` will be created in the output directory,
which will include either FASTA or FASTQ files named as follows:

* ``subset_files/virus_all_sequences1.fasta`` FASTA output
* ``subset_files/virus_all_sequences1.fastq`` FASTQ output

And, for paired end datasets,

* ``subset_files/virus_all_sequences2.fasta`` FASTA output
* ``subset_files/virus_all_sequences2.fastq`` FASTQ output

However, if the ``-t 12345`` option is used rather than ``virus_all_...``
the files will be named ``virus_12345_...`` instead.
"""


def main():
    """Run kodoka retrieve."""
    parser = argparse.ArgumentParser(description=help_text,
                                     epilog=help_epilog,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--version',
                        action='version',
                        version='Kodoja v' + version)
    parser.add_argument('-o', '--file_dir', type=str, required=True,
                        help='Path to directory of kodoja_search results, required')
    parser.add_argument('-r1', '--read1', type=str, required=True,
                        help='Read 1 file path, required')
    parser.add_argument('-r2', '--read2', type=str, default=False,
                        help='Read 2 file path, default: False')
    parser.add_argument('-f', '--user_format', type=str, default='fastq',
                        help='Sequence data format, default: fastq')
    parser.add_argument('-t', '--taxID', type=int, default=False,
                        help='Virus tax ID for subsetting, default: All viral sequences')
    parser.add_argument('-g', '--genus', action='store_true',
                        help='Include sequences classified at genus')
    parser.add_argument('-s', '--stringent', action='store_true',
                        help='Only subset sequences identified by both tools')
    args = parser.parse_args()

    table_summary = pd.read_csv(os.path.join(args.file_dir, "virus_table.txt"), sep="\t", header=0,
                                index_col=False)
    kodoja_vrl = pd.read_csv(os.path.join(args.file_dir, "kodoja_VRL.txt"), sep="\t", header=0,
                             index_col=False).fillna('')
    args.file_dir += check_path(args.file_dir)
    output_dir = os.path.join(args.file_dir, "subset_files/")

    # Create directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if args.taxID:
        TaxId_out = [args.taxID]
        label = 'virus_' + str(args.taxID)
        if args.genus:
            more_taxids = []
            with open(os.path.join(args.file_dir, 'genus_taxid.pkl'), 'rb') as id_dict:
                genus_taxid = pickle.load(id_dict)
            for sp_taxid in TaxId_out:
                genus_name = table_summary.Genus[table_summary['Species TaxID'] == sp_taxid].values[0]
                if genus_name in genus_taxid:
                    for items in genus_taxid[genus_name]:
                        if items not in more_taxids:
                            more_taxids.append(items)
            for new_taxid in more_taxids:
                TaxId_out.append(new_taxid)
    else:
        kraken_taxid = list(kodoja_vrl.kraken_tax_ID[kodoja_vrl['kraken_seq_tax'].str.contains("Viruses")])
        kraken_taxid += list(kodoja_vrl.kraken_tax_ID[kodoja_vrl['kraken_seq_tax'].str.contains("Viroids")])
        kaiju_taxid = list(kodoja_vrl.kaiju_tax_ID[kodoja_vrl['kaiju_seq_tax'].str.contains("Viruses")])
        TaxId_out = sorted(set(kraken_taxid + kaiju_taxid))
        # TaxId_out = list(table_summary['Species TaxID'])
        label = 'virus_all'

    if args.stringent:
        rows_wanted = kodoja_vrl['combined_result'].isin(TaxId_out)
    else:
        rows_wanted = (kodoja_vrl['kraken_tax_ID'].isin(TaxId_out) |
                       kodoja_vrl['kaiju_tax_ID'].isin(TaxId_out))
    # Since kodoja v0.0.8 the Seq_ID column has been just the ID,
    # but on earlier versions would be full description line -
    # thus splitting on the first white space:
    seqID_wanted = set(_.rstrip("\n").split(None, 1)[0] for _ in kodoja_vrl.loc[rows_wanted, 'Seq_ID'])

    sequence_subset(output_dir, args.read1, label + "_sequences1.", args.user_format,
                    seqID_wanted)

    if args.read2:
        # Cope with case where inputs are example/1 and example/2 style:
        def rename2to1(identifier):
            if identifier.endswith("/1"):
                return identifier[:-1] + "2"
            else:
                return identifier
        assert rename2to1("example/1") == "example/2"
        seqID_wanted = set(rename2to1(_) for _ in seqID_wanted)
        sequence_subset(output_dir, args.read2, label + "_sequences2.", args.user_format,
                        seqID_wanted)


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
