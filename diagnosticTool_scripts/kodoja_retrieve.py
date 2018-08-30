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


def main():
    """Run kodoka retrieve."""
    parser = argparse.ArgumentParser(description='Retrieve viral sequences of interest')
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
    seqID_wanted = list(kodoja_vrl.loc[rows_wanted, 'Seq_ID'])

    sequence_subset(output_dir, args.read1, label + "_sequences1.", args.user_format,
                    seqID_wanted, label + '_sequences1.txt')

    if args.read2:
        with open(os.path.join(args.file_dir, 'ids1.pkl'), 'rb') as id_dict:
            ids1 = pickle.load(id_dict)
        with open(os.path.join(args.file_dir, 'ids2.pkl'), 'rb') as id_dict:
            ids2 = pickle.load(id_dict)
        iv_ids1 = dict((v, k) for k, v in ids1.items())
        kodoja_vrl["Seq_ID"] = kodoja_vrl["Seq_ID"].map(iv_ids1)
        kodoja_vrl["Seq_ID"] = kodoja_vrl["Seq_ID"].map(ids2)
        seqID_wanted = list(kodoja_vrl.loc[rows_wanted, 'Seq_ID'])
        sequence_subset(output_dir, args.read2, label + "_sequences2.", args.user_format,
                        seqID_wanted, label + '_sequences2.txt')


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
