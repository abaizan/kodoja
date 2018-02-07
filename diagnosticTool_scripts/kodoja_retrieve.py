#!/usr/bin/env python
"""Retrieve viral sequences of interest."""
from __future__ import print_function

import os
import pandas as pd
from diagnostic_modules import version
from diagnostic_modules import sequence_subset
from diagnostic_modules import check_path
import argparse
import pickle

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

table_summary = pd.read_csv(args.file_dir + "virus_table.txt", sep="\t", header=0,
                            index_col=False)
kodoja_vrl = pd.read_csv(args.file_dir + "kodoja_VRL.txt", sep="\t", header=0,
                         index_col=False)
args.file_dir += check_path(args.file_dir)
output_dir = args.file_dir + 'subset_files/'

# Create directory
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

if args.taxID:
    TaxId_out = [args.taxID]
    label = 'virus_' + str(args.taxID)
else:
    TaxId_out = list(table_summary['Species TaxID'])
    label = 'virus_all'

if args.genus:
    more_taxids = []
    with open(args.file_dir + 'genus_taxid.pkl', 'rb') as id_dict:
        genus_taxid = pickle.load(id_dict)
    for sp_taxid in TaxId_out:
        genus_name = table_summary.Genus[table_summary['Species TaxID'] == sp_taxid].values[0]
        if genus_name in list(genus_taxid.keys()):
            for items in genus_taxid[genus_name]:
                if items not in more_taxids:
                    more_taxids.append(items)
    for new_taxid in more_taxids:
        TaxId_out.append(new_taxid)

if args.stringent:
    rows_wanted = kodoja_vrl['combined_result'].isin(TaxId_out)
else:
    rows_wanted = (kodoja_vrl['kraken_tax_ID'].isin(TaxId_out) |
                   kodoja_vrl['kaiju_tax_ID'].isin(TaxId_out))
seqID_wanted = list(kodoja_vrl.loc[rows_wanted, 'Seq_ID'])

sequence_subset(output_dir, args.read1, label + "_sequences1.", args.user_format,
                seqID_wanted, label + '_sequences1.txt')
if args.read2:
    with open(args.file_dir + 'ids1.pkl', 'rb') as id_dict:
        ids1 = pickle.load(id_dict)
    with open(args.file_dir + 'ids2.pkl', 'rb') as id_dict:
        ids2 = pickle.load(id_dict)
    iv_ids1 = dict((v, k) for k, v in ids1.iteritems())
    kodoja_vrl["Seq_ID"] = kodoja_vrl["Seq_ID"].map(iv_ids1)
    kodoja_vrl["Seq_ID"] = kodoja_vrl["Seq_ID"].map(ids2)
    seqID_wanted = list(kodoja_vrl.loc[rows_wanted, 'Seq_ID'])
    sequence_subset(output_dir, args.read2, label + "_sequences2.", args.user_format,
                    seqID_wanted, label + '_sequences2.txt')
