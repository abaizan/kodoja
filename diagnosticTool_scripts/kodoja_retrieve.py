#!/usr/bin/env python
"""Retrieve viral sequences of interest."""
from __future__ import print_function

import pandas as pd
from diagnostic_modules import version
from diagnostic_modules import sequence_subset
import argparse
import pickle

parser = argparse.ArgumentParser(description='Retrieve viral sequences of interest')
parser.add_argument('--version',
                    action='version',
                    version='Kodoja v' + version)
parser.add_argument('-r1', '--read1', type=str, required=True,
                    metavar='', help='Read 1 file path, required')
parser.add_argument('-o', '--out_dir', type=str, required=True,
                    metavar='', help='Output directory path, required')
parser.add_argument('-r2', '--read2', type=str, default=False,
                    metavar='', help='Read 2 file path, default: False')
parser.add_argument('-f', '--user_format', type=str, default='fastq',
                    metavar='', help='Sequence data format, default: fastq')
parser.add_argument('-t', '--taxID', type=int, default=False,
                    metavar='', help='Virus tax ID for subsetting, default: All viral sequences')
args = parser.parse_args()

table_summary = pd.read_csv(args.out_dir + "virus_table.txt", sep="\t", header=0,
                            index_col=False)
kodoja_vrl = pd.read_csv(args.out_dir + "kodoja_VRL.txt", sep="\t", header=0,
                         index_col=False)

if args.taxID:
    TaxId_out = [args.taxID]
    label = 'virus_' + str(args.taxID)

else:
    TaxId_out = list(table_summary.Tax_ID)
    label = 'virus_all'

rows_wanted = (kodoja_vrl['kraken_tax_ID'].isin(TaxId_out) |
               kodoja_vrl['kaiju_tax_ID'].isin(TaxId_out))
seqID_wanted = list(kodoja_vrl.loc[rows_wanted, 'Seq_ID'])

sequence_subset(args.out_dir, args.read1, label + "_sequences1.", args.user_format,
                seqID_wanted, label + '_sequences1.txt')
if args.read2:
    print('We got to importing dict')
    with open(args.out_dir + 'ids1.pkl', 'rb') as id_dict:
        ids1 = pickle.load(id_dict)
    with open(args.out_dir + 'ids2.pkl', 'rb') as id_dict:
        ids2 = pickle.load(id_dict)
    iv_ids1 = dict((v, k) for k, v in ids1.iteritems())
    kodoja_vrl["Seq_ID"] = kodoja_vrl["Seq_ID"].map(iv_ids1)
    kodoja_vrl["Seq_ID"] = kodoja_vrl["Seq_ID"].map(ids2)
    seqID_wanted = list(kodoja_vrl.loc[rows_wanted, 'Seq_ID'])
    sequence_subset(args.out_dir, args.read2, label + "_sequences2.", args.user_format,
                    seqID_wanted, label + '_sequences2.txt')
