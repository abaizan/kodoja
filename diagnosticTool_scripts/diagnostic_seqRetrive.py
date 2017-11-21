import os
import sys
sys.path.insert(0, '/home/ae42909/viral_diagnostics/diagnosticTool_scripts/')
from diagnostic_modules import *

file1 = "/home/ae42909/Scratch/smallTest_data/100_Potato_withViruses_1.fastq"
file2 = False
user_format = "fastq"
out_dir = "/home/ae42909/Scratch/100Seq_krakenDB_viral/"
userList_TaxID = False

table_summary = pd.read_csv(out_dir + "virus_table.txt", sep="\t", header = 0,
                       index_col=False)
kodoja_vrl = pd.read_csv(out_dir + "kodoja_VRL.txt", sep="\t", header = 0,
                       index_col=False)

if userList_TaxID:
    TaxId_out = userList_TaxID
    
else:
    TaxId_out = list(table_summary.Tax_ID)

rows_wanted = (kodoja_vrl['kraken_tax_ID'].isin(TaxId_out) |
               kodoja_vrl['kaiju_tax_ID'].isin(TaxId_out))
seqID_wanted = list(kodoja_vrl.loc[rows_wanted, 'Seq_ID'])

sequence_subset(out_dir, file1, "viral_sequences1.", user_format,
                seqID_wanted, 'viral_sequences1.txt')
if file2:
    with open(out_dir + 'ids1.pkl', 'rb') as id_dict:
        ids1 = pickle.load(id_dict)
    with open(out_dir + 'ids2.pkl', 'rb') as id_dict:
        ids2 = pickle.load(id_dict)
    iv_ids1 = dict((v,k) for k,v in ids1.iteritems())
    kodoja_vrl = kodoja_vrl.replace({"Seq_ID": iv_ids1})
    kodoja_vrl = kodoja_vrl.replace({"Seq_ID": ids2})
    seqID_wanted = list(kodoja_vrl.loc[rows_wanted, 'Seq_ID'])
    sequence_subset(out_dir, file2, "viral_sequences2.", user_format,
                    seqID_wanted, 'viral_sequences2.txt')
