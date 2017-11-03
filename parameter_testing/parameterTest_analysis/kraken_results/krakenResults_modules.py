import os
import math
import pandas as pd

os.chdir('/home/ae42909/viral_diagnostics/diagnosticTool_scripts/')
from diagnostic_modules import *

# to test k-mer size (kraken)
def kraken_testSummary_table(kraken_data, tag):
    kraken_class = dict(kraken_data['Tax_ID'].value_counts())
    kraken_levels = pd.Series(kraken_data.Seq_tax.values,index=kraken_data.Tax_ID).to_dict()

    levels_dict = kraken_levels.copy()
    levels_dict.pop(0, None)
    levels_tax = {key: list(map(str, value.split('|'))) for key, value in levels_dict.items()}
    LCA_tax = {}
    for key, tax in levels_tax.items():
        if tax[-1][0] != 's':
            LCA_tax[key] = tax[-1]
            levels_tax.pop(key, tax) # remove entry

    species_dict = {}
    for key in levels_tax:
        species_dict[key] = " ".join(levels_tax[key][-1][3:].split("_"))

    associated_tax = {}
    for key_species in levels_tax:
        associated_tax_list = [key_species]
        for key_lca in LCA_tax:
            if LCA_tax[key_lca] in levels_tax[key_species]:
                associated_tax_list.append(key_lca)
        associated_tax[key_species] = associated_tax_list

    table_summary = pd.DataFrame(columns=['Species', 'Tax_ID', 'kraken'])
    table_summary['Tax_ID'] =  map(int, levels_tax.keys())
    table_summary['kraken'] = table_summary['Tax_ID'].map(kraken_class)
    table_summary['Species'] = table_summary['Tax_ID'].map(species_dict)
    table_summary = table_summary.sort_values('kraken', ascending=False)
    table_summary = table_summary.reset_index(drop=True)
    table_summary.to_csv('virus_table_' + tag, sep = '\t', index = False)

    return table_summary

# To run
def kmerResults(wdir, kraken_colNames, expected_viruses):
    kmer_results = pd.DataFrame(columns=expected_viruses.keys()) 
    for root, subdirs, files in os.walk(wdir):
        for filename in files:
            if 'log_file_' in filename:
                with open(os.path.join(root, filename), 'r') as in_file:
                    for lines in in_file:
                        if 'kraken database =' in lines:
                            kmer = int(lines.partition('krakenDB_k')[2].rpartition('_m')[0])      
                tag = filename.split('log_file_',1)[1]
                kraken_table = 'kraken_table_' + tag
                kraken_label = 'kraken_labels_' + tag
                kraken_fullTable = format_result_table(wdir, kraken_table, kraken_label, kraken_colNames)
                result_table = kraken_testSummary_table(kraken_fullTable, tag)

                list_idx = []
                for keys in expected_viruses:
                    idx = result_table.index[result_table['Tax_ID'] == expected_viruses[keys]]
                    assert len(idx) <= 1, 'Tax ID is in more than one entry in result table'
                    if len(idx) < 1:
                        table_idx = float('nan')
                    elif len(idx) == 1:
                        table_idx = idx[0]
                    list_idx.append(table_idx)
                kmer_results.loc[kmer] = list_idx
                kmer_results.to_csv('results_table.txt', sep = '\t')
    

