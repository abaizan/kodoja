import os
import math
import pandas as pd
import sys

sys.path.insert(0,'/home/ae42909/viral_diagnostics/diagnosticTool_scripts/')
from diagnostic_modules import *

# to test parameter
def tool_testSummary_table(tool, tool_data, tag):
    tool_class = dict(tool_data['Tax_ID'].value_counts())
    tool_levels = pd.Series(tool_data.Seq_tax.values,index=tool_data.Tax_ID).to_dict()

    levels_dict = tool_levels.copy()
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

    table_summary = pd.DataFrame(columns=['Species', 'Tax_ID', tool])
    table_summary['Tax_ID'] =  map(int, levels_tax.keys())
    table_summary[tool] = table_summary['Tax_ID'].map(tool_class)
    table_summary['Species'] = table_summary['Tax_ID'].map(species_dict)
    table_summary = table_summary.sort_values(tool, ascending=False)
    table_summary = table_summary.reset_index(drop=True)
    table_summary.to_csv('virus_table_' + tag, sep = '\t', index = False)

    return table_summary

# To sumamrise kraken k-mer results
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
                result_table = tool_testSummary_table('kraken', kraken_fullTable, tag)

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


# To sumarise results
def tool_Results(tool, wdir, tool_colNames, expected_viruses):
    for root, subdirs, files in os.walk(wdir):
        for filename in files:
            if 'log_file_' in filename:
                tag = filename.split('log_file_',1)[1]
                tool_table = tool + '_table_' + tag
                tool_label = tool + '_labels_' + tag
                tool_fullTable = format_result_table(wdir, tool_table, tool_label, tool_colNames)
                result_table = tool_testSummary_table(tool, tool_fullTable, tag)
                total = sum(result_table[tool])
                result_table['Percent'] = (result_table[tool]/total)*100
                expectedV_table = result_table[result_table['Tax_ID'].isin(expected_viruses.values())]
                expectedV_table.to_csv(wdir + 'expectedVirus_table_' + tag, sep = '\t', index = False)
                sum_expectedV = sum(expectedV_table.Percent)
                with open(os.path.join(root, filename), 'a') as in_file:
                    in_file.write('expected virus score = ' + str(sum_expectedV) + '\n')








