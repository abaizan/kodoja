#python 2.7
import os
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use('ggplot')
import time
#from Bio import SeqIO

#directory = os.environ["OUTDIR"]
directory = '/home/ae42909/Scratch/synthPotato_pipeline/'

# Import initial results from sequence_reanalysis.py
kraiju = pd.read_csv(directory + 'kraiju_table', header = 0, sep='\t')

# Summary table - list tax_ids identified in classification, count how many seuqneces fall under each by kraken, kaiju and when the two are in agreement, 
def summary_table(kraiju_data):
    ncbi_tax = pd.read_csv('/home/ae42909/Scratch/kraken/kraken_analysis/customDatabase/NCBI_taxonomy.csv', sep=",")
    kraken_class = dict(kraiju_data['kraken_Tax_ID'].value_counts())
    kaiju_class = dict(kraiju_data['kaiju_Tax_ID'].value_counts())
    combined_class = dict(kraiju_data['combined_result'].value_counts())

    class_list = list(set(list(kraken_class.keys()) + list(kaiju_class.keys())))
    class_list.remove(0)

    table_summary = pd.DataFrame(columns=['Name', 'Tax_ID', 'kraken', 'kaiju','combined'])
    table_summary['Tax_ID'] =  map(int, class_list)
    table_summary['kraken'] = table_summary['Tax_ID'].map(kraken_class)
    table_summary['kaiju'] = table_summary['Tax_ID'].map(kaiju_class)
    table_summary['combined'] = table_summary['Tax_ID'].map(combined_class)


    def name_lookup(each_row):
        try:
            tax_name = ncbi_tax[ncbi_tax['Tax_ID'] == each_row['Tax_ID']]['Name'].item()
        except:
            tax_name = 'Not found'
        return tax_name  

    table_summary['Name'] = table_summary.apply(lambda row: name_lookup(row), axis=1)
    table_summary = table_summary.sort_values(['combined'], ascending=False)

    return table_summary

kraiju_summary = summary_table(kraiju)

ncbi_nodes = pd.read_csv('/home/ae42909/Scratch/kaiju/kaijudb/nodes.dmp', sep="\t", header = None)
# nodes.dmp file consists of taxonomy nodes. The description for each node includes the following
# fields:
# 	tax_id					-- node id in GenBank taxonomy database
#  	parent tax_id				-- parent node id in GenBank taxonomy database
#  	rank					-- rank of this node (superkingdom, kingdom, ...) 
#  	embl code				-- locus-name prefix; not unique
#  	division id				-- see division.dmp file
#  	inherited div flag  (1 or 0)		-- 1 if node inherits division from parent
#  	genetic code id				-- see gencode.dmp file
#  	inherited GC  flag  (1 or 0)		-- 1 if node inherits genetic code from parent
#  	mitochondrial genetic code id		-- see gencode.dmp file
#  	inherited MGC flag  (1 or 0)		-- 1 if node inherits mitochondrial gencode from parent
#  	GenBank hidden flag (1 or 0)            -- 1 if name is suppressed in GenBank entry lineage
#  	hidden subtree root flag (1 or 0)       -- 1 if this subtree has no sequence data yet
#  	comments				-- free-text comments and citations


for row in range(0,len(kraiju_summary)):
    row_node = ncbi_nodes[ncbi_nodes[0] == kraiju_summary['Tax_ID'][row]]
    while row_rank[4].item() !=

    if row_rank[4].item() == 'species':
        parent_node = row_node[2].item()
        row_node = ncbi_nodes[ncbi_nodes[0] == kraiju_summary['Tax_ID'][row]]
 
row_node = ncbi_nodes[ncbi_nodes[0] == kraiju_summary['Tax_ID'][1012]]

# Test if 'genus' and 'phylum' are always at the same location after species
species = ncbi_nodes[ncbi_nodes[4] == 'species']
after_species = []

for row in range(0, len(species)):
    row_node = species.iloc[row]
    parent_tax = row_node[2].item()
    parent_row = ncbi_nodes[ncbi_nodes[0] == parent_tax]
    after_species.append(parent_row[4].item())
