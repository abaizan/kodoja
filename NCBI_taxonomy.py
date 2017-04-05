# I would like to know what each species etc has a taxomony ID, particularly as Kaiju doesn't have an output version with spcies names. ALso I need to identify which IDs are virus so I can extract them to put into Kaiju (or other for protein analysis)

import os
import pandas as pd
import numpy as np

os.chdir('/home/ae42909/Scratch/kraken/taxonomy')

tax_files={'nodes':'nodes.dmp','meged':'merged.dmp', 'division':'division.dmp','names':'names.dmp',}

for file in tax_files.keys():
    TaxFile = pd.read_csv(tax_files[file], sep="\t", header = None)
    exec(file + " = TaxFile")

# All file columns are separated by | (which makes a new column, so lways twice as many columns)
# nodes.dmp file consists of taxonomy nodes. The description for each node includes the following
#fields:
#  0      tax_id                                  -- node id in GenBank taxonomy database
#  2      parent tax_id                           -- parent node id in GenBank taxonomy database
#  4      rank                                    -- rank of this node (superkingdom, kingdom, ...)
#  6      embl code                               -- locus-name prefix; not unique
#  8      division id                             -- see division.dmp file
# 10      inherited div flag  (1 or 0)            -- 1 if node inherits division from parent
#        genetic code id                         -- see gencode.dmp file
#        inherited GC  flag  (1 or 0)            -- 1 if node inherits genetic code from parent
#        mitochondrial genetic code id           -- see gencode.dmp file
#        inherited MGC flag  (1 or 0)            -- 1 if node inherits mitochondrial gencode from parent
#        GenBank hidden flag (1 or 0)            -- 1 if name is suppressed in GenBank entry lineage
#        hidden subtree root flag (1 or 0)	-- 1 if this subtree has no sequence data yet
#        comments                                -- free-text comments and citations


# Merged nodes file (merged.dmp):
#        old_tax_id                              -- id of nodes which has been merged
#        new_tax_id                              -- id of nodes which is result of merging
#
# Divisions file (division.dmp):
#        division id                             -- taxonomy database division id
#        division cde                            -- GenBank division code (three characters)
#        division name                           -- e.g. BCT, PLN, VRT, MAM, PRI...
#        comments


# Taxonomy names file (names.dmp):
#        tax_id                                  -- the id of node associated with this name
#        name_txt                                -- name itself
#        unique name                             -- the unique variant of this name if name not unique
#        name class                              -- (synonym, common name, ...)

# You could take parent ID to try and work out a way of classifing everything like Kraken does (ranks)
names_need = [0,2,6]
tax_names = names[names.columns[names_need]]
tax_names.columns = ['Tax_ID', 'Name', 'Class']

nodes_need = [0,8]
tax_nodes = nodes[nodes.columns[nodes_need]]
tax_nodes.columns = ['Tax_ID', 'Div_ID']

tax_merge = pd.merge(tax_nodes, tax_names, on='Tax_ID', how='outer')
tax_reduce = tax_merge.loc[tax_merge['Class'] == 'scientific name'] #otherwise you get the same ID with differnt names
Div_dict = division.set_index([0])[2].to_dict()
tax_complete = tax_reduce.replace({'Div_ID': Div_dict}) #Replace the codes with BCT PLN etc

# Export to kraken directory
path = '/home/ae42909/Scratch/kraken/kraken_analysis/customDatabase'
tax_complete.to_csv(os.path.join(path,r'NCBI_taxonomy.csv'))
