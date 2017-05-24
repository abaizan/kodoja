
#Analysing data frm kraken output
#python 2.7
import os
import pandas as pd
import numpy as np
#import re
#import matplotlib.pyplot as plt
#import matplotlib
#matplotlib.style.use('ggplot')
from Bio import SeqIO

#directory = os.environ["OUTDIR"]
#os.chdir(directory)
os.chdir('/home/ae42909/Scratch/test_pipeline')


# Name of the data to be analysed
result_data='kraken_nt_results'
label_data='kraken_nt_labels'
result_names=["Classified", "Seq_ID","Tax_ID", "length", "k-mer classification" ]
label_names=["Seq_ID", "Seq_tax"] #(superkingdom, kingdom, phylum, class, order, family, genus, species)

# Read tables in with pandas
kraken_sequences = pd.read_csv(result_data, sep="\t", header = None, names= result_names)
kraken_label = pd.read_csv(label_data, sep="\t", header = None, names= label_names)

# Merge result tables for Kraken - to include taxonomc names to each sequence that was classified - http://chrisalbon.com/python/pandas_join_merge_dataframe.html
kraken_result = pd.merge(kraken_sequences, kraken_label, on='Seq_ID', how='outer')

#  To be able to identify what broad category (plant, virus bacteria, etc) a sequence has been categorised into and NCBI taxonomy table was made with the broad categories next to each taxID and merged with the kraken results table

# Possible categories for NCBI_taxonomy table abreviations and long names:
Div_tax = {'UNA':'Unannotated', 'BCT':'Bacteria', 'ENV':'Environmental samples', 'SYN':'Synthetic', 'PLN':'Plants', 'INV':'Invertebrates', 'VRT':'Other vertebrates', 'MAM':'Other mammals', 'PRI':'Primates', 'ROD':'Rodents', 'VRL':'Viruses', 'PHG':'Phages'}

# Read in NCBI table - made using script /home/scripts_inProcess/NCBI_taxonomy.py
ncbi_tax = pd.read_csv('/home/ae42909/Scratch/kraken/kraken_analysis/customDatabase/NCBI_taxonomy.csv', sep=",")
ncbi_slice = ncbi_tax.iloc[:,[1,2]]

# Merge tables NCBI taxonomy categories with kraken data
kraken_all =  pd.merge(kraken_result, ncbi_slice, on='Tax_ID', how='outer')

# Remove entries from the NCBI table that do not correspond to any in the kraken results table
kraken_all = kraken_all.dropna(subset = ['Seq_ID'])
kraken_all[['Seq_ID']] = kraken_all[['Seq_ID']].astype(int).astype(str) 


# Make a list of "Seq_ID" column value if sequence is unclassified in "Classified" column or classified as VRL (virus) in column "Div_ID". This list will be used to determine which sequences will be further analysed by Kaiju and retranslation.py

unclassified_IDs = kraken_all.loc[(kraken_all.Classified == 'U'), ['Seq_ID']]
VRL_IDs = kraken_all.loc[(kraken_all.Div_ID == 'VRL'), ['Seq_ID']]
reanalyse_IDs = unclassified_IDs['Seq_ID'].tolist() + VRL_IDs['Seq_ID'].tolist()

# Import replaced sequence IDs
# id1_values = open("ID1", "r").read().split('\n')
# id2_values = open("ID2", "r").read().split('\n')
# id_keys = list(range(1,len(id1_values)))
# id1_dic = dict(zip(id_keys, id1_values))
# id2_dic = dict(zip(id_keys, id2_values))

# reanalyse_ID1 = map(id1_dic.get, reanalyse_IDs)
# reanalyse_ID2 = map(id2_dic.get, reanalyse_IDs)

reanalyse_ID1 = [s + "/1" for s in reanalyse_IDs]
reanalyse_ID2 = [s + "/2" for s in reanalyse_IDs]


# Use biopython to make new fastq files of sequences to be reanalysed. 

def reanalyse_subset (input_file, output_file, id_list):
    outfile = open("reanalyse_ID.txt", "w")
    print >> outfile, "\n".join(str(i) for i in id_list)
    outfile.close()
    id_file = "reanalyse_ID.txt"

    wanted = set(line.rstrip("\n").split(None,1)[0] for line in open(id_file))
    print "Found %i unique identifiers in %s" % (len(wanted), id_file)
    
    records = (r for r in SeqIO.parse(input_file, "fastq") if r.id in wanted)
    count = SeqIO.write(records, output_file, "fastq")
    print "Saved %i records from %s to %s" % (count, input_file, output_file)
    if count < len(wanted):
        print "Warning %i IDs not found in %s" % (len(wanted)-count, input_file)

reanalyse_subset(input_file = "PE1_renamed",output_file = "PE1_subset", id_list=reanalyse_ID1)
reanalyse_subset(input_file = "PE2_renamed",output_file = "PE2_subset", id_list=reanalyse_ID2)
      


