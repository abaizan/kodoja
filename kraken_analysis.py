#Analysing data frm kraken output
#python 2.7
import os
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use('ggplot')

os.chdir('/mnt/shared/scratch/ae42909/201609_BBSRC_Diagnostics/kraken/kraken_analysis/customDatabase')

# Name of the data to be analysed
result_data='PotatoViruses_CD_results'
label_data='Potato_withViruses_classified.labels'
result_names=["Classified", "Seq_ID","Tax_ID", "length", "k-mer classification" ]
label_names=["Seq_ID", "Seq_tax"] #(superkingdom, kingdom, phylum, class, order, family, genus, species)

# Read tables in with pandas
kraken_result = pd.read_csv(result_data, sep="\t", header = None, names= result_names)
kraken_label = pd.read_csv(label_data, sep="\t", header = None, names= label_names)

# Merge result tables for Kraken - to include taxonomc names to each sequence that was classified - http://chrisalbon.com/python/pandas_join_merge_dataframe.html
kraken_all = pd.merge(kraken_result, kraken_label, on='Seq_ID', how='outer')

### To be able to identify what broad category (plant, virus bacteria, etc) a sequence has been categorised into and NCBI taxonomy table was made with the broad categories next to each taxID and merged with the kraken results table

# Possible categories for NCBI_taxonomy table abreviations and long names:
Div_tax = {'UNA':'Unannotated', 'BCT':'Bacteria', 'ENV':'Environmental samples', 'SYN':'Synthetic', 'PLN':'Plants', 'INV':'Invertebrates', 'VRT':'Other vertebrates', 'MAM':'Other mammals', 'PRI':'Primates', 'ROD':'Rodents', 'VRL':'Viruses', 'PHG':'Phages'}

# Read in NCBI table - made using script /home/scripts_inProcess/NCBI_taxonomy.py
ncbi_tax = pd.read_csv('NCBI_taxonomy.csv', sep=",")
ncbi_slice = ncbi_tax.iloc[:,[1,2]]

# Merge tablesNCBI taxonomy categories with kraken data
kraken_all =  pd.merge(kraken_all, ncbi_slice, on='Tax_ID', how='outer')


#### Testing synthetic data ####

### When concatinating libraries to make synthetic data, each had an ID that was kept for identification (first value dictionary). Second value in dict is the correct taxID from NCBI - to determine if the sequences were classified correctlyfor the correct species. The Third value in dict is the correct genus of the species 
library_ID =  {}
library_ID['Tabacco_etch']=['SRR3466597', 12227, 'Potyvirus']
library_ID['UBSV']=['ERR996011', 946046, 'Ipomovirus']
library_ID['BSV']=['ERR996013', 137758, 'Ipomovirus']
library_ID['Potato']=['SRR1207289', 4113, 'Solanum']

# Summarise kraken data: percent of sequences (from each key in library_ID) which are unclassified,  classified correctly or classified incorrectly at the species level
kraken_summary = pd.DataFrame(columns=['Unclassified', 'Classified_correctly','Classified_incorrectly'], index=library_ID.keys())
 
for key in library_ID.keys():
    unclassified_total = kraken_all[(kraken_all.Seq_ID.str.contains(library_ID[key][0])) & (kraken_all['Classified'] == 'U')].count()

    classified_correct = kraken_all[(kraken_all.Seq_ID.str.contains(library_ID[key][0])) & (kraken_all['Classified'] == 'C') & (kraken_all['Tax_ID'] == library_ID[key][1])].count()

    classified_incorrect = kraken_all[(kraken_all.Seq_ID.str.contains(library_ID[key][0])) & (kraken_all['Classified'] == 'C') & ~(kraken_all['Tax_ID'] == library_ID[key][1])].count()

    total_reads = float(unclassified_total[0] + classified_correct[0] + classified_incorrect[0])

    kraken_summary.loc[key] = pd.Series({'Unclassified': (float(unclassified_total[0])/total_reads)*100, 'Classified_correctly':(float(classified_correct[0])/total_reads)*100, 'Classified_incorrectly':(float(classified_incorrect[0])/total_reads)*100})

# Figure: /home/ae42909/python_figures/synthPotato_kraken1.png
kraken_summary.plot.barh(stacked=True)
plt.xlabel('Percent of reads in group')
plt.title('Kraken classification \nof synthetics RNA-seq library')
plt.show()




############################# WORKING ON IT #############################

#Find frequency of incorrectly classified species in Potato reads
potato_classInc = kraken_all[(kraken_all.Seq_ID.str.contains(library_ID[key][0])) & (kraken_all['Classified'] == 'C') & ~(kraken_all['Tax_ID'] == library_ID[key][1])]

potato_freq = potato_classInc[potato_classInc['Seq_tax'].value_counts()

#Regular expression - find genus, species etc
line = 'd__Eukaryota|k__Viridiplantae|p__Streptophyta|o__Solanales|f__Solanaceae|g__Solanum|s__Solanum_tuberosum'
regex = r'\|[' + convertKey + ']\_\_(.*)?\|'
re.match(regex,line).group()

line = 'd__Eukaryota|k__Viridiplantae|p__Streptophyta|o__Solanales|f__Solanaceae|g__Solanum|s__Solanum_tuberosum'
regex = r'\|[' + convertKey + ']\_\_([^|]*)'
re.findall(regex,line)


#Regular expression - find genus, species etc for each sequence. Make dictionary with level options
tax_levels = {'superkingdom':'d', 'kingdom':'k', 'phylum':'p', 'klass':'c', 'order':'o', 'family':'f','genus':'g','species':'s'}

# Swap keys and values in dictionary
#tax_levels = dict((v,k) for k,v in tax_levels.iteritems())

# Select the levels that you want specifically (can be as many as you like)
tax_need = {k: tax_levels[k] for k in ('genus', 'species')}

# Make a column for each level choosen in tax_need and fill in for each sequence, calculate sensitivity and precision for each level.
for level in tax_need.keys():
    convertKey = tax_need[level]
    regex = r'\|[' + convertKey + ']\_\_([^|]*)'  # Find |letter__anystring and negation based regex [^|]* which means anything but pipe
    kraken_all[level] = kraken_all['Seq_tax'].str.findall(regex)

# Calculate sensitivity and precision fo levels choosen in tax_need. Sensitivity was calculated as the percentage of reads assigned to the correct genus/phylum out of the total number of reads in the input. Precision was calculated as the percentage of reads assigned to the correct genus/phylum out of the number of classified reads, excluding reads classified correctly to a rank above genus/phylum-level \cite{Menzel2015}

kraken_calculation = pd.DataFrame(columns=['Unclassified', 'Sensitivity','Precision'], index=library_ID.keys())

unclassified_total = kraken_all[kraken_all['Classified'] == 'U'].count()[0]

genus = ['Potyvirus', 'Ipomovirus', 'Solanum']
species = []

for keys in library_ID.keys():
    sensitivity_each = kraken_all[kraken_all['genus'] == [library_ID[key][2]])].count()[0]
    sensitivity =

classified_incorrect = kraken_all[(kraken_all.Seq_ID.str.contains(library_ID[key][0])) & (kraken_all['Classified'] == 'C') & ~(kraken_all['Tax_ID'] == library_ID[key][1])].count()

total_reads = float(unclassified_total[0] + classified_correct[0] + classified_incorrect[0])

kraken_summary.loc[key] = pd.Series({'Unclassified': (float(unclassified_total[0])/total_reads)*100, 'Classified_correctly':(float(classified_correct[0])/total_reads)*100, 'Classified_incorrectly':(float(classified_incorrect[0])/total_reads)*100})
