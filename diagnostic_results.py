#Analysing data from the diagnostic_pipeline.sh script
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
#directory = '/home/ae42909/Scratch/synthPotato_pipeline/'
directory = '/home/ae42909/Scratch/fullPipeline_krakenDB_viral/'

# Import initial results from sequence_reanalysis.py
# kraken_nt_table = pd.read_csv(directory + 'kraken_nt_table', header = 0, sep='\t')
kraken_nt_table = pd.read_csv(directory + 'kraken_FormatedTable.txt', header = 0, sep='\t')

# Import Kaiju results and add ncbi div_ids
kaiju_data= directory + 'kaiju_ouput'
kaiju_names=["kaiju_classified", "Seq_ID","Tax_ID", "length_bestmatch", "Tax_AN","accession_multiple", "Fragment" ]
kaiju_result = pd.read_csv(kaiju_data, sep="\t", header = None, names= kaiju_names)
#kaiju_small = kaiju_result[['kaiju_classified', 'Seq_ID','kaiju_TaxID']]

ncbi_tax = pd.read_csv('/home/ae42909/Scratch/kraken/kraken_analysis/customDatabase/NCBI_taxonomy.csv', sep=",")
ncbi_slice = ncbi_tax.iloc[:,[1,2]]
kaiju_final =  pd.merge(kaiju_result, ncbi_slice, on='Tax_ID', how='outer')
kaiju_final = kaiju_final.dropna(subset = ['Seq_ID']) # Remove entries from the NCBI table that do not correspond to result

# Merge Kaiju and kraken results, sort, reindex and rename two columns
kraiju = pd.merge(kraken_nt_table, kaiju_final, on='Seq_ID', how='outer')
kraiju = kraiju.sort_values(['Seq_ID'])
kraiju = kraiju.reset_index(drop=True)
kraiju.rename(columns={'Div_ID_x':'kraken_Div_ID', 'Div_ID_y':'kaiju_Div_ID', 'Tax_ID_x':'kraken_Tax_ID', 'Tax_ID_y':'kaiju_Tax_ID'}, inplace=True)

# Get results where they are in agreement
kraiju['combined_result'] = None
kraiju['combined_result'] = kraiju.kraken_Tax_ID[kraiju['kraken_Tax_ID'] == kraiju['kaiju_Tax_ID']]
kraiju.to_csv(directory  + 'kraiju_table', sep='\t', index= False)

# Separate reads which are classified as VRL, make a table with all identified viruses and count number of intances for each
kraiju_vrl = kraiju[(kraiju.kraken_Div_ID == 'VRL') | (kraiju.kaiju_Div_ID == 'VRL')]

# Create a summary table 

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

vrl_summary = summary_table(kraiju_vrl)
vrl_summary.to_csv('virus table.txt', sep = '\t', index = False)

# Add colums for backtranslate results
# kraiju["Backtranslate_classification"] = np.nan
# kraiju["Backtranslate_TaxID"] = np.nan
# kraiju["Backtranslate_k-mer"] = np.nan

#Import results from aa kraken analysis (backtranslated sequences), add labels and NCBi taxa

# def kraken_analysis(result_data, label_data):
#     result_names=["Classified", "Seq_ID","Tax_ID", "length", "k-mer" ]
#     label_names=["Seq_ID", "Seq_tax"] #(superkingdom, kingdom, phylum, class, order, family, genus, species)

#     # Read tables in with pandas
#     kraken_result = pd.read_csv(directory + result_data, sep="\t", header = None, names= result_names)
#     kraken_label = pd.read_csv(directory + label_data, sep="\t", header = None, names= label_names)

#     # Merge result tables for Kraken - to include taxonomc names to each sequence that was classified - http://chrisalbon.com/python/pandas_join_merge_dataframe.html
#     kraken_result = pd.merge(kraken_result, kraken_label, on='Seq_ID', how='outer')

#     # Add taxa information from NCBI
#     Div_tax = {'UNA':'Unannotated', 'BCT':'Bacteria', 'ENV':'Environmental samples', 'SYN':'Synthetic', 'PLN':'Plants', 'INV':'Invertebrates', 'VRT':'Other vertebrates', 'MAM':'Other mammals', 'PRI':'Primates', 'ROD':'Rodents', 'VRL':'Viruses', 'PHG':'Phages'}

#     # Read in NCBI table - made using script /home/scripts_inProcess/NCBI_taxonomy.py
#     ncbi_tax = pd.read_csv('/home/ae42909/Scratch/kraken/kraken_analysis/customDatabase/NCBI_taxonomy.csv', sep=",")
#     ncbi_slice = ncbi_tax.iloc[:,[1,2]]

#     # Merge tables NCBI taxonomy categories with kraken data
#     kraken_all =  pd.merge(kraken_result, ncbi_slice, on='Tax_ID', how='outer')

#     # Remove entries from the NCBI table that do not correspond to any in the kraken results table, convert sequence IDs back to str and export results table to append further data
#     kraken_all = kraken_all.dropna(subset = ['Seq_ID'])
#     kraken_all[['Seq_ID']] = kraken_all[['Seq_ID']].astype(str) # Not the same as the def in 'sequence_reanalysis.py
#     return kraken_all

# kraken_aa_table = kraken_analysis('kraken_aa_results', 'kraken_aa_labels')

# ### Find which of the 6 frames is most likely
# # Make a new column for Seq_ID and frame (1:6) separetely

# frame_df = pd.DataFrame(kraken_aa_table.Seq_ID.str.split('_',1).tolist(),columns = ['Seq_ID','aa_frame'])

#  # remove colmn with Seq_ID + frame number

# kraken_frame = pd.concat([kraken_aa_table.drop('Seq_ID', axis=1), frame_df], axis=1)
# kraken_frame['Seq_ID'] = kraken_frame['Seq_ID'].apply(int)
# kraken_frame = kraken_frame.sort_values(['Seq_ID', 'aa_frame'])
# kraken_frame = kraken_frame.reset_index(drop=True) # reset the index as all classifed were out of order (so index was confusing)



# t0 = time.time()
# #for seqIDs in range(0,len(kraken_frame)/6):
# for i in range(0,200):
#     start = i*6
#     end  = start + 6
    
#     df_seqID = kraken_frame.iloc[start:end]
#     class_seqID = np.sum(class_arr[start:end])
#     table_index = np.where(ix_arr == int(kraken_frame.Seq_ID[start]))[0][0]
#     if class_seqID < 1.0:
#         kraiju.ix[table_index, 'Backtranslate_classification'] = 'U'
#     if class_seqID == 1.0:
#         kraiju.ix[table_index, 'Backtranslate_classification'] = 'C'
#         kraiju.ix[table_index, 'Backtranslate_TaxID'] = df_seqID[df_seqID['Classified'] == 'C'].Tax_ID.item()
#         kraiju.ix[table_index, 'Backtranslate_k-mer'] = df_seqID[df_seqID['Classified'] == 'C'].Seq_tax.item()
#     if class_seqID > 1.0:
#         kraiju.ix[table_index, 'Backtranslate_classification'] = 'A'
# t1 = time.time()
# print t1-t0



# def backtranslate_analysis (kraiju_data, frame_data):
#     kraiju_reanalysis = kraiju_data.loc[(kraiju_data.Classified == 'U')].append(kraiju_data.loc[(kraiju_data.Div_ID == 'VRL')])
#     kraiju_left = kraiju_data.loc[(kraiju_data.Classified != 'U') & (kraiju_data.Div_ID != 'VRL')]

#     class_arr = np.asarray(frame_data['Classified'] == 'C')
#     ix_arr = np.asarray(kraiju_data['Seq_ID'])
#     df_together = pd.DataFrame(data={'kraiju_ix':[],
#                              'Backtranslate_classification':[],
#                              'Backtranslate_TaxID':[],
#                              'Backtranslate_k-mer':[]})

#     #for i in range(0,len(frame_data)/6):
#     for i in range(0,2000):
#         start = i*6
#         end  = start + 6

#         df_seqID = frame_data.iloc[start:end]
#         class_seqID = np.sum(class_arr[start:end])
#         table_index = np.where(ix_arr== int(frame_data.Seq_ID[start]))[0][0]


#         if class_seqID < 1:
#             df_fill = pd.DataFrame(data={'kraiju_ix':table_index,
#                                          'Backtranslate_classification':['U']})
#         if class_seqID == 1:
#             df_fill = pd.DataFrame(data={'kraiju_ix':table_index,
#                                  'Backtranslate_classification':['C'],
#                                  'Backtranslate_TaxID':[df_seqID[df_seqID['Classified'] == 'C'].Tax_ID.item()],
#                                  'Backtranslate_k-mer':[df_seqID[df_seqID['Classified'] == 'C'].Seq_tax.item()]})

#         elif class_seqID > 1:
#             df_fill = pd.DataFrame(data={'kraiju_ix':table_index,
#                                  'Backtranslate_classification':['A']})


#         df_together = df_together.append(df_fill)

#     df_together = df_together.set_index('kraiju_ix')
#     kraiju_reanalysis[['Backtranslate_classification', 'Backtranslate_TaxID', 'Backtranslate_k-mer']] = kraiju_reanalysis[['Backtranslate_classification', 'Backtranslate_TaxID', 'Backtranslate_k-mer']].fillna(df_together)

#     kraijate = pd.concat([kraiju_reanalysis, kraiju_left]).sort_values(['Seq_ID'])
#     return kraijate

# %mprun kraijate = backtranslate_analysis(kraiju, kraken_frame)

# kraijate.to_csv(directory  + 'kraijate_table', sep='\t', index= False)



# Get the real sequence ID    

def reanalyse_data (list_sequenceData, list_subsetDataNames):
    
    id_file = "reanalyse_SeqIO_IDs.txt"
    wanted = set(line.rstrip("\n").split(None,1)[0] for line in open(id_file))
    print "Found %i unique identifiers in %s" % (len(wanted), id_file)
    
    for libraries in list_sequenceData:
        input_file = list_sequenceData[libraries]
        output_file = list_subsetDataNames[libraries]
        records = (r for r in SeqIO.parse(input_file, "fastq") if r.id in wanted)
        count = SeqIO.write(records, output_file, "fastq")
        print "Saved %i records from %s to %s" % (count, input_file, output_file)
        if count < len(wanted):
            print "Warning %i IDs not found in %s" % (len(wanted)-count, input_file)

sequenceData = ["/home/ae42909/synthetic_RNAseq/mappingRNAseq/concatenated_fastaFiles/Potato_withViruses_1.fastq", "/home/ae42909/synthetic_RNAseq/mappingRNAseq/concatenated_fastaFiles/Potato_withViruses_2.fastq"]
subsetDataNames = ["subset_Potato_withViruses_1.fastq", "subset_Potato_withViruses_2.fastq"]








#### Testing synthetic data ####

### When concatinating libraries to make synthetic data, each had an ID that was kept for identification (first value dictionary). Second value in dict is the correct taxID from NCBI - to determine if the sequences were classified correctlyfor the correct species. The Third value in dict is the correct genus of the species 
library_ID =  {}
library_ID['Tabacco_etch']=['SRR3466597', 12227, 'Potyvirus']
library_ID['UBSV']=['ERR996011', 946046, 'Ipomovirus']
library_ID['BSV']=['ERR996013', 137758, 'Ipomovirus']
library_ID['Potato']=['SRR1207289', 4113, 'Solanum']

# Convert 'Tax_ID' column to str
kraken_all['Tax_ID'] = kraken_all['Tax_ID'].astype(str)

# Summarise kraken data: percent of sequences (from each key in library_ID) which are unclassified,  classified correctly or classified incorrectly at the species level
kraken_summary = pd.DataFrame(columns=['Unclassified', 'Classified_correctly','Classified_incorrectly'], index=library_ID.keys())
 
for key in library_ID.keys():
    unclassified_total = kraken_all[(kraken_all.Seq_ID.contains(library_ID[key][0])) & (kraken_all['Classified'] == 'U')].count()

    classified_correct = kraken_all[(kraken_all.Seq_ID.str.contains(library_ID[key][0])) & (kraken_all['Classified'] == 'C') & (kraken_all['Tax_ID'] == library_ID[key][1])].count()

    classified_incorrect = kraken_all[(kraken_all.Seq_ID.str.contains(library_ID[key][0])) & (kraken_all['Classified'] == 'C') & ~(kraken_all['Tax_ID'] == library_ID[key][1])].count()

    total_reads = float(unclassified_total[0] + classified_correct[0] + classified_incorrect[0])

    kraken_summary.loc[key] = pd.Series({'Unclassified': (float(unclassified_total[0])/total_reads)*100, 'Classified_correctly':(float(classified_correct[0])/total_reads)*100, 'Classified_incorrectly':(float(classified_incorrect[0])/total_reads)*100})

# Figure: /home/ae42909/python_figures/synthPotato_kraken1.png
kraken_summary.plot.barh(stacked=True)
plt.xlabel('Percent of reads in group')
plt.title('Kraken classification \nof synthetics RNA-seq library')
plt.show()

library_ID =  {}
library_ID['Tabacco_etch']=['SRR3466597', 12227, 'Potyvirus']
library_ID['UBSV']=['ERR996011', 946046, 'Ipomovirus']
library_ID['BSV']=['ERR996013', 137758, 'Ipomovirus']
library_ID['Potato']=['SRR1207289', 4113, 'Solanum']

def percents (obj):
    p = float(obj[0])/float(total[0])*100
    return round(p,2)
kraken_summary = pd.DataFrame(columns=['Unclassified', 'Classified_correctly','Classified_incorrectly'], index=library_ID.keys())
kaiju_summary = pd.DataFrame(columns=['Unclassified', 'Classified_correctly','Classified_incorrectly'], index=library_ID.keys())
together_summary = pd.DataFrame(columns=['Unclassified', 'Classified_correctly','Classified_incorrectly'], index=library_ID.keys())

for key in library_ID.keys():
    subset_key = Ks[(Ks.Seq_ID.str.contains(library_ID[key][0]))]
    total = subset_key.count()

    kraken_correct = subset_key[(subset_key['Classified_x'] == 'C') & (subset_key['Tax_ID_x'] == library_ID[key][1])].count()
    kraken_incorr = subset_key[(subset_key['Classified_x'] == 'C') & ~(subset_key['Tax_ID_x'] == library_ID[key][1])].count()
    kraken_unclass = subset_key[(subset_key['Classified_x'] == 'U')].count()
    
    kaiju_correct = subset_key[(subset_key['Classified_y'] == 'C') & (subset_key['Tax_ID_y'] == library_ID[key][1])].count()
    kaiju_incorr = subset_key[(subset_key['Classified_y'] == 'C') & ~(subset_key['Tax_ID_y'] == library_ID[key][1])].count()
    kaiju_unclass = subset_key[(subset_key['Classified_y'] == 'U')].count() 

#    either_correct = subset_key[(subset_key['Tax_ID_x'] == library_ID[key][1]) | (subset_key['Tax_ID_y'] == library_ID[key][1])].count()

    both_classified = subset_key[(subset_key['Classified_x'] == 'C') & (subset_key['Classified_y'] == 'C')]
    both_correct = both_classified[(both_classified['Tax_ID_x'] == library_ID[key][1]) & (both_classified['Tax_ID_y'] == library_ID[key][1])].count()
    both_incorrect = both_classified[(both_classified['Classified_x'] == 'C') & (['Classified_y'] == 'C') & ~(both_classified['Tax_ID_x'] == library_ID[key][1]) & ~(both_classified['Tax_ID_y'] == library_ID[key][1])].count()
    both_unclassified = subset_key[(subset_key['Classified_x'] == 'U') & (subset_key['Classified_y'] == 'U')].count()

    kraken_summary.loc[key] = pd.Series({'Unclassified':percents(kraken_unclass), 'Classified_correctly':percents(kraken_correct), 'Classified_incorrectly':percents(kraken_incorr)})
    kaiju_summary.loc[key] = pd.Series({'Unclassified':percents(kaiju_unclass) , 'Classified_correctly': percents(kaiju_correct), 'Classified_incorrectly': percents(kaiju_incorr)})
    together_summary.loc[key] = pd.Series({'Unclassified':percents(both_unclassified), 'Classified_correctly':percents(both_correct), 'Classified_incorrectly':percents(both_incorrect)})


# Figure: /home/ae42909/python_figures/synthPotato_kraken1.png
kraken_summary.plot.barh(stacked=True)
plt.xlabel('Percent of reads in group')
plt.title('Kraken classification \nof synthetics RNA-seq library')
plt.show()

kaiju_summary.plot.barh(stacked=True)
plt.xlabel('Percent of reads in group')
plt.title('Kaiju classification \nof synthetics RNA-seq library')
plt.show()

together_summary.plot.barh(stacked=True)
plt.xlabel('Percent of reads in group')
plt.title('Kraken and Kaiju combined classification \nof synthetics RNA-seq library')
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


###### Making a sequence identifier library after these have been replaced with numberscontaining the previous identifiers and the new

with open("/home/ae42909/synthetic_RNAseq/fastq_IDs") as f:
    ID_list = f.read().splitlines()

ID_dict = {}  

for i in range(0,len(ID_list)):
    ID_dict[i + 1] =ID_list[i]

example = pd.DataFrame({'col2': {0: 'a', 1: 2, 2: np.nan}, 'col1': {0: 'w', 1: 1, 2: 2}})

example.replace({"col1": ID_dict})

# Copy the sequence ID numbers to another column and replace with second sequence ID

# Import replaced sequence IDs
id1_values = open("ID1", "r").read().split('\n')
id2_values = open("ID2", "r").read().split('\n')
id_keys = list(range(1,len(id1_values)))
id1_dic = dict(zip(id_keys, id1_values))
id2_dic = dict(zip(id_keys, id2_values))

reanalyse_ID1 = map(id1_dic.get, reanalyse_IDs)
reanalyse_ID2 = map(id2_dic.get, reanalyse_IDs)
    
