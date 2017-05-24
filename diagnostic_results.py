
#Analysing data from the diagnostic_pipeline.sh script
#python 2.7
import os
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use('ggplot')
from Bio import SeqIO

directory = os.environ["OUTDIR"]

# Name of the data to be analysed

result_data='kraken_aa_results'
label_data= 'kraken_nt_labels'

def kraken_analysis(result_data, label_data):
    kraken_data = directory + '/' + result_data
    kraken_labelData = directory + '/' + label_data

    result_names=["Classified", "Seq_ID","Tax_ID", "length", "k-mer classification" ]
    label_names=["Seq_ID", "Seq_tax"] #(superkingdom, kingdom, phylum, class, order, family, genus, species)

    # Read tables in with pandas
    kraken_result = pd.read_csv(result_data, sep="\t", header = None, names= result_names)
    kraken_label = pd.read_csv(label_data, sep="\t", header = None, names= label_names)

    # Merge result tables for Kraken - to include taxonomc names to each sequence that was classified - http://chrisalbon.com/python/pandas_join_merge_dataframe.html
    kraken_table = pd.merge(kraken_result, kraken_label, on='Seq_ID', how='outer')

### To be able to identify what broad category (plant, virus bacteria, etc) a sequence has been categorised into and NCBI taxonomy table was made with the broad categories next to each taxID and merged with the kraken results table

# Possible categories for NCBI_taxonomy table abreviations and long names:
Div_tax = {'UNA':'Unannotated', 'BCT':'Bacteria', 'ENV':'Environmental samples', 'SYN':'Synthetic', 'PLN':'Plants', 'INV':'Invertebrates', 'VRT':'Other vertebrates', 'MAM':'Other mammals', 'PRI':'Primates', 'ROD':'Rodents', 'VRL':'Viruses', 'PHG':'Phages'}

# Read in NCBI table - made using script /home/scripts_inProcess/NCBI_taxonomy.py
ncbi_tax = pd.read_csv('/home/ae42909/Scratch/kraken/kraken_analysis/customDatabase/NCBI_taxonomy.csv', sep=",")
ncbi_slice = ncbi_tax.iloc[:,[1,2]]

# Merge tables NCBI taxonomy categories with kraken data
kraken_all =  pd.merge(kraken_all, ncbi_slice, on='Tax_ID', how='outer')

# Remove entries from the NCBI table that do not correspond to any in the kraken results table
kraken_all = kraken_all.dropna(subset = ['Seq_ID', 'Classified'])

# Make a list of "Seq_ID" column value if sequence is unclassified in "Classified" column or classified as VRL (virus) in column "Div_ID". This list will be used to determine which sequences will be further analysed by Kaiju and retranslation.py

unclassified_IDs = kraken_all.loc[(kraken_all.Classified == 'U'), ['Seq_ID']]
VRL_IDs = kraken_all.loc[(kraken_all.Div_ID == 'VRL'), ['Seq_ID']]
reanalyse_IDs = unclassified_IDs['Seq_ID'].tolist() + VRL_IDs['Seq_ID'].tolist()

# Export the list seq_reanalyse so that each line is a Seq_ID
outfile = open("test_reanalyse_IDs1.txt", "w")
print >> outfile, "\n".join(str(i) for i in test)
outfile.close()


# Use biopython to make new fastq files of sequences to be reanalysed. The IDs from Kraken ('Seq_ID') don't correspond to the IDs when you parse data with SeqIO, they are shortened (don't have the "/1" and "/2" that identifies them as paired end 1 and 2). It is easy to add the "/1" (as is done in the section "Test ID" below), but it could be differennt for different paired end files so need to find a better solution.

####### Test ID: Having a problem with te IDs (they need to be exact so tha seqio can identify them in "records"). This test is to check if indeed taht is the case (only requires a "/1" or "/2"). The output of this test ("subset_Potato_withViruses_1/2.fastq") was fed into Kaiju for protein analysis and it worked.

# Test for first pairend fastq file:
#test1 = [s + "/1" for s in reanalyse_IDs]

#outfile = open("test_reanalyse_IDs1.txt", "w")
#print >> outfile, "\n".join(str(i) for i in test1)
#outfile.close()

#id_file = "test_reanalyse_IDs1.txt"
#wanted = set(line.rstrip("\n").split(None,1)[0] for line in open(id_file))
#print "Found %i unique identifiers in %s" % (len(wanted), id_file)

#input_file = "Potato_withViruses_1.fastq"
#output_file = "subset_Potato_withViruses_1.fastq"
#records = (r for r in SeqIO.parse(input_file, "fastq") if r.id in wanted)
#count = SeqIO.write(records, output_file, "fastq")
#print "Saved %i records from %s to %s" % (count, input_file, output_file)
#if count < len(wanted):
#    print "Warning %i IDs not found in %s" % (len(wanted)-count, input_file)

# Test for second pairend fastq file:
#test2 = [s + "/2" for s in reanalyse_IDs]

#outfile = open("test_reanalyse_IDs2.txt", "w")
#print >> outfile, "\n".join(str(i) for i in test2)
#outfile.close()
    
#id_file = "test_reanalyse_IDs2.txt"
#wanted = set(line.rstrip("\n").split(None,1)[0] for line in open(id_file))
#print "Found %i unique identifiers in %s" % (len(wanted), id_file)

#input_file = "Potato_withViruses_2.fastq"
#output_file = "subset_Potato_withViruses_2.fastq"
#records = (r for r in SeqIO.parse(input_file, "fastq") if r.id in wanted)
#count = SeqIO.write(records, output_file, "fastq")
#print "Saved %i records from %s to %s" % (count, input_file, output_file)
#if count < len(wanted):
#    print "Warning %i IDs not found in %s" % (len(wanted)-count, input_file)


# Find what the SeqIO identifiers have in common - could find extra pattern to add to the Kraken IDs

def longestSubstringFinder(string1, string2):
    answer = ""
    len1, len2 = len(string1), len(string2)
    for i in range(len1):
        match = ""
        for j in range(len2):
            if (i + j < len1 and string1[i + j] == string2[j]):
                match += string2[j]
            else:
                if (len(match) > len(answer)): answer = match
                match = ""
    return answer

print longestSubstringFinder("apple pie available", "apple pies")
print longestSubstringFinder("apples", "appleses")
print longestSubstringFinder("bapples", "cappleses")

input_file = "Potato_withViruses_1.fastq"

all_IDs = []
for record in SeqIO.parse(input_file, "fastq"):
    all_IDs.append(record.id)



reanalyse_SeqIO_IDs = []
for ids in reanalyse_IDs:
    possible_IDs = [s for s in all_IDs if ids in s]
    if len(possible_IDs) > 1:
        reanalyse_SeqIO_IDs.append(min(possible_IDs, key=len))
    else:
        reanalyse_SeqIO_IDs.append(possible_IDs[0])
        

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



############ Merge Kaiju and kraken results
os.chdir('/home/ae42909/Scratch/kaiju')

# Name of the data to be analysed
#result_data='SynthPotato_ouput'
result_data='subsetSynthPotato_ouput'
result_names=["Classified", "Seq_ID","Tax_ID", "length_bestmatch", "Tax_AN","accession_multiple", "Fragment" ]

kaiju_result = pd.read_csv(result_data, sep="\t", header = None, names= result_names)

kaiju_small = kaiju_result[['Classified','Seq_ID','Tax_ID']]

# Merge Kraken and Kaiju results
Ks = pd.merge(kraken_all, kaiju_small, on='Seq_ID', how='outer')




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
    
