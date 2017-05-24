#Analysing data frm kaiju output
#python 2.7
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use('ggplot')

os.chdir('/home/ae42909/Scratch/kaiju')

# Name of the data to be analysed
#result_data='SynthPotato_ouput'
result_data='subsetSynthPotato_ouput'
result_names=["Classified", "Seq_ID","Tax_ID", "length_bestmatch", "Tax_AN","accession_multiple", "Fragment" ]

kaiju_result = pd.read_csv(result_data, sep="\t", header = None, names= result_names)

# Reads classified in total for each ID, correctly, incorrectly and unclassified.
library_ID =  {}
library_ID['Tabacco_etch']=['SRR3466597', 12227]
library_ID['UBSV']=['ERR996011', 946046]
library_ID['BSV']=['ERR996013', 137758]
library_ID['Potato']=['SRR1207289', 4113]

kaiju_summary = pd.DataFrame(columns=['Unclassified','Classified_correctly','Classified_incorrectly'], index=library_ID.keys())

for key in library_ID.keys():
    unclassified_total = kaiju_result[(kaiju_result.Seq_ID.str.contains(library_ID[key][0])) & (kaiju_result['Classified'] == 'U')].count()

    classified_correct = kaiju_result[(kaiju_result.Seq_ID.str.contains(library_ID[key][0])) & (kaiju_result['Classified'] == 'C') & (kaiju_result['Tax_ID'] == library_ID[key][1])].count()

    classified_incorrect = kaiju_result[(kaiju_result.Seq_ID.str.contains(library_ID[key][0])) & (kaiju_result['Classified'] == 'C') & ~(kaiju_result['Tax_ID'] == library_ID[key][1])].count()

    total_reads = float(unclassified_total[0] + classified_correct[0] + classified_incorrect[0])

    kaiju_summary.loc[key] = pd.Series({'Unclassified': (float(unclassified_total[0])/total_reads)*100, 'Classified_correctly':(float(classified_correct[0])/total_reads)*100, 'Classified_incorrectly':(float(classified_incorrect[0])/total_reads)*100})

kaiju_summary.plot.barh(stacked=True)
plt.xlabel('Percent of reads in group')
plt.title('Kaiju classification (5 missmatches) \nof synthetics RNA-seq library')
plt.show()

# Figure: /home/ae42909/python_figures/synthPotato_kaiju1.png
