import pandas as pd
import re

mash_data = '/home/ae42909/Scratch/mash_test/distances.tab'

column_labels=['Reference_ID', 'Query-ID', 'Mash_distance', 'P_value', 'Matching_hashes']

mash_result = pd.read_csv(mash_data, sep="\t", header = None, names=column_labels)

# Get the matching hash number alone
mash_result['Hits'] = mash_result.Matching_hashes.str.split('/').str[0].astype(float)

total_hashes = float(mash_result['Matching_hashes'][0].split('/')[1])

# Extract the names of the species in mash_relevant
def find_TaxID(Ref_ID):
#    result = re.findall('(.*)\.fna$', Ref_ID.split('.-')[-1])
    result = re.findall('-AC-([0-9]*)-', Ref_ID)
    if not result:
        result = re.findall('-NC-([0-9]*)-', Ref_ID)
    return "".join(result)

mash_result['Tax_ID'] = mash_result['Reference_ID'].apply(find_TaxID)

# Take results which have a p-value of 0.05 or lower
mash_relevant = mash_result[mash_result.P_value <= 0.05]

# Check which organisms are being detected
library_ID =  {}
library_ID['Tabacco_etch']=['SRR3466597', '12227', 'Potyvirus']
library_ID['UBSV']=['ERR996011', '946046', 'Ipomovirus']
library_ID['BSV']=['ERR996013', '137758', 'Ipomovirus']
library_ID['Potato']=['SRR1207289', '4113', 'Solanum']

for keys in library_ID:
    if library_ID[keys][1] in list(mash_relevant.Tax_ID):
        print keys



