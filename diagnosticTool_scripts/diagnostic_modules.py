# python 2.7.13
import subprocess
import pandas as pd
import numpy as np
import os
from Bio import SeqIO

# Test input data format
def test_format(file1, user_format):
    with open(file1) as myfile:
        small_file = [next(myfile) for x in xrange(8)]

    file_format = "not identified"

    if small_file[0][0] == "@" and small_file[4][0] == "@":
        file_format = "fastq"
    if small_file[0][0] == ">" and small_file[2][0] == ">":
        file_format = "fasta"

    assert (file_format == "fasta") | (file_format == "fastq"), "Cannot proceed with file as it is not in fasta or fastq format."
    assert user_format == file_format, "File has been detected to be in " + file_format + " format rather than " + user_format + " format"


# QC and trim data
def fastqc_trim(out_dir, file1, trim_minlen, threads, adapter_file, file2 = False):

    subprocess.call("fastqc " + file1 + " -o " + out_dir, shell=True)

    trimAdapt_command = "ILLUMINACLIP:" + adapter_file + ":2:30:10 LEADING:20 TRAILING:20 MINLEN:" + trim_minlen

    if file2:
        subprocess.call("fastqc " + file2 + " -o " + out_dir, shell=True)
        subprocess.call("java -jar /mnt/apps/trimmomatic/0.32/trimmomatic.jar PE -threads " + threads + " " + file1 + " " + file2 + " PE_trimmed_data_1P PE_trimmed_data_1U PE_trimmed_data_2P PE_trimmed_data_2U " + trimAdapt_command, shell=True)
        subprocess.call("rm PE_trimmed_data_1U PE_trimmed_data_2U", shell=True)
    else:
        subprocess.call("java -jar /mnt/apps/trimmomatic/0.32/trimmomatic.jar SE -threads " + threads + " " + file1 + " SE_trimmed_data " + trimAdapt_command, shell=True)


# Order and replace sequence IDs with numberic IDs
def rename_seq(trim_file, out_dir, user_format, paired=False):
    # Make a new file with for alphabetically ordered sequence IDs
    ids = sorted(rec.id for rec in SeqIO.parse(trim_file, user_format))
    if not paired:
        id_file = open(out_dir + "ID.txt", 'w')
    else:
        id_file = open(out_dir + "ID" + paired + ".txt", 'w')

    for item in ids:
        id_file.write("%s\n" % item)

    # Reorder the sequence file based on odered IDs
    record_index = SeqIO.index(trim_file, user_format)
    records = (record_index[id] for id in ids)
    SeqIO.write(records, out_dir + "sorted", user_format)

    sorted_file = out_dir + "sorted"
    if not paired:
        renamed_file = open(out_dir + "renamed_file." + user_format, "w")
        paired_str = ""
    else:
        renamed_file = open(out_dir + "renamed_file" + str(paired) + "." + user_format, "w")
        paired_str = "/" + str(paired)

    # Write a new file where seqc sequence ID is replaced with the correct symbol (">" or "@")
    # followed by a number (1::N), and if it's a paired read followd by "/1" or "/2"
    with open(sorted_file, 'r') as in_file:
        seqNum = 1
        if user_format == "fasta":
            symb = ">"
            for line in in_file:
                if line[0] == symb:
                    renamed_file.write(symb + str(seqNum) + paired_str + "\n")
                    seqNum = seqNum + 1
                else:
                    renamed_file.write(line)
        elif user_format == "fastq":
            symb = "@"
            for lineNum, line in enumerate(in_file):
                if lineNum % 4 == 0:
                    renamed_file.write(symb + str(seqNum) + paired_str + "\n")
                    seqNum = seqNum + 1
                else:
                    renamed_file.write(line)
                
        # Once you have the renamed file you don't need the sorted file - delete
        subprocess.call("rm " + out_dir + "sorted", shell=True)


# Kraken classification
def kraken_classify(renamed_file1, threads, user_format, kraken_db, renamed_file2 = False):
    if user_format == "fastq":
        format_switch = " --fastq-input"
    elif user_format == "fasta":
        format_switch = " --fasta-input"
        assert (format_switch == " --fastq-input") | (format_switch == " --fasta-input"), "Incorrect format - check correct format assigned."
        
    kraken_command = "kraken --preload --threads " + threads + " --db " + kraken_db + format_switch
    if renamed_file2:
        kraken_command += " --paired " +  renamed_file1 + user_format + " " + renamed_file2 + user_format + " > kraken_table.txt"
    else:
        kraken_command += " " + renamed_file1 + user_format + " > kraken_table.txt"
    subprocess.call(kraken_command, shell = True)
    subprocess.call("kraken-translate --mpa-format --db " + kraken_db + " kraken_table.txt > kraken_labels.txt", shell = True)

def format_result_table(out_dir, data_table, data_labels, table_colNames):
    label_colNames=["Seq_ID", "Seq_tax"] #Seq_tax: d__superkingdom, k__kingdom, p__phylum, c__class, o__order, f__family, g__genus, s__species
    seq_data = pd.read_csv(out_dir + data_table, sep="\t", header = None, names= table_colNames, index_col=False)
    seq_labelData = pd.read_csv(out_dir + data_labels, sep="\t", header = None, names= label_colNames)
    seq_result = pd.merge(seq_data, seq_labelData, on='Seq_ID', how='outer')

    ncbi_tax = pd.read_csv('/home/ae42909/Scratch/kraken/kraken_analysis/customDatabase/NCBI_taxonomy.csv', sep=",")
    ncbi_slice = ncbi_tax.iloc[:,[1,2]]
    seq_final =  pd.merge(seq_result, ncbi_slice, on='Tax_ID', how='outer')
    seq_final = seq_final.dropna(subset = ['Seq_ID']) # Remove entries from the NCBI table that do not correspond to result

    return seq_final


# Subset viral and unclassified sequences
def seq_reanalysis(kraken_table, kraken_labels, out_dir, user_format, renamed_file1, renamed_file2 = False):
    kraken_colNames = ["kraken_classified", "Seq_ID","Tax_ID", "kraken_length", "kraken_k-mer"]
    kraken_results = format_result_table(out_dir, "kraken_table.txt", "kraken_labels.txt", kraken_colNames)
    kraken_results[['Seq_ID']] = kraken_results[['Seq_ID']].astype(int).astype(str)
    kraken_results.to_csv(out_dir  + 'kraken_FormattedTable.txt', sep='\t', index= False)

    # Delete unformatted kraken table
    subprocess.call("rm " + "kraken_table.txt kraken_labels.txt", shell = True)

    # Make a list of "Seq_ID" column value if sequence is unclassified in "Classified" column or
    #  classified as VRL (virus) in column "Div_ID". This list will be used to determine which sequences
    #  will be further analysed by Kaiju 
    unclassified_IDs = kraken_results.loc[(kraken_results.kraken_classified == 'U'), ['Seq_ID']]
    VRL_IDs = kraken_results.loc[(kraken_results.Div_ID == 'VRL'), ['Seq_ID']]
    reanalyse_IDs = unclassified_IDs['Seq_ID'].tolist() + VRL_IDs['Seq_ID'].tolist()


    # Use biopython to make new fastq files of sequences to be reanalysed. 
    def reanalyse_subset (input_file, output_file, id_list):
        outfile = open(out_dir  + 'reanalyse_ID.txt', 'w')
        print >> outfile, "\n".join(str(i) for i in id_list)
        outfile.close()
        id_file = out_dir  + "reanalyse_ID.txt"

        wanted = set(line.rstrip("\n").split(None,1)[0] for line in open(id_file))
        print "Found %i unique identifiers in %s" % (len(wanted), id_file)

        records = (r for r in SeqIO.parse(out_dir  + input_file, user_format) if r.id in wanted)
        count = SeqIO.write(records, out_dir  + output_file, user_format)
        print "Saved %i records from %s to %s" % (count, input_file, output_file)
        if count < len(wanted):
            print "Warning %i IDs not found in %s" % (len(wanted)-count, input_file)

    if renamed_file2:
        reanalyse_ID1 = [s + "/1" for s in reanalyse_IDs]
        reanalyse_ID2 = [s + "/2" for s in reanalyse_IDs]
        reanalyse_subset(renamed_file1 + user_format, "subset_file1." + user_format, reanalyse_ID1)
        reanalyse_subset(renamed_file2 + user_format , "subset_file2." + user_format, reanalyse_ID2)
        # Delete "renamed" files
        subprocess.call("rm renamed_file1." + user_format + " renamed_file2." + user_format, shell=True)
    else:
        reanalyse_subset(renamed_file1 + user_format, "subset_file1." + user_format, reanalyse_IDs)
        subprocess.call("rm renamed_file." + user_format, shell=True)

    # Deleted "reanalyse.txt"
    subprocess.call("rm reanalyse_ID.txt", shell=True)

# Kaiju classification of subset sequences
def kaiju_classify(subset_file1, user_format, threads, kaiju_nodes, kaiju_fmi, kaiju_names, kaiju_missmatch, kraken_db, subset_file2 = False):
    if subset_file2:
        kaiju_command = "/home/ae42909/Programs/Kaiju/kaiju-v1.5.0-linux-x86_64-static/bin/kaiju -z " + threads + " -t " + kaiju_nodes + " -f " + kaiju_fmi + " -i " + subset_file1 + user_format + " and -j " + subset_file2 + user_format + " -o kaiju_table.txt -x -v -a greedy -e" + kaiju_missmatch
    else:
        kaiju_command = "/home/ae42909/Programs/Kaiju/kaiju-v1.5.0-linux-x86_64-static/bin/kaiju -z " + threads + " -t " + kaiju_nodes + " -f " + kaiju_fmi + " -i "  + subset_file1 + user_format + " -o kaiju_table.txt -x -v -a greedy -e" + kaiju_missmatch

    subprocess.call(kaiju_command, shell = True)
    subprocess.call("kraken-translate --mpa-format --db " + kraken_db + " " + "kaiju_table.txt > kaiju_labels.txt", shell = True)

    # Delete subset files
    subprocess.call("rm subset_file1.fastq", shell=True)
    if subset_file2:
        subprocess.call("rm subset_file2.fastq", shell=True)


# Import kraken and kaiju results, merge and summarise
def result_analysis(out_dir, kraken_Formattedtable, kaiju_table, kaiju_label, file1_IDs, file2_IDs = False):
    # Import kraken table
    kraken_results = pd.read_csv(out_dir + kraken_Formattedtable, header = 0, sep='\t',
                                dtype={"Classified":str, "Seq_ID": int,"Tax_ID": float, "length": float, "k-mer classification":str, "Seq_tax": str, "Div_ID": str})

    kaiju_colNames =["kaiju_classified", "Seq_ID","Tax_ID", "kaiju_lenBest", "kaiju_tax_AN","kaiju_accession", "kaiju_fragment"]
    kaiju_results = format_result_table(out_dir, "kaiju_table.txt", "kaiju_labels.txt", kaiju_colNames) 

    # Merge Kaiju and kraken results, sort, reindex and rename two columns
    kraiju = pd.merge(kraken_results, kaiju_results, on='Seq_ID', how='outer')
    kraiju = kraiju.sort_values(['Seq_ID'])
    kraiju = kraiju.reset_index(drop=True)
    kraiju.rename(columns={'Div_ID_x':'kraken_div_ID', 'Div_ID_y':'kaiju_div_ID', 'Tax_ID_x':'kraken_tax_ID', 'Tax_ID_y':'kaiju_tax_ID', "Seq_tax_x":"kraken_seq_tax", "Seq_tax_y":"kaiju_seq_tax"}, inplace=True)

    # Delete unnecessary files
    subprocess.call("rm kaiju_table.txt kaiju_labels.txt kraken_FormattedTable.txt ", shell=True)

    # Get results where they are in agreement
    kraiju['combined_result'] = kraiju.kraken_tax_ID[kraiju['kraken_tax_ID'] == kraiju['kaiju_tax_ID']]

    # Replace numeric IDs with real IDs
    num_ids = list(kraiju.Seq_ID)

    if file2_IDs:
        with open('ID.txt', 'w') as ID_file, open('ID1.txt') as f1, open('ID2.txt') as f2:
            for line1, line2 in zip(f1, f2):
                ID_file.write("{} {}\n".format(line1.rstrip(), line2.rstrip()))
        subprocess.call("rm ID1.txt ID2.txt", shell=True)

    real_ids = []
    with open(out_dir + "ID.txt", 'r') as f:
        for line in f:
            real_ids.append(line.rstrip())

    id_dict = dict(zip(num_ids, real_ids))
    seq_col = kraiju["Seq_ID"]
    id_col = seq_col.map(id_dict)
    kraiju["Seq_ID"] = id_col

    kraiju.to_csv(out_dir  + 'kraiju_FormattedTable.txt', sep='\t', index= False)

    # Separate reads which are classified as VRL, make a table with all identified viruses and count number of intances for each
    kraiju_vrl = kraiju[(kraiju.kraken_div_ID == 'VRL') | (kraiju.kaiju_div_ID == 'VRL')]

# Create a summary table - Species, genus, phylum, number of reads classified kraken/kaiju/agreed by both 

    def summary_table(kraiju_data):
        kraken_class = dict(kraiju_data['kraken_tax_ID'].value_counts())
        kraken_levels = pd.Series(kraiju_data.kraken_seq_tax.values,index=kraiju_data.kraken_tax_ID).to_dict()
        kaiju_class = dict(kraiju_data['kaiju_tax_ID'].value_counts())
        kaiju_levels = pd.Series(kraiju_data.kaiju_seq_tax.values,index=kraiju_data.kaiju_tax_ID).to_dict()
        combined_class = dict(kraiju_data['combined_result'].value_counts())

        levels_dict = kraken_levels.copy()
        levels_dict.update(kaiju_levels)
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

        table_summary = pd.DataFrame(columns=['Species', 'Tax_ID', 'kraken', 'kaiju','combined', 'associated_seq'])
        table_summary['Tax_ID'] =  map(int, levels_tax.keys())
        table_summary['kraken'] = table_summary['Tax_ID'].map(kraken_class)
        table_summary['kaiju'] = table_summary['Tax_ID'].map(kaiju_class)
        table_summary['combined'] = table_summary['Tax_ID'].map(combined_class)
        table_summary['Species'] = table_summary['Tax_ID'].map(species_dict)
        table_summary.to_csv('virus_table.txt', sep = '\t', index = False)

        # table_summary = pd.read_csv(out_dir + "virus_table.txt", header = 0, sep='\t')

        # Make a table of sequences of each species (including LCA potential sequences)
        def species_tables(tax_id):
            sp_table = kraiju_data.loc[kraiju_data['kraken_tax_ID'].isin(associated_tax[tax_id]) | kraiju_data['kaiju_tax_ID'].isin(associated_tax[tax_id])]
            return sp_table

        for index, row in table_summary.iterrows():
            sp_tax = species_tables(row.Tax_ID)
            sp_tax.to_csv("sp_" + str(row.Tax_ID) + ".txt", sep = '\t', index = False)

    summary_table(kraiju_vrl)

# def retrive_viralSequences(summary_table):

