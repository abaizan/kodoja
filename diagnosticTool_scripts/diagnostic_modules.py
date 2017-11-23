# python 2.7.13
import subprocess
import pandas as pd
from Bio import SeqIO
import random
import os
import pickle
from math import isnan

ncbi_file = '/home/ae42909/Scratch/kraken/kraken_analysis/customDatabase/NCBI_taxonomy.csv'


def check_path(dirs):
    """Check if directory path has '/' at the end.

    Return value is either '/' or empty string ''.
    """
    if dirs[-1] != "/":
        return "/"
    else:
        return ""


def test_format(file1, user_format):
    """Check if data is in the fasta or fastq format and
    assert the user has specified the correct format for
    the data provided.

    Return an assert stament and stop or continue.
    """
    with open(file1) as myfile:
        small_file = [next(myfile) for x in xrange(8)]

    file_format = "not identified"

    if small_file[0][0] == "@" and small_file[4][0] == "@":
        file_format = "fastq"
    if small_file[0][0] == ">":
        file_format = "fasta"

    assert (file_format == "fasta") | (file_format == "fastq"), \
        "Cannot proceed with file as it is not in fasta or fastq format."
    assert user_format == file_format, \
        "File has been detected to be in " + file_format + \
        " format rather than " + user_format + " format."


def check_file(file1, out_dir, user_format, file2=False):
    """Rename sequnce ids for SE or PE files to ensure
    consistency between kraken and kaiju (which modify
    id names). Create dictionaries containing real IDs and 
    renamed version and pickle. If data is PE, assert 
    paired files have the same number of entries and if 
    the paired reads are mached by choosing a random 
    entry from the first list of ids (minus meatdata) and 
    the same entry for the second list of ids (can have 
    one character different  as could be named 1 or 2).
    """
    def str_overlap(str1,str2):
        """ Determine number of matching characters between
        two strings.

        Returns int with number of matching characters
        """

        count = 0
        for i in range(min(len(str1), len(str2))):
            if str1[i] == str2[i]:
                count = count + 1
        return count

    def rename_seqIDs(input_file, out_dir, user_format, paired=False):
        """ Write a new file where each sequence ID is replaced with 
        the the first character (">" or "@") followed by a number 
        (1::N), and if it's a paired read followd by "/1" or "/2"
        (called 'renamed_file'), and a dictionary composed of 
        dict[newID] = oldID.

        """
        output_file = out_dir + "renamed_file_"
        if paired == 2:
            output_file += str(paired) + "." + user_format
        else:
            output_file += "1." + user_format
        id_dict = {}
        with open(input_file, 'r') as in_file, open(output_file, 'w') as out_file:
            seqNum = 1
            seqid_line = (4 if user_format == 'fastq' else 2)

            for lineNum, line in enumerate(in_file):
                if lineNum %  seqid_line == 0:
                    new_line = line[0] + str(seqNum)
                    id_dict[seqNum] = line[1:].strip()
                    if paired:
                        new_line += '/' + str(paired)
                    new_line += "\n"
                    out_file.write(new_line)
                    seqNum += 1
                else:
                    out_file.write(line)
        return id_dict

    if file2:
        ids1 = rename_seqIDs(file1, out_dir, user_format, paired=1)
        ids2 = rename_seqIDs(file2, out_dir, user_format, paired=2)
        with open(out_dir + 'ids2.pkl', 'wb') as pkl_dict:
            pickle.dump(ids2, pkl_dict, protocol=pickle.HIGHEST_PROTOCOL)

        assert len(ids1) == len(ids2), \
            "Paired files have different number of reads"

        for values in range(1,50):
            random_id = random.randint(1, len(ids1)-1)
            id_1 = ids1[random_id].split()[0]
            id_2 = ids2[random_id].split()[0]
            id_overlap = str_overlap(id_1, id_2)
            assert id_overlap == len(id_1) or id_overlap == len(id_1)- 1, \
                "Paired-end sequences don't match"
    else:
        ids1 = rename_seqIDs(file1, out_dir, user_format, paired=False)

    with open(out_dir + 'ids1.pkl', 'wb') as pkl_dict:
        pickle.dump(ids1, pkl_dict, protocol=pickle.HIGHEST_PROTOCOL)


def fastqc_trim(out_dir, file1, trim_minlen, threads, adapter_file, file2=False):
    """Takes fastq data (either single or paired), trims sequences using trimmomatic
    (in the case of paried end reads, it deletes extra files) and uses fastqc to
    show the user what the sequence quality looks like after trimming.

    Returns trimmed sequence files and fastq analysis files
    """
    trimAdapt_command = "ILLUMINACLIP:" + adapter_file + \
                        ":2:30:10 LEADING:20 TRAILING:20 MINLEN:" + \
                        str(trim_minlen)

    if file2:
        subprocess.call("java -jar /mnt/apps/trimmomatic/0.32/trimmomatic.jar PE -threads " + \
                        str(threads) + " " + file1 + " " + file2 + \
                        " trimmed_read1 PE_trimmed_data_1U trimmed_read2 PE_trimmed_data_2U " + \
                        trimAdapt_command, shell=True)
        subprocess.call("rm PE_trimmed_data_1U PE_trimmed_data_2U", shell=True)
        subprocess.call("fastqc trimmed_read1 -o " + out_dir, shell=True)
        subprocess.call("fastqc trimmed_read2 -o " + out_dir, shell=True)


    else:
        subprocess.call("java -jar /mnt/apps/trimmomatic/0.32/trimmomatic.jar SE -threads " + \
                        str(threads) + " " + file1 + " trimmed_read1 " + \
                        trimAdapt_command, shell=True)
        subprocess.call("fastqc trimmed_read1 -o " + out_dir, shell=True)


def kraken_classify(kraken_file1, threads, user_format, kraken_db, kraken_file2=False,
                    quick_minhits=False, preload=False):
    """Uses kraken to classify sequences. Add appropiate switches for kraken
    command (format, preload, minimum hits, if paired or single end) and call
    kraken command, followed by kraken-translate to get full taxonomy for each
    sequence based on thir sequence id (Seq_tax: d__superkingdom, k__kingdom,
    p__phylum, c__class, o__order, f__family, g__genus, s__species).

    Return kraken_table file with a row for each sequence and kraken classification
    (or unclassified) and kraken_labels file witha row for each sequence that was
    classified by kraken with full taxonomy.
    """
    if user_format == "fastq":
        format_switch = " --fastq-input"
    elif user_format == "fasta":
        format_switch = " --fasta-input"

    if preload:
        kraken_command = "kraken --preload "
    else:
        kraken_command = "kraken "

    kraken_command += "--threads " + str(threads) + " --db " + kraken_db + format_switch

    if quick_minhits:
        kraken_command += " --quick --min-hits " + str(quick_minhits)

    if kraken_file2:
        kraken_command += " --paired " +  kraken_file1 + " " + \
                          kraken_file2 + " > kraken_table.txt"
    else:
        kraken_command += " " + kraken_file1 + " > kraken_table.txt"

    subprocess.call(kraken_command, shell=True)
    subprocess.call("kraken-translate --mpa-format --db " + kraken_db + \
                    " kraken_table.txt > kraken_labels.txt", shell=True)

def format_result_table(out_dir, data_table, data_labels, table_colNames, ncbi_file=ncbi_file):
    """Merge the classification data (either kraken or kaiju) with the 'label'
    data which has full taxonomy for the classified sequence. Also adds a column with
    the sequence type (e.g. 'VRL' for virus, 'PLN' for plant)

    Return merged table
    """
    label_colNames=["Seq_ID", "Seq_tax"]
    seq_data = pd.read_csv(out_dir + data_table, sep="\t", header=None, names=table_colNames,
                           index_col=False)
    seq_labelData = pd.read_csv(out_dir + data_labels, sep="\t", header=None,
                                names=label_colNames)
    seq_result = pd.merge(seq_data, seq_labelData, on='Seq_ID', how='outer')

    ncbi_tax = pd.read_csv(ncbi_file, sep=",")
    ncbi_slice = ncbi_tax.iloc[:, [1, 2]]
    seq_final =  pd.merge(seq_result, ncbi_slice, on='Tax_ID', how='outer')
    seq_final = seq_final.dropna(subset=['Seq_ID'])

    return seq_final

def sequence_subset (out_dir, input_file, output_file, user_format, id_list, id_list_name):
    """Create a subset of sequences based on seuqnce IDs in id_list.

    Return subset of sequences file and a text file with list of sequence IDs.
    """
    outfile = open(out_dir  + id_list_name, 'w')
    print >> outfile, "\n".join(str(i) for i in id_list)
    outfile.close()
    id_file = out_dir + id_list_name

    wanted = set(line.rstrip("\n").split(None, 1)[0] for line in open(id_file))
    print "Found %i unique identifiers in %s" % (len(wanted), id_file)

    records = (r for r in SeqIO.parse(input_file, user_format) \
               if r.id in wanted)
    count = SeqIO.write(records, out_dir  + output_file + user_format, user_format)
    print "Saved %i records from %s to %s" % (count, input_file, output_file)
    if count < len(wanted):
        print "Warning %i IDs not found in %s" % (len(wanted)-count, input_file)


def seq_reanalysis(kraken_table, kraken_labels, ncbi_file, out_dir, user_format, forSubset_file1,
                   subset = False, forSubset_file2 = False):
    """Merge kraken_table and kraken_labels using format_result_table() and write to disk
    (delete kraken_table and kraken_label).
    If subset = True, make a list of "Seq_ID" column value if sequence is unclassified
    in "Classified" column or classified as VRL (virus) in column "Div_ID". This list will be
    used to subset sequences using sequence_subset(). This should be used when the host plant
    genome is used to classify sequences.

    Return merged kraken tableresult tables and subsetted sequence files (i subset=True).
    """
    kraken_colNames = ["kraken_classified", "Seq_ID","Tax_ID", "kraken_length", 
                        "kraken_k-mer"]
    kraken_fullTable = format_result_table(out_dir, "kraken_table.txt",
                                           "kraken_labels.txt", kraken_colNames, ncbi_file)
    kraken_fullTable['Seq_ID'] = kraken_fullTable['Seq_ID'].astype(float)
    kraken_fullTable['Seq_ID'] = kraken_fullTable['Seq_ID'].astype(int)

    kraken_results = kraken_fullTable[["kraken_classified", "Seq_ID","Tax_ID", "Seq_tax", 
                                        "Div_ID"]]
    kraken_results.to_csv(out_dir  + 'kraken_VRL.txt', sep='\t', index= False)
    
    with open('ids1.pkl', 'rb') as id_dict:
        ids1 = pickle.load(id_dict)
    kraken_fullTable = kraken_fullTable.replace({"Seq_ID": ids1})
    kraken_fullTable.to_csv(out_dir  + "kraken_FormattedTable.txt", sep='\t', index= False)
    subprocess.call("gzip " + out_dir  + "kraken_FormattedTable.txt", shell =True)
    subprocess.call("rm " + "kraken_table.txt kraken_labels.txt", shell=True)

    if subset:
        unclassified_IDs = kraken_results.loc[(kraken_results.kraken_classified == 'U'), ['Seq_ID']]
        VRL_IDs = kraken_results.loc[(kraken_results.Div_ID == 'VRL'), ['Seq_ID']]
        reanalyse_IDs = unclassified_IDs['Seq_ID'].tolist() + VRL_IDs['Seq_ID'].tolist()

        subset_file1 = forSubset_file1
        delete_file = False
        if forSubset_file2:
            subset_file2 = forSubset_file2
        else:
            subset_file2 = False

        for dirs, sub_dirs, files in os.walk(out_dir):
            if forSubset_file1 in files:
                subset_file1 = out_dir + forSubset_file1
                delete_file = True
                if forSubset_file2:
                    subset_file2 = out_dir + forSubset_file2

        if forSubset_file2:
            reanalyse_ID1 = [str(s) + "/1" for s in reanalyse_IDs]
            reanalyse_ID2 = [str(s) + "/2" for s in reanalyse_IDs]
            sequence_subset(out_dir, subset_file2, "subset_file2.", user_format,
                            reanalyse_ID2, 'reanalyse_ID.txt')
        else:
            reanalyse_ID1 = [str(s) for s in reanalyse_IDs]

        sequence_subset(out_dir, subset_file1, "subset_file1.", user_format,
                        reanalyse_ID1, 'reanalyse_ID.txt')

        if delete_file:
            subprocess.call("rm " + forSubset_file1, shell=True)
            if subset_file2:
                subprocess.call('rm ' + forSubset_file2, shell=True)

        subprocess.call("rm reanalyse_ID.txt", shell=True)


def kaiju_classify(kaiju_file1, threads, out_dir, kaiju_db, kaiju_minlen, kraken_db,
                   kaiju_file2 = False, kaiju_mismatch = False, kaiju_score = False):
    """Run kaiju command for kaiju classification of sequences. It ensures if
    mismatches are allowed that a score has also been provided. Once classification
    is complete, it uses kraken-translate (as in kraken_classify()) to get full
    taxonomy names for each sequence that has been classified. It deletes the files
    used for this analysis.

    """
    kaiju_nodes = kaiju_db + "nodes.dmp"
    kaiju_fmi = kaiju_db + "kaiju_library.fmi"
    # kaiju_names = kaiju_db + "names.dmp"

    if kaiju_mismatch:
        assert(kaiju_score), "Set kaiju_score for greedy mode"
        mode = "greedy -e " + str(kaiju_mismatch) + " -s " + str(kaiju_score)
    else:
        mode = "mem"

    kaiju_command = "kaiju -z " + str(threads) + " -t " + kaiju_nodes + " -f " + kaiju_fmi + \
                    " -i " + kaiju_file1 + " -o kaiju_table.txt -x -v -a " + mode + " -m " + \
                    str(kaiju_minlen)

    if kaiju_file2:
        kaiju_command += " -j " + kaiju_file2

    subprocess.call(kaiju_command, shell=True)
    subprocess.call("kraken-translate --mpa-format --db " + kraken_db + " " +
                    "kaiju_table.txt > kaiju_labels.txt", shell = True)


    for dirs, sub_dirs, files in os.walk(out_dir):
        # Only delete file when it's in out_put
        for filenames in files:
            if kaiju_file1 == filenames:
                subprocess.call("rm " + kaiju_file1, shell=True)
                if kaiju_file2:
                    subprocess.call("rm " + kaiju_file2, shell=True)

def result_analysis(out_dir, kraken_VRL, kaiju_table, kaiju_label, ncbi_file):
    """Imports kraken results table, formats kaiju_table and kaiju_labels and merges
    kraken and kaiju results into one table (kodoja). It then separates reads which
    are classified as VRL, makes a table with all identified viruses and count number
    of intances for each usin virusSummary().
    """
    kraken_results = pd.read_csv(out_dir + kraken_VRL, header=0, sep='\t',
                                 dtype={"kraken_classified": str, "Seq_ID": int,
                                        "Tax_ID": float, "Seq_tax": str, "Div_ID": str})

    kaiju_colNames =["kaiju_classified", "Seq_ID", "Tax_ID", "kaiju_lenBest",
                     "kaiju_tax_AN", "kaiju_accession", "kaiju_fragment"]
    kaiju_fullTable = format_result_table(out_dir, "kaiju_table.txt", "kaiju_labels.txt",
                                          kaiju_colNames, ncbi_file)
    kaiju_fullTable['Seq_ID'] = kaiju_fullTable['Seq_ID'].astype(float)
    kaiju_fullTable['Seq_ID'] = kaiju_fullTable['Seq_ID'].astype(int)
    kaiju_results = kaiju_fullTable[["kaiju_classified", "Seq_ID", "Tax_ID", "Seq_tax", "Div_ID"]]

    with open('ids1.pkl', 'rb') as id_dict:
        ids1 = pickle.load(id_dict)
    kaiju_fullTable = kaiju_fullTable.replace({"Seq_ID": ids1})
    kaiju_fullTable.to_csv(out_dir  + 'kaiju_FormattedTable.txt', sep='\t', index= False)  
    subprocess.call('gzip ' + out_dir  + 'kaiju_FormattedTable.txt', shell =True)

    kodoja = pd.merge(kraken_results, kaiju_results, on='Seq_ID', how='outer')
    assert len(kraken_results) == len(kodoja), \
        'ERROR: Kraken and Kaiju results not merged properly'
    kodoja = kodoja.sort_values(['Seq_ID'])
    kodoja = kodoja.reset_index(drop=True)
    kodoja.rename(columns={'Div_ID_x': 'kraken_div_ID', 'Div_ID_y': 'kaiju_div_ID', 
                           "Seq_tax_x": "kraken_seq_tax", "Seq_tax_y": "kaiju_seq_tax", 
                           'Tax_ID_x': 'kraken_tax_ID', 'Tax_ID_y': 'kaiju_tax_ID'}, inplace=True)

    with open('ids1.pkl', 'rb') as id_dict:
        ids1 = pickle.load(id_dict)
    kodoja = kodoja.replace({"Seq_ID": ids1})

    subprocess.call("rm kaiju_table.txt kaiju_labels.txt kraken_VRL.txt ", shell=True)

    kodoja['combined_result'] = kodoja.kraken_tax_ID[kodoja['kraken_tax_ID'] == kodoja['kaiju_tax_ID']]
    kodoja.to_csv(out_dir  + 'kodoja_VRL.txt', sep='\t', index= False)

    kodoja_vrl = kodoja[(kodoja['kraken_div_ID'] == 'VRL')|(kodoja['kaiju_div_ID'] == 'VRL')]

    def virusSummary(kodoja_data):
        """Creates a summary table with virus species names, tax id, count of
        sequences by kraken, kaiju and sequences that were identified by both
        tools as belonging to that species.

        For each tax id, a sequence count for kraken, kaiju and the combined
        is made. '_levels' dict have all tax ids present in th table with the
        taxanomic 'labels' given by kraken-traslate.

        'associated_tax' dict, has tax ids which would be related to a species
        tax id, as they belong to taxa which are higher, and therefore if
        they could belong to a species but cannot be identified specifically
        (i.e. a sequence whih has been given the following label
        'd__Viruses|f__Closteroviridae|g__Ampelovirus' could be an unspecifically
        identified 'Grapevine_leafroll-associated_virus_4' the label for which is
        'd__Viruses|f__Closteroviridae|g__Ampelovirus|s__Grapevine_leafroll-associated_virus_4').
        """
        kraken_class = dict(kodoja_data['kraken_tax_ID'].value_counts())
        kraken_levels = pd.Series(kodoja_data.kraken_seq_tax.values,
                                  index=kodoja_data.kraken_tax_ID).to_dict()
        kaiju_class = dict(kodoja_data['kaiju_tax_ID'].value_counts())
        kaiju_levels = pd.Series(kodoja_data.kaiju_seq_tax.values,
                                 index=kodoja_data.kaiju_tax_ID).to_dict()
        combined_class = dict(kodoja_data['combined_result'].value_counts())

        levels_dict = kraken_levels.copy()
        levels_dict.update(kaiju_levels)
        print levels_dict
        levels_dict.pop(0, None)
        levels_dict = {k: levels_dict[k] for k in levels_dict if not isnan(k)}
        print levels_dict
        levels_tax = {key: list(map(str, value.split('|')))
                      for key, value in levels_dict.items()}
        print levels_tax
        LCA_tax = {}
        for key, tax in levels_tax.items():
            if tax[-1][0] != 's':
                LCA_tax[key] = tax[-1]
                levels_tax.pop(key, tax)

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

        table_summary = pd.DataFrame(columns=['Species', 'Tax_ID', 'kraken', 'kaiju', 'combined'])
        table_summary['Tax_ID'] =  map(int, levels_tax.keys())
        table_summary['kraken'] = table_summary['Tax_ID'].map(kraken_class)
        table_summary['kaiju'] = table_summary['Tax_ID'].map(kaiju_class)
        table_summary['combined'] = table_summary['Tax_ID'].map(combined_class)
        table_summary['Species'] = table_summary['Tax_ID'].map(species_dict)
        table_summary = table_summary.sort_values('kraken', ascending=False)
        table_summary = table_summary.sort_values('combined', ascending=False)
        table_summary.to_csv(out_dir + 'virus_table.txt', sep='\t', index=False)


    virusSummary(kodoja_vrl)
