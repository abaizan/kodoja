# python 2.7.13
import subprocess
import pandas as pd
import numpy as np
import os
from Bio import SeqIO
from itertools import izip

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
def fastqc_trim(out_dir, file1, trimlen, threads, file2 = False):

    subprocess.call("fastqc " + file1 + " -o " + out_dir, shell=True)

    if file2:
        subprocess.call("fastqc " + file2 + " -o " + out_dir, shell=True)
        subprocess.call("java -jar /mnt/apps/trimmomatic/0.32/trimmomatic.jar PE -threads " + threads + " " + file1 + " " + file2 + 
                        " PE_trimmed_data_1P PE_trimmed_data_1U PE_trimmed_data_2P PE_trimmed_data_2U " +
                        "ILLUMINACLIP:$ADAPTERS:2:30:10 LEADING:20 TRAILING:20 MINLEN:" + trimlen, shell=True)
        subprocess.call("rm PE_trimmed_data_1U PE_trimmed_data_2U", shell=True)
    else:
        subprocess.call("java -jar /mnt/apps/trimmomatic/0.32/trimmomatic.jar SE -threads " + threads + " " + file1 +
                        " SE_trimmed_data " +
                        "ILLUMINACLIP:$ADAPTERS:2:30:10 LEADING:20 TRAILING:20 MINLEN:" + trimlen, shell=True)


# Order and replace sequence IDs with numberic IDs
def rename_seq(trim_file, out_dir, user_format, paired=False):
    if user_format == "fasta":
        form_line = 2
        symb = ">"
    if user_format == "fastq":
        form_line = 4
        symb = "@"

    ids = sorted(rec.id for rec in SeqIO.parse(trim_file, user_format))
    record_index = SeqIO.index(trim_file, user_format)
    records = (record_index[id] for id in ids)
    SeqIO.write(records, out_dir + "sorted", user_format)
    subprocess.call("rm " + trim_file, shell=True)

    def replace_ids(sorted_file, out_dir, user_format, pairedNum=False):
        if not pairedNum:
            renamed_file = open(out_dir + "renamed_file." + user_format, "w")
        else:
            renamed_file = open(out_dir + "renamed_file" + str(pairedNum) + "." + user_format, "w")

        with open(sorted_file, 'r') as in_file:
            seqNum = 1
            for lineNum, line in enumerate(in_file):
                if lineNum % form_line == 0:
                    renamed_file.write(symb + str(seqNum) + "/" + str(pairedNum) + "\n")
                    seqNum = seqNum + 1
                else:
                    renamed_file.write(line)


    replace_ids(out_dir + "sorted", out_dir, user_format, pairedNum=paired)
    subprocess.call("rm " + out_dir + "sorted", shell=True)


# IDs = [] 
# with open(trim_file, 'r') as in_file:
#     for lineNum, line in enumerate(in_file):
#         if lineNum % form_line == 0:
#             IDs.append(all_lines)
       
# np.argsort(IDs)
        
#     for line in in_file:
#         if lineNum % form_line == 0:
#             lines_toJoin = [ for line in in_file][4:8]
#             lineNum = lineNum+1
#             print lines

    
#     for line in f:
#         if lineNum % form_line == 0:
#             new_line = "".join(f[lineNum:(lineNum+form_line)]).replace("\n","\t", form_line-1) # One line per sequence (only leave newline (\n) at end)
#             seqLine_temp.write(new_line)
#         lineNum = lineNum + 1

# seqLine_temp.flush()
# with open(seqLine_temp.name) as temp_file:
#     sorted_file = sorted(temp_file.readlines())
#     id_values = []
#     line = 0
#     for seqs in sorted_file:
#         id_values.append(seqs.split("\t", 1)[0])




   
#     for line in xrange(len(sorted_file)):
#         if line % form_line == 0:
#             id_values.append(sorted_file[line].split("\n", 1)[0])

# id_dict = dict(list(enumerate(id_values)))




# take the first string (sep \t) and put into dictionary
# if paired = 1: change first string (sep \t) to numbers + "/1", if paried = 2: change first string (sep \t) to numbers + "/2"; else: change first string (sep \t) to numbers

# with open(trim_file1) as f:
#     trim_file = [next(f) for x in xrange(40)]
	    
# trim_file = "/home/ae42909/Scratch/test/PE_test"  



# Kraken classification
# def kraken_class()

# Subset viral and unclassified sequences
# seq_reanalysis()

# Kaiju classification of subset sequences
# kaiju_class()

# Merege results
# result_analysis()
