import os
import sys
sys.path.insert(0, '/home/ae42909/viral_diagnostics/diagnosticTool_scripts/')
from diagnostic_modules import *

os.chdir(fileDir)

# types of data to test: '@seqid 1','@seqid 1:more', '@seqid/1', ',1', '|1', '\1', '^1', '!1'

if user_format == 'fastq':
    seqid_line = 4
else:
    seqid = 2

input_file = '20_Potato_withViruses_1.fastq'
output_file = '20.Nspace1_Potato_withViruses_1.fastq'
seqid_task = line.split('.')[0] + '.' + str(count) + ' 1\n'
count = 1

with open(output_file, 'w') as out_file, open(input_file, 'r') as in_file:
    for lineNum, line in enumerate(in_file):
        if lineNum % seqid_line == 0:
            print seqid_task
            out_file.write(seqid_task)
            count += 1
        else:
            out_file.write(line)


def paired_ids(fname, user_format):
    """Get list of sequence identifires for each paired file.

    Return list_ids and write list_ids to a text file in working
    directory (each id on a new line).
    """
    list_ids = []
    format_num = 4
    if user_format == "fasta":
        format_num = 2
    with open(fname, 'r') as in_file:
        for lineNum, line in enumerate(in_file):
            if lineNum % format_num == 0:
                list_ids.append(line)

    return list_ids

ids1 = paired_ids(file1, user_format)
ids2 = paired_ids(file2, user_format)



fileDir = '/home/ae42909/Scratch/smallTest_data/'
file1 = fileDir + '20_Potato_withViruses_1.fastq'
file2 = False
user_format = "fastq"
out_dir = "/home/ae42909/Scratch/parameter_test/kraken/krakenDB_k31_m15/"
threads = 1
user_format = "fastq"
kraken_db = "/home/ae42909/Scratch/parameter_test/kraken/databases/"

kraken_classify(file1, threads, user_format, kraken_db, quick_minhits = quick_minhits, preload = preload)
