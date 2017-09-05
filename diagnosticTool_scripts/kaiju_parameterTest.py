import os
import time
import subprocess
import uuid
import sys
from diagnostic_modules import *
from kaiju_parameterValues import kaiju_parameters


# General parameters
file1 = "/home/ae42909/synthetic_RNAseq/mappingRNAseq/concatenated_fastaFiles/Potato_withViruses_1."
file2 = "/home/ae42909/synthetic_RNAseq/mappingRNAseq/concatenated_fastaFiles/Potato_withViruses_2."
# file1 = "/home/ae42909/Scratch/100_Potato_withViruses_1."
# file2 = "/home/ae42909/Scratch/100_Potato_withViruses_2."
user_format = "fastq"
out_dir = "/home/ae42909/Scratch/parameter_test/kaiju/"
threads = 4
kraken_db = "/home/ae42909/Scratch/krakenDB_viral"
kaiju_db = "/home/ae42909/Scratch/kaijuDB_viral/"

run_dir = out_dir + "run2"
if os.path.exists(run_dir):
    sys.exit("run number already exists")
else:    
    os.makedirs(run_dir)

os.chdir(out_dir)

def kaiju_param():    
    t0 = time.time()
    if file2:
        kaiju_classify(file1, user_format, threads, kaiju_db, kaiju_minlen, kraken_db, file2, kaiju_mismatch, kaiju_score)
        text = "input1 = " + file1 + user_format + "\n" + "input2 = " + file2 + user_format + "\n"
    else:
        kaiju_classify(file1, user_format, threads, kaiju_db, kaiju_minlen, kraken_db, subset_file2, kaiju_mismatch, kaiju_score)
        text = "input1 = " + file1 + "\n"


    subprocess.call("mv " + out_dir + "kaiju_table.txt " + out_dir + "kaiju_table_" + tag + ".txt", shell = True)
    subprocess.call("mv " + out_dir + "kaiju_labels.txt " + out_dir + "kaiju_labels_" + tag + ".txt", shell = True)
    t1 = time.time()
    with open(out_dir + tag + ".txt", 'w') as out_file:
        text += "kaiju database = " + kaiju_db + "\n" + "threads = " + str(threads) + "\n"
        text += "kaiju_minlen = " + str(kaiju_minlen) + "\n" + "kaiju_mismatch =" + str(kaiju_mismatch) + "\n" + "kaiju_score =" + str(kaiju_score) + "\n"
        text += "time (min) = " + str((t1-t0)/60) + "\n"
        out_file.write(text)

    subprocess.call("mv *.txt " + run_dir, shell = True)



for minlens in kaiju_parameters["kaiju_minlen"]:
    for mismatches in kaiju_parameters["kaiju_mismatch"]:
        tag = str(uuid.uuid4())
        kaiju_minlen = minlens
        kaiju_mismatch = mismatches
        if kaiju_mismatch:
            for scores in kaiju_parameters["kaiju_score"]:
                kaiju_score = scores
                kaiju_param()      
        else:
            kaiju_score = False
            kaiju_param()






