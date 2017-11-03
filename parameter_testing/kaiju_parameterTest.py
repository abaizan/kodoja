import os
import time
import subprocess
import uuid
import sys
from kaiju_parameterValues import kaiju_parameters

sys.path.insert(0, '/home/ae42909/viral_diagnostics/diagnosticTool_scripts/')
from diagnostic_modules import *


# General parameters
file1 = '/mnt/shared/projects/virology/201609_BBSRC_Diagnostics/Data/synthetic/review-paper-test-datasets/SRR1269627Pear.fastq'
file2 = False
out_dir = "/home/ae42909/Scratch/parameter_test/kaiju/"
threads = 4
kraken_db = "/home/ae42909/Scratch/krakenDB_viral"
kaiju_db = "/home/ae42909/Scratch/kaijuDB_viral/"

run_dir = out_dir + "run10" + "/"
if os.path.exists(run_dir):
    sys.exit("run number already exists")
else:    
    os.makedirs(run_dir)

os.chdir(run_dir)

def kaiju_param():    
    t0 = time.time()
    if file2:
        kaiju_classify(file1, threads, kaiju_db, kaiju_minlen, kraken_db, file2, kaiju_mismatch, kaiju_score)
        text = "input1 = " + file1 + user_format + "\n" + "input2 = " + file2 + user_format + "\n"
    else:
        kaiju_classify(file1, threads, kaiju_db, kaiju_minlen, kraken_db, kaiju_file2 = False, kaiju_mismatch = kaiju_mismatch, kaiju_score = kaiju_score)
        text = "input1 = " + file1 + "\n"


    subprocess.call("mv " + run_dir + "kaiju_table.txt " + run_dir + "kaiju_table_" + tag + ".txt", shell = True)
    subprocess.call("mv " + run_dir + "kaiju_labels.txt " + run_dir + "kaiju_labels_" + tag + ".txt", shell = True)
    t1 = time.time()
    with open(run_dir + 'log_file_' + tag + ".txt", 'w') as out_file:
        text += "kaiju database = " + kaiju_db + "\n" + "threads = " + str(threads) + "\n"
        text += "kaiju_minlen = " + str(kaiju_minlen) + "\n" + "kaiju_mismatch =" + str(kaiju_mismatch) + "\n" + "kaiju_score =" + str(kaiju_score) + "\n"
        text += "time (min) = " + str((t1-t0)/60) + "\n"
        out_file.write(text)

#    subprocess.call("mv *.txt " + run_dir, shell = True)



for minlens in kaiju_parameters["kaiju_minlen"]:
    for mismatches in kaiju_parameters["kaiju_mismatch"]:
        kaiju_minlen = minlens
        kaiju_mismatch = mismatches
        if kaiju_mismatch:
            for scores in kaiju_parameters["kaiju_score"]:
                kaiju_score = scores
                tag = str(uuid.uuid4())
                kaiju_param()      
        else:
            kaiju_score = False
            tag = str(uuid.uuid4())
            kaiju_param()
        





