
import os
import time
import subprocess
import uuid
import sys
from diagnostic_modules import *
from kraken_parameterValues import kraken_parameters

# General parameters
file1 = "/home/ae42909/synthetic_RNAseq/mappingRNAseq/concatenated_fastaFiles/Potato_withViruses_1."
file2 = "/home/ae42909/synthetic_RNAseq/mappingRNAseq/concatenated_fastaFiles/Potato_withViruses_2."
# file1 = "/home/ae42909/Scratch/100_Potato_withViruses_1."
# file2 = "/home/ae42909/Scratch/100_Potato_withViruses_2."
user_format = "fastq"
out_dir = "/home/ae42909/Scratch/parameter_test/kraken/"
threads = 4
user_format = "fastq"
kraken_db_dir = "/home/ae42909/Scratch/parameter_test/kraken/"

run = "run2"
if os.path.exists(out_dir + run):
    sys.exit("run number already exists")
else:    
    os.makedirs(out_dir + run)

os.chdir(out_dir)

for root, subdirs, files in os.walk(kraken_db_dir):
    for dirs in subdirs:
        if dirs[:6] == "kraken":
            for minhits in kraken_parameters["quick_minhits"]:
                for preload_switch in kraken_parameters["preload"]:
                    kraken_db = os.path.join(root, dirs)
                    quick_minhits = minhits
                    preload = preload_switch
                    tag = str(uuid.uuid4())

                    t0 = time.time()
                    if file2:
                        # Kraken classification
                        kraken_classify(file1, threads, user_format, kraken_db, file2, quick_minhits, preload)
                        text = "input1 = " + file1 + user_format + "\n" + "input2 = " + file2 + user_format + "\n"
                    else:
                        # Kraken classification
                        kraken_classify(file1, threads, user_format, kraken_db, quick_minhits, preload)
                        text = "input1 = " + file1 + "\n"

                        
                    subprocess.call("mv " + out_dir + "kraken_table.txt " + out_dir + "kraken_table_" + tag + ".txt", shell = True)
                    subprocess.call("mv " + out_dir + "kraken_labels.txt " + out_dir + "kraken_labels_" + tag + ".txt", shell = True)
                    t1 = time.time()
                    with open(out_dir + tag + ".txt", 'w') as out_file:
                        text += "kraken database = " + kraken_db + "\n" + "threads = " + str(threads) + "\n"
                        text += "quick_minhits = " + str(minhits) + "\n" + "preload =" + str(preload) + "\n"
                        text += "time (min) = " + str((t1-t0)/60) + "\n"
                        out_file.write(text)


subprocess.call("mv *.txt " + out_dir + run, shell = True)


