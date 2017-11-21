import os
import time
import subprocess
import uuid
import sys
sys.path.insert(0, '/home/ae42909/viral_diagnostics/parameter_testing/')
from kaiju_parameterValues import kaiju_parameters
sys.path.insert(0, '/home/ae42909/viral_diagnostics/diagnosticTool_scripts/')
from diagnostic_modules import *

# General parameters
# CHECK KAIJU RUN NUMBERS! ll ~/Scratch/parameter_test/kaiju/
run_list = ['run28']
fileDir = '/home/ae42909/data_forTesting/GrapevineJo_data/set4/'
data_type = 'PE'

# Static parameters
user_format = "fastq"
out_dir = "/home/ae42909/Scratch/parameter_test/kaiju/"
threads = 4
kraken_db = "/home/ae42909/Scratch/krakenDB_viral"
kaiju_db = "/home/ae42909/Scratch/kaijuDB_viral/"

def kaiju_param():    
    t0 = time.time()
    kaiju_classify(file1, threads, out_dir, kaiju_db, kaiju_minlen, kraken_db, file2,
                   kaiju_mismatch = kaiju_mismatch, kaiju_score = kaiju_score)
    subprocess.call("mv " + run_dir + "kaiju_table.txt " + run_dir + "kaiju_table_" + tag + ".txt", shell = True)
    subprocess.call("mv " + run_dir + "kaiju_labels.txt " + run_dir + "kaiju_labels_" + tag + ".txt", shell = True)
    t1 = time.time()
    with open(run_dir + 'log_file_' + tag + ".txt", 'w') as out_file:
        text = "input1 = " + str(file1) + "\n" + "input2 = " + str(file2) + "\n"
        text += "kaiju database = " + kaiju_db + "\n" + "threads = " + str(threads) + "\n"
        text += "kaiju_minlen = " + str(kaiju_minlen) + "\n" + "kaiju_mismatch =" + str(kaiju_mismatch) + "\n" + "kaiju_score =" + str(kaiju_score) + "\n"
        text += "time (min) = " + str((t1-t0)/60) + "\n"
        out_file.write(text)


for dirs, sub_dirs, files in os.walk(fileDir):
    if data_type == 'PE':
        assert((len(files)/2) == len(run_list)), \
               "Don't have enough runs for number of files"
        files = sorted(files)
        for run_idx, runs in enumerate(run_list):
            run_dir = out_dir + runs + "/"
            file1 = fileDir + files[run_idx*2]
            file2 = fileDir + files[run_idx*2+1]

            if os.path.exists(run_dir):
                sys.exit("run number already exists")
            else:    
                os.makedirs(run_dir)
            os.chdir(run_dir)

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
            log_file = open(out_dir + runs + "_file.txt", "w")
            log_file.write("input1 = " + str(file1) + "\n" + "input2 = " + str(file2) + "\n"
                           "kaiju_parameters = " + str(kaiju_parameters) + "\n")
            log_file.close()
    else:
        assert(len(files) == len(run_list)), \
               "Don't have enough runs for number of files"
        for run_idx, runs in enumerate(run_list):
            run_dir = out_dir + runs + "/"
            file1 = fileDir + files[run_idx]
            file2 = False

            if os.path.exists(run_dir):
                sys.exit("run number already exists")
            else:    
                os.makedirs(run_dir)
            os.chdir(run_dir)

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
            log_file = open(out_dir + runs + "_file.txt", "w")
            log_file.write("input1 = " + str(file1) + "\n" + "input2 = " + str(file2) + "\n" +
                           "kaiju_parameters = " + str(kaiju_parameters) + "\n")
            log_file.close()





        





