import os
import sys

sys.path.insert(0,'/home/ae42909/viral_diagnostics/parameter_testing/parameterTest_analysis/')
from parameterResults_modules import *

initial_dir = '/home/ae42909/Scratch/parameter_test/kaiju/'
kaiju_colNames =["kaiju_classified", "Seq_ID","Tax_ID", "kaiju_lenBest", "kaiju_tax_AN","kaiju_accession", "kaiju_fragment"]

expected_virus = {'PVA':, 12215, 'PVS':12169, 'PVX':12183, 'PLRV':12045, 'PVT':36403, 'APMMV':1296569}
runNums = [47, 48, 49, 52, 54]

for run in runNums:
    wdir = initial_dir + 'run' + str(run) + '/'
    os.chdir(wdir)
    tool_Results('kaiju', wdir, kaiju_colNames, expected_viruses)
