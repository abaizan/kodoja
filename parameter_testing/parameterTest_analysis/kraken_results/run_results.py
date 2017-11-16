# NOTE: diagnostic modules needs a defult ncbi_file for this script to work (ncbi_file defined within modules)
import os
import sys

sys.path.insert(0,'/home/ae42909/viral_diagnostics/parameter_testing/parameterTest_analysis/')
from parameterResults_modules import *

wdir = '/home/ae42909/Scratch/parameter_test/kraken/run15/'
os.chdir(wdir)

expected_viruses = {'PPV': 12211, 'CVA': 42882}


kraken_colNames = ["kraken_classified", "Seq_ID","Tax_ID", "kraken_length", "kraken_k-mer"]
tool_Results('kraken', wdir, kraken_colNames, expected_viruses)
