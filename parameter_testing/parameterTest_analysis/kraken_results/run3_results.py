# run on  ipython
import os

os.chdir('/home/ae42909/viral_diagnostics/parameter_testing/parameterTest_analysis/')
from parameterResults_modules import *

wdir = '/home/ae42909/Scratch/parameter_test/kraken/run3/'
os.chdir(wdir)

# PB64-S1
expected_viruses = {'GRSPaV': 196400, 'GVB':35289, 'GFkV':103722, 'GLRaV-3':55951, 'HSviroid': 12893}

# PB64-S2
# expected_viruses = {'PNRV':37733}

# PB64-S3
# expected_viruses = {'RBDV':12451, 'RYNV':198310}

# PB64-S5
# expected_viruses = {'SPSMV-1':603333}

# PB64-S7
# expected_viruses = {'SMoV': 167161}


kraken_colNames = ["kraken_classified", "Seq_ID","Tax_ID", "kraken_length", "kraken_k-mer"]
# kmerResults(wdir, kraken_colNames, expected_viruses)
tool_Results('kraken', wdir, kraken_colNames, expected_viruses)
