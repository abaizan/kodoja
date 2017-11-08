# NOTE: diagnostic modules needs a defult ncbi_file for this script to work (ncbi_file defined within modules)
import os
import sys

sys.path.insert(0,'/home/ae42909/viral_diagnostics/parameter_testing/parameterTest_analysis/')
from parameterResults_modules import *

wdir = '/home/ae42909/Scratch/parameter_test/kraken/run12/'
os.chdir(wdir)

# PB64-S1 (Grape)
# expected_viruses = {'GRSPaV': 196400, 'GVB':35289, 'GFkV':103722, 'GLRaV-3':55951, 'HSviroid': 12893}

# PB64-S2 (Peach)
# expected_viruses = {'PNRV':37733}

# PB64-S3 (Raspberry)
# expected_viruses = {'RBDV':12451, 'RYNV':198310}

# PB64-S4 (Cabbage)
# expected_viruses = {'PCV':1183241, 'ACV':1476980} # ACV not in refseq?

# PB64-S5 (Sweet potato)
# expected_viruses = {'SPSMV-1':603333}

# PB64-S7 (Strawberry)
# expected_viruses = {'SMoV': 167161}

# SRR1269627 (Pear)
# expected_viruses = {'ASGV': 28347, 'AGCAV': 1211388, 'ASPV': 35350, 'PrVT': 1472425}

# SRR112893 (Pepper)
expected_viruses = {'ALPV': 209529, 'BPEV': 354328, 'CLCuV':53010 , 'CYVMVA':79236,
                    'PepLCV': 83839,'PepLCBV': 223305, 'PeSV': 157777,
                    'ToLCRnV': 938276,'ToLCBDB': 223344, 'ToLCJoV': 28350, 'GaILV': 1468172,
                    'ToICGV':219299,'TVCV':107324 }
                     # PepLCV and ToLCBDB - betasatellite?


kraken_colNames = ["kraken_classified", "Seq_ID","Tax_ID", "kraken_length", "kraken_k-mer"]
tool_Results('kraken', wdir, kraken_colNames, expected_viruses)
