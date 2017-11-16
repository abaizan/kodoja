import os
import sys

sys.path.insert(0,'/home/ae42909/viral_diagnostics/parameter_testing/parameterTest_analysis/')
from parameterResults_modules import *

wdir = '/home/ae42909/Scratch/parameter_test/kaiju/run15/'
os.chdir(wdir)

# PB64-S1 (Grape) sRNA-seq
# expected_viruses = {'GRSPaV': 196400, 'GVB':35289, 'GFkV':103722, 'GLRaV-3':55951, 'HSviroid': 12893}

# PB64-S2 (Peach) sRNA-seq
# expected_viruses = {'PNRV':37733}

# PB64-S3 (Raspberry) sRNA-seq
# expected_viruses = {'RBDV':12451, 'RYNV':198310}

# PB64-S4 (Cabbage) sRNA-seq
# expected_viruses = {'PCV':1183241, 'ACV':1476980} # ACV not in refseq?

# PB64-S5 (Sweet potato) sRNA-seq
# expected_viruses = {'SPSMV-1':603333}

# PB64-S7 (Strawberry) RNA-seq
# expected_viruses = {'SMoV': 167161}

# SRR1269627 (Pear) RNA-seq
# expected_viruses = {'ASGV': 28347, 'AGCAV': 1211388, 'ASPV': 35350, 'PrVT': 1472425}

# SRR112893 (Pepper)
expected_viruses = {'ALPV': 209529, 'BPEV': 354328, 'CLCuV':53010 , 'CYVMVA':79236,
                    'PepLCV': 83839,'PepLCBV': 223305, 'PeSV': 157777,'ToLCRnV': 938276,
                     'ToLCBDB': 223344, 'ToLCJoV': 28350, 'GaILV': 1468172,'ToICGV':219299,
                     'TVCV':107324 }
                     # PepLCV and ToLCBDB - betasatellite?

# SRR1291170/71 & SRR1289655, SRR1291169 (Plum)
# expected_viruses = {'PPV': 12211, 'CVA': 42882}


kaiju_colNames =["kaiju_classified", "Seq_ID","Tax_ID", "kaiju_lenBest", "kaiju_tax_AN","kaiju_accession", "kaiju_fragment"]
tool_Results('kaiju', wdir, kaiju_colNames, expected_viruses)
