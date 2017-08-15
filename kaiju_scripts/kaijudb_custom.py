import pandas as pd

viral_file = '/home/ae42909/Scratch/kaiju/kaijuDB_virus/protein_files/ncbi_downloads/viral.1.protein.faa'

with open(viral_file, 'r') as in_file:
     record = 0
     for line in in_file:
          if line[0] == ">":
               renamed_file.write()
               record = record + 1
          else:
               renamed_file.write(line)
          



# >ref|YP_008320337.1| terminase small subunit [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320338.1| terminase large subunit [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320339.1| portal protein [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320340.1| Clp protease-like protein [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320341.1| major capsid protein [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320342.1| hypothetical protein IBBPl23_06 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320343.1| head-tail connector protein [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320344.1| head-tail joining protein [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320345.1| hypothetical protein IBBPl23_09 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320346.1| hypothetical protein IBBPl23_10 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320347.1| major tail protein [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320348.1| hypothetical protein IBBPl23_12 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320349.1| hypothetical protein IBBPl23_13 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320350.1| tail length tape-measure protein [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320351.1| tail protein [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320352.1| hypothetical protein IBBPl23_16 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320353.1| tail protein [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320354.1| hypothetical protein IBBPl23_18 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320355.1| hypothetical protein IBBPl23_19 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320356.1| putative membrane protein [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320357.1| N-acetylmuramoyl-L-alanine amidase [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320358.1| hypothetical protein IBBPl23_21A [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320359.1| putative membrane protein [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320360.1| hypothetical protein IBBPl23_23 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320361.1| hypothetical protein IBBPl23_24 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320362.1| hypothetical protein IBBPl23_25 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320363.1| toxin 1 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320364.2| hypothetical protein IBBPl23_27 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320365.1| hypothetical protein IBBPl23_27A [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320366.1| hypothetical protein IBBPl23_27B [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320367.1| hypothetical protein IBBPl23_28 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320368.1| hypothetical protein IBBPl23_29 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320369.1| recombination protein [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320370.1| SOS-response repressor and protease LexA [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320371.1| hypothetical protein IBBPl23_32 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320372.1| hypothetical protein IBBPl23_33 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320373.1| hypothetical protein IBBPl23_34 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320374.1| hypothetical protein IBBPl23_35 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320375.1| hypothetical protein IBBPl23_36 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320376.1| hypothetical protein IBBPl23_37 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320377.1| hypothetical protein IBBPl23_38 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320378.1| antirepressor protein [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320379.1| hypothetical protein IBBPl23_40 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320380.1| hypothetical protein IBBPl23_41 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320381.1| hypothetical protein IBBPl23_42 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320382.1| hypothetical protein IBBPl23_43 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320383.1| hypothetical protein IBBPl23_44 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320384.1| hypothetical protein IBBPl23_45 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320385.1| recombinational DNA repair protein RecT [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320386.1| beta-lactamase [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320387.1| hypothetical protein IBBPl23_48 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320388.1| hypothetical protein IBBPl23_49 [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320389.1| replicative DNA helicase [Paenibacillus phage phiIBB_Pl23]

# >ref|YP_008320390.1| hypothetical protein IBBPl23_51 [Paenibacillus phage phiIBB_Pl23]


          


