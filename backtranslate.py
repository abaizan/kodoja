# Python 2.7

# NCBI standard codon table (acquired http://www.bioinformatics.org/JaMBW/2/3/TranslationTables.html)
codon_table = {
    'A': ('GCT', 'GCC', 'GCA', 'GCG'),
    'C': ('TGT', 'TGC'),
    'D': ('GAT', 'GAC'),
    'E': ('GAA', 'GAG'),
    'F': ('TTT', 'TTC'),
    'G': ('GGT', 'GGC', 'GGA', 'GGG'),
    'I': ('ATT', 'ATC', 'ATA'),
    'H': ('CAT', 'CAC'),
    'K': ('AAA', 'AAG'),
    'L': ('TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'),
    'M': ('ATG',),
    'N': ('AAT', 'AAC'),
    'P': ('CCT', 'CCC', 'CCA', 'CCG'),
    'Q': ('CAA', 'CAG'),
    'R': ('CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
    'S': ('TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'),
    'T': ('ACT', 'ACC', 'ACA', 'ACG'),
    'V': ('GTT', 'GTC', 'GTA', 'GTG'),
    'W': ('TGG',),
    'Y': ('TAT', 'TAC'),
    '*': ('TAA', 'TAG', 'TGA'),
}

# Take the first item for each key in the codon_table dictionary to produce a reduced codon library
reduced_codon = {}

for keys in codon_table.keys():
    reduced_codon[keys] = codon_table[keys][0]


# Replace each aa with codon from reduced_codon. Transeq adds letters (eg. 'B' and 'X' for unknown codons) which will be back translated to 'N' as kraken doesn't identify them

def retranslate(protein_file, outfile_name):
    with open(protein_file, 'r') as f:
        protein_list = f.read().splitlines()
        protein_list = filter(None,  protein_list)
    for line in range(0,len(protein_list)):
        if protein_list[line][0] == '>':
            continue
        list_line = list(protein_list[line])
        list_line = [reduced_codon.get(aa, 'N') for aa in list_line]
        protein_list[line] = ''.join(list_line)

    outfile = open(outfile_name, 'w')
    for item in protein_list:
        outfile.write("%s\n" % item)


# Back translate the reference sequence (viral)
#retranslate(protein_file = "/home/ae42909/Scratch/emboss_transeq/twolines_concatinatedRefseqViral_genomic.tax.fna",
#            outfile_name = '/home/ae42909/Scratch/emboss_transeq/retranslated_RefSeqViruses.fna')

# Back translate the fastq files
retranslate(protein_file = "PE1_translated",
            outfile_name = 'PE1_backtranslated')

retranslate(protein_file = "PE2_translated",
            outfile_name = 'PE2_backtranslated')








