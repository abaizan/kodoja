fastq_file1 = "/home/ae42909/Scratch/100_Potato_withViruses_1.fastq"
fastq_file2 = "/home/ae42909/Scratch/100_Potato_withViruses_2.fastq"
fastaFile_name1 = "/home/ae42909/Scratch/100_Potato_withViruses_1.fasta"
fastaFile_name2 = "/home/ae42909/Scratch/100_Potato_withViruses_2.fasta"

def fastq_to_fasta(fastq_file, fastaFile_name):
    fasta_file = open(fastaFile_name, "w")
    with open(fastq_file, 'r') as in_file:
        seq_num = 0
        for lineNum, line in enumerate(in_file):
                if lineNum % 4 == 0:
                    list_line = list(line)
                    seq_num += 1
                    fasta_file.write(">" + "".join(list_line[1:]))
                if (lineNum - 1) % 4 == 0:
                    fasta_file.write(line)
    return seq_num


fastq1 = fastq_to_fasta(fastq_file1, fastaFile_name1)
fastq2 = fastq_to_fasta(fastq_file2, fastaFile_name2)
 


