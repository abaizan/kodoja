import pandas as pd

out_dir = "/home/ae42909/Scratch/diagnostic_results/Barrero_data_results/PB64-S1_2/"
file1 = "/home/ae42909/data_forTesting/Barerro_data/PB64-S1_clean.fq.gz_trim.fastq"
user_format = "fastq"
userList_TaxID = [71186]

table_summary = pd.read_csv(out_dir + "virus_table.txt", header = 0, sep='\t')
kodoja_vrl = pd.read_csv(out_dir + "", header = 0, sep='\t')

if userList_TaxID:
    TaxId_out = userList_TaxID
    
else:
    TaxId_out = list(table_summary.Tax_ID)


rows_wanted = (kodoja_vrl['kraken_tax_ID'].isin(TaxId_out) | kodoja_vrl['kaiju_tax_ID'].isin(TaxId_out))
seqID_wanted = list(kodoja_vrl.loc[rows_wanted, 'Seq_ID'])

def sequence_subset (out_dir, input_file, output_file, user_format, id_list, id_list_name):
    """Create a subset of sequences based on seuqnce IDs in id_list.

    Return subset of sequences file and a text file with list of sequence IDs.
    """
    outfile = open(out_dir  + id_list_name, 'w')
    print >> outfile, "\n".join(str(i) for i in id_list)
    outfile.close()
    id_file = out_dir + id_list_name

    wanted = set(line.rstrip("\n").split(None,1)[0] for line in open(id_file))
    print "Found %i unique identifiers in %s" % (len(wanted), id_file)

    records = (r for r in SeqIO.parse(input_file, user_format) \
               if r.id in wanted)
    count = SeqIO.write(records, out_dir  + output_file + user_format, user_format)
    print "Saved %i records from %s to %s" % (count, input_file, output_file)
    if count < len(wanted):
        print "Warning %i IDs not found in %s" % (len(wanted)-count, input_file)

sequence_subset(out_dir, file1, "viral_sequences.", user_format, seqID_wanted, 'viral_sequences_list.txt')


 
        # Use biopython to make new fastq files of sequences to be reanalysed
        # def reanalyse_subset (input_file, output_file, user_format, id_list):
        #     """Create a subset of sequences based on seuqnce IDs in id_list.

        #     Return subset of sequences.
        #     """
        #     outfile = open(out_dir  + 'reanalyse_ID.txt', 'w')
        #     print >> outfile, "\n".join(str(i) for i in id_list)
        #     outfile.close()
        #     id_file = out_dir  + "reanalyse_ID.txt"

        #     wanted = set(line.rstrip("\n").split(None,1)[0] for line in open(id_file))
        #     print "Found %i unique identifiers in %s" % (len(wanted), id_file)

        #     records = (r for r in SeqIO.parse(out_dir + input_file, user_format) \
        #                if r.id in wanted)
        #     count = SeqIO.write(records, out_dir  + output_file + user_format, user_format)
        #     print "Saved %i records from %s to %s" % (count, input_file, output_file)
        #     if count < len(wanted):
        #         print "Warning %i IDs not found in %s" % (len(wanted)-count, input_file)
