def check_file(file1, out_dir, user_format, file2=False):
    """Rename sequnce ids for SE or PE files to ensure
    consistency between kraken and kaiju (which modify
    id names). Create dictionaries containing real IDs and 
    renamed version and pickle. If data is PE, assert 
    paired files have the same number of entries and if 
    the paired reads are mached by choosing a random 
    entry from the first list of ids (minus meatdata) and 
    the same entry for the second list of ids (can have 
    one character different  as could be named 1 or 2).
    """
    def str_overlap(str1,str2):
        """ Determine number of matching characters between
        two strings.

        Returns int with number of matching characters
        """

        count = 0
        for i in range(min(len(str1), len(str2))):
            if str1[i] == str2[i]:
                count = count + 1
        return count

    def rename_seqIDs(input_file, out_dir, user_format, paired=False):
        """ Write a new file where each sequence ID is replaced with 
        the the first character (">" or "@") followed by a number 
        (1::N), and if it's a paired read followd by "/1" or "/2"
        (called 'renamed_file'), and a dictionary composed of 
        dict[newID] = oldID.

        """
        output_file = out_dir + "renamed_file_"
        if paired == 2:
            output_file += str(paired) + "." + user_format
        else:
            output_file += "1." + user_format
        id_dict = {}
        with open(input_file, 'r') as in_file, open(output_file, 'w') as out_file:
            seqNum = 1
            if user_format == 'fastq':
                seqid_line = 4
            else:
                seqid_line = 2

            for lineNum, line in enumerate(in_file):
                if lineNum %  seqid_line == 0:
                    new_line = line[0] + str(seqNum)
                    id_dict[seqNum] = line.strip()
                    if paired:
                        new_line += '/' + str(paired)
                    new_line += "\n"
                    out_file.write(new_line)
                    seqNum += 1
                else:
                    out_file.write(line)
        return id_dict

    if file2:
        ids1 = rename_seqIDs(file1, out_dir, user_format, paired=1)
        ids2 = rename_seqIDs(file2, out_dir, user_format, paired=2)
        with open(out_dir + 'ids2.pkl', 'wb') as pkl_dict:
            pickle.dump(ids2, pkl_dict, protocol=pickle.HIGHEST_PROTOCOL)

        assert len(ids1) == len(ids2), \
            "Paired files have different number of reads"

        for values in range(1,50):
            random_id = random.randint(1, len(ids1)-1)
            id_1 = ids1[random_id].split()[0]
            id_2 = ids2[random_id].split()[0]
            id_overlap = str_overlap(id_1, id_2)
            assert  id_overlap == len(id_1) or id_overlap == len(id_1)- 1, \
                "Paired-end sequences don't match"
    else:
        ids1 = rename_seqIDs(file1, out_dir, user_format, paired=False)

    with open(out_dir + 'ids1.pkl', 'wb') as pkl_dict:
        pickle.dump(ids1, pkl_dict, protocol=pickle.HIGHEST_PROTOCOL)




