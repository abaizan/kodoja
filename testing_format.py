
import os

user_format = os.environ["FORMAT"]
input_file  = os.environ["FILE1"]

with open(input_file) as myfile:
    small_file = [next(myfile) for x in xrange(8)]

file_format = 'empty'

if small_file[0][0] == '@' and small_file[3][0] == '@':
    file_format = 'fastq'
if small_file[0][0] == '>' and small_file[2][0] == '>':
    file_format = 'fasta'


assert (file_format == "fasta") | (file_format == "fastq"), "Cannot proceed with file as it is not in fasta or fastq format."

if not user_format == file_format:
    print 'File has been detected to be in ' + file_format + ' format rather than ' + user_format + 'format. Continuing as ' + file_format


