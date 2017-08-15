#!/bin/bash

#$ -cwd
#$ -j yes

# For fastq data only - need to add for fasta format

if [ -z ${FILE2+x} ]
then 
    cat SE_trimmed_data | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > SE_sorted #order
    awk 'NR == 1 || NR % 4 == 1' SE_sorted > ID
    cat SE_sorted | awk '{print (NR%4 == 1) ? "@" ++i "" : $0 }'  > SE_renamed # rename to N
    
else
    cat PE_trimmed_data_1P | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > PE1_sorted #order
    cat PE_trimmed_data_2P | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > PE2_sorted

    awk 'NR == 1 || NR % 4 == 1' PE1_sorted > ID1
    awk 'NR == 1 || NR % 4 == 1' PE2_sorted > ID2

    cat PE1_sorted | awk '{print (NR%4 == 1) ? "@" ++i "/1" : $0 }'  > PE1_renamed # rename to N/1
    cat PE2_sorted | awk '{print (NR%4 == 1) ? "@" ++i "/2" : $0 }'  > PE2_renamed # rename to N/1
   
fi


