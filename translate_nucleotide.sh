#!/bin/bash

sed 's/\/1$//' PE1_subset > PE1_short # remove the '/1' to be added later for kraken analysis, as kraken recognises "/1" as a paired read
sed 's/\/2$//' PE2_subset > PE2_short

transeq PE1_short PE1_trans -frame=6
transeq PE2_short PE2_trans -frame=6


cat PE1_trans |  awk '{print (NR%2 == 1) ? $0 "/1"  :$0}' > PE1_translated
cat PE2_trans |  awk '{print (NR%2 == 1) ? $0 "/2"  :$0}' > PE2_translated
