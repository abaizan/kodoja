#!/bin/bash

#$ -cwd
#$ -j yes
THREADS=4
TRIMLEN=$1
FILE1=$2

#trimmomatic requires both PE files to be called the same so file2 can b recognised by file1

fastqc $FILE1 -o $OUTDIR

# If $FILE2 is empty (i.e. single-end data) the do trimmomatic as SE, otherwise do paired-end
if [ -z ${FILE2+x} ]
then 
    trimmomatic SE -threads $THREADS $FILE1 SE_trimmed_data ILLUMINACLIP:$ADAPTERS:2:30:10 LEADING:20 TRAILING:20 MINLEN:$TRIMLEN
     
else
    FILE2=$3
    fastqc $FILE2 -o $OUTDIR
    trimmomatic PE -basein $FILE1 -baseout PE_trimmed_data -threads $THREADS  ILLUMINACLIP:$ADAPTERS:2:30:10 LEADING:20 TRAILING:20 MINLEN:$TRIMLEN
fi
