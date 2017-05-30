#!/bin/bash

if [[ "$FORMAT" =~ ^(fastq|fasta)$ ]]; then
    echo "$FORMAT is in the list"
else
    exit "Incorrect format"
fi
