#!/bin/bash

#$ -cwd
#$ -j yes
#$ -pe smp 4
THREADS=4
TREATMENT=bam_data/RAGko_2_sorted.bam
CONTROL=bam_data/RAGko_1_sorted.bam
TREAT_OUT=bootstrap/RAG2_FP
CON_OUT=bootstrap/RAG1_FP

wellington_bootstrap.py -fdrlimit -10 -A $TREATMENT $CONTROL $TREAT_OUT $CON_OUT

#usage: wellington_bootstrap.py [-h] [-fp FOOTPRINT_SIZES] [-fdr FDR_CUTOFF]
 #                              [-fdriter FDR_ITERATIONS] [-fdrlimit FDR_LIMIT]
 #                              [-p PROCESSES] [-A]
 #                              treatment_bam control_bam bedsites
 #                              treatment_only_output control_only_output
#
#Scores Differential Footprints using Wellington-Bootstrap.

#positional arguments:
#  treatment_bam         BAM file for treatment
#  control_bam           BAM file for control
#  bedsites              BED file of genomic locations to search in
#  treatment_only_output
#                        File to write treatment specific fooprints scores to
#  control_only_output   File to write control specific footprint scores to

#optional arguments:
#  -h, --help            show this help message and exit
#  -fp FOOTPRINT_SIZES, --footprint-sizes FOOTPRINT_SIZES
#                        Range of footprint sizes to try in format
#                        "from,to,step" (default: 11,26,2)
#  -fdr FDR_CUTOFF, --FDR_cutoff FDR_CUTOFF
#                        Detect footprints using the FDR selection method at a
#                        specific FDR (default: 0.01)
#  -fdriter FDR_ITERATIONS, --FDR_iterations FDR_ITERATIONS
#                        How many randomisations to use when performing FDR
#                        calculations (default: 100)
#  -fdrlimit FDR_LIMIT, --FDR_limit FDR_LIMIT
#                        Minimum p-value to be considered significant for FDR
#                        calculation (default: -20)
#  -p PROCESSES, --processes PROCESSES
#                        Number of processes to use (default: uses all CPUs)
#  -A                    ATAC-seq mode (default: False)