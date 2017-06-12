#!/bin/bash

#$ -cwd
#$ -j yes
#$ -pe smp 4
THREADS=4
PEAKS=peak_intersect/pyDNase_EBF1
BAM=bam_data/RAGko_1_sorted.bam
DIROUT=wellington/RAGko1_EBF1

#optional arguments:
#  -h, --help            show this help message and exit
#  -b, --bonferroni      Performs a bonferroni correction (default: False)
#  -sh SHOULDER_SIZES, --shoulder-sizes SHOULDER_SIZES
#                        Range of shoulder sizes to try in format
#                        "from,to,step" (default: 35,36,1)
#  -fp FOOTPRINT_SIZES, --footprint-sizes FOOTPRINT_SIZES
#                        Range of footprint sizes to try in format
#                        "from,to,step" (default: 11,26,2)
#  -d, --one_dimension   Use Wellington 1D instead of Wellington (default:
#                        False)
#  -fdr FDR_CUTOFF, --FDR_cutoff FDR_CUTOFF
#                        Write footprints using the FDR selection method at a
#                        specific FDR (default: 0.01)
#  -fdriter FDR_ITERATIONS, --FDR_iterations FDR_ITERATIONS
#                        How many randomisations to use when performing FDR
#                        calculations (default: 100)
#  -fdrlimit FDR_LIMIT, --FDR_limit FDR_LIMIT
#                        Minimum p-value to be considered significant for FDR
#                        calculation (default: -20)
#  -pv PV_CUTOFFS, --pv_cutoffs PV_CUTOFFS
#                        Select footprints using a range of pvalue cutoffs
#                        (default: -10,-20,-30,-40,-50,-75,-100,-300,-500,-700
#  -dm, --dont-merge-footprints
#                        Disables merging of overlapping footprints (Default:
#                        False)
#  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
#                        The prefix for results files (default:
#                        <reads.regions>)
#  -p P                  Number of processes to use (default: uses all CPUs)
#  -A                    ATAC-seq mode (default: False)

wellington_footprints.py -A $PEAKS $BAM $DIROUT
