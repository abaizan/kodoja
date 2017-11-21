#!/bin/bash

#$ -cwd
#$ -j yes
#$ -pe smp 1

python /home/ae42909/viral_diagnostics/diagnosticTool_scripts/diagnostic_seqRetrive.py
