#!/bin/bash

#$ -cwd
#$ -j yes
#$ -pe smp 1

python /home/ae42909/viral_diagnostics/parameter_testing/parameterTest_analysis/kraken_results/run_results.py
