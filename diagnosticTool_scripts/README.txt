Diagnostic pipeline for the detection of viral sequences in RNA-seq data
------------------------------------------------------------------------

This pipeline was made for the detection of virus sequences in a metagenomic plant RNA-seq dataset. 
It takes the raw data (either fasta or fastq) and uses k-mer-based tools (as well as other intermediate tools)
to detect viral sequences within the sample.

There are three main scripts to run this pipeline on the cluster:
diagnostic_master.sh - this script is to qsub the master python script
diagnostic_master.py - this script contains the parameters and calls the modules to run the pipeline
diagnostic_module.py - this contains all the functions called by diagnostic_master.py

There are also scripts to make a custom k-mer database for kraken:
database_master.sh - this script is to qsub the master python script
database_master.py - this script contains the parameters and calls the modules to create a new database
database_modules.py - contains all the functions called by database_master.py
