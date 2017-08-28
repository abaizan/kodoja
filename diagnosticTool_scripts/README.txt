Diagnostic pipeline for the detection of viral sequences in RNA-seq data
------------------------------------------------------------------------

This pipeline was made for the detection of virus sequences in a metagenomic 
plant RNA-seq dataset. 

It takes the raw data (either fasta or fastq) and uses k-mer-based tools to 
detect viral sequences within the sample. If the raw data is in fastq format, 
a QC and trimming script will be implemented. There are also scripts to make 
a custom k-mer database for kraken and kaiju, which downloads and uses all 
ncbi refseq viral data (genomic for kraken and protein for kaiju). 

Tools used as part of the pipeline and for database construction: fastqc 
(0.10.0), trimmomatic (0.32), kraken (0.10.5-beta), kaiju (1.5.0)
Linux tools: gzip/gunzip, tar, 

Scripts are written in python (run on python 2.7.13 here)
Python packages needed: subprocess, pandas, numpy, os, urllib, re, 
ncbi-genome-download and SeqIO (from biopython).

There are three main scripts to run this pipeline on the cluster:
diagnostic_master.sh - this script is to qsub the master python script
diagnostic_master.py - this script contains the parameters and calls the 
                       modules to run the pipeline
diagnostic_module.py - this contains all the functions called by 
                       diagnostic_master.py

The scripts to make a custom k-mer database for kraken and kaiju:
database_master.sh - this script is to qsub the master python script
database_master.py - this script contains the parameters and calls the 
                     modules to create a new database
database_modules.py - contains all the functions called by database_master.py

# Diagnostic pipeline parameters:
General:
    file1 - path to the single-end or first paired-end file
    user_format - specify the file-type for file1 ("fasta" or "fastq")
    threads - number of cluster nodes (e.g. "4")
    out_dir - path to the folder where result files will end up (doesn't need 
             to exits, will get made as part of pipeline)
    file2 - False unless changed to contain path to second paired-end data file
Trimmomatic:
    trim_minlen - minimun length read after trimming (drop reads below this 
		number e.g. "50")
    adapter_file - fasta file with Illumina adaptor sequences to allow trimming 
	       thereof
    ILUMINACLIP 2:30:10 (this parameter can't be changed)
		<seed mismatches>:<palindrome threshold>:<simple clip threshold> 
		seedMismatches: specifies the maximum mismatch count which will
			still allow a full match to be performed
    		palindromeClipThreshold: specifies how accurate the match 
			between the two 'adapter ligated' reads must be for PE 
			palindrome read alignment.
		simpleClipThreshold: specifies how accurate the match between 
		any adapter etc. sequence must be against a read.
    LEADING:20 (this parameter can't be changed)
		Specifies the minimum quality required to keep a base
    TRAILING:20 (this parameter can't be changed)
		Specifies the minimum quality required to keep a base
Kraken:
    kraken_db - path to kraken database
Kaiju:
    kaiju_nodes - path to nodes.dmp
    kaiju_names - path to names.dmp
    kaiju_fmi - path to kaiju database (.fmi file)
    kaiju_missmatch - number of mismatches allowed by kaiju (e.g. "5")

# Database construction parameters:
    tool - specify the tool for which a datbase is being made (either "kraken" 
	   or "kaiju")
ncbi-genome-download:
    ncbi_download_parallel - how many files to download in parallel (e.g. "4")
Kraken:
    threads -  number of cluster nodes (e.g. "10")
    kraken_kmer - k-mer size for kraken database fragments 
    kraken_minimizer - 
    kraken_db_dir - path to folder for database files
    jellyfish_hash_size - 
Kaiju:
    kaiju_db_dir - path to folder for database files
