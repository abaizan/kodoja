Diagnostic pipeline for the detection of viral sequences in RNA-seq data
------------------------------------------------------------------------

This pipeline was made for the detection of virus sequences in a metagenomic
plant RNA-seq dataset.

It takes the raw data (either fasta or fastq) and uses Kraken, a k-mer-based tool,
and Kaiju which used the Burrows–Wheeler transform, to
detect viral sequences within the sample. If the raw data is in fastq format,
a QC and trimming script will be implemented. There are also scripts to make
a custom database for kraken and kaiju, which downloads and uses all ncbi refseq
viral data (genomic for kraken and protein for kaiju). You can also add extra files
to the database by downloading the reference database you want to use and put it into
a folder called "extra" in your genome_download_dir location (see database construction
below). These files need to be  and if it
doesn't have already, change the file extension to either ".fna" for genomic or
".faa" for protein

Tools used as part of the pipeline and for database construction: fastqc
(0.10.0), trimmomatic (0.32), kraken (0.10.5-beta), kaiju (1.5.0), jellyfish
(1.1.11)
Linux tools: gzip/gunzip, tar

Scripts are written in python (run on python 2.7.13 here)
Python packages needed: subprocess, pandas, numpy, os, urllib, re,
ncbi-genome-download and SeqIO (from biopython).

IMPORTANT: When executing the pipeline, do not put original data in the result
file before executing script

There are three main scripts to run this pipeline on the cluster:
diagnostic_master.sh - this script is to qsub the master python script
diagnostic_master.py - this script contains the parameters and calls the
                       modules to run the pipeline
diagnostic_modules.py - this contains all the functions called by
                       diagnostic_master.py

The scripts to make a custom database for kraken and kaiju:
database_master.sh - this script is to qsub the master python script
database_master.py - this script contains the parameters and calls the
                     modules to create a new database
database_modules.py - contains all the functions called by database_master.py

#############################################################################
# Diagnostic pipeline parameters:
General:
    file1 - path to the single-end or first paired-end file
    user_format - specify the file-type for file1 ("fasta" or "fastq")
    threads - number of cluster nodes (e.g. "4")
    out_dir - path to the folder where result files will end up (doesn't need
             to exits, will get made as part of pipeline)
    file2 - False unless changed to contain path to second paired-end data file
Trimmomatic:
    trim_minlen - minimum length read after trimming (drop reads below this
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
    quick_minhits - number specifying number of hits (n). Quick operation mode
                    of Kraken, where instead of querying all k-mers in the
                    database, it stops at nth k-mer hit preload - binary on/off
Kaiju:
    kaiju_db - path to kaiju database, nodes.dmp and names.dmp files
    kaiju_minlen - minimun required fragment length length (e.g. 11)
    kaiju_mismatch - number of mismatches allowed by kaiju (e.g. 5)
    kaiju_score - minimum required match if mismatches introduced (e.g. 65)
    -x (can't be changed) is used to enable filtering of query sequences
    containing low-complexity regions by using the SEG algorithm from the blast+
    package. Enabling this option is always recommended in order to avoid false
    positive matches caused by spurious matches due to simple repeat patterns or
    other sequencing noise.

#############################################################################
# Database construction parameters:
    tool - specify the tool for which a datbase is being made (either "kraken"
	   or "kaiju")
	genome_download_dir - path to the file where all reference genomes will be
        downloaded
	extra_files - list of the names (strings) of the files you want to add to
        either the kraken or kaiju databases. They have to be in the fasta
        format and their extensions have to be either ".fna" or ".faa" and they
        have to be compressed (eg. ".fna.gz").
	extra_taxid - list of each ncbi taxonomic id (numeric) for the species in
        extra files. They have to be in order (ie. if
        extra_file = ["species1.fna.gz", "species2.fna.gz"] then
        extra_taxid = [taxid_sp1, taxid_sp2]))
ncbi-genome-download:
    ncbi_download_parallel - how many files to download in parallel (e.g. 4)
Kraken:
    threads -  number of cluster nodes (e.g. 4)
    kraken_kmer - k-mer size for kraken database fragments
    kraken_minimizer -
    kraken_db_dir - path to folder for refseq database files
    jellyfish_hash_size -
Kaiju:
    kaiju_db_dir - path to folder for refseq database files

#############################################################################
# Parameter testing:
See README.txt in /home/ae42909/Scratch/parameter_test
