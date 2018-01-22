[![Build Status](https://travis-ci.org/abaizan/kodoja.svg?branch=master)](https://travis-ci.org/abaizan/kodoja)

# Kodoja: A workflow for virus detection in plants using k-mer analysis of RNA-sequencing data

Kodoja takes the raw data (either fasta or fastq) and uses Kraken, a k-mer-based tool,
and Kaiju, which used the Burrows–Wheeler transform, to detect viral sequences in RNA-seq or sRNA-seq data.

## Overview
There are three main scripts:

* kodoja_search.py - to classify RNA-seq data.
* kodoja_build.py - to download viral/host genomes and create new Kraken and Kaiju databases.
* kodoja_retrieve.py - pull out sequences of interest from kodoja_search.py results file.

diagnostic_modules.py and database_modules.py contain the fuctions called by diagnostic_master and database_master.
.sh files are example script for submission to SGE cluster

## License

Kodoja is released under the MIT licence, see file ``LICENSE.txt`` for details.

## Dependencies

* FastQC v0.11.5,
* Trimmomatic v0.36,
* Kraken v1.0,
* Kaiju v1.5.0

Python packages:
* numpy
* biopython
* pandas
* ncbi-genome-download

## Installation

A conda package has been prepared on the BioConda channel which will install Kodoja and a known
working combination of the dependencies, all with just:

```
$ conda install -c bioconda kodoja
```

For manual installation, you must install all the dependencies by hand and then add the main
scripts folder to your ``$PATH`` so that you can run ``kodoja_search.py`` etc at the command
line.

## Usage
```
here
```

IMPORTANT: When executing the pipeline, do not put original data in the result
file before executing script


## Classification pipeline parameters:
### General:
--read1 - path to the single-end or first paired-end file (required)
--data_format - specify the file-type for file1 ("fasta" or "fastq" - default='fastq')
--output_dir - path to the results folder (required)
--threads - number of threads on cluster (default=1)
--read2 - path to second paired-end file (default=False)

### Kraken:
--kraken_db - path to kraken database (required)
--kraken_quick - Quick operation mode of Kraken, where instead of querying all k-mers in the
                    database, it stops at nth k-mer hit preload (default=False)
### Kaiju:
--kaiju_db - path to kaiju database, nodes.dmp and names.dmp files (required)
--kaiju_minlen - minimun required fragment length length (default=15)
--kaiju_mismatch - number of mismatches allowed by kaiju (default=1)
--kaiju_score - minimum required match if mismatches introduced (default=85)
Set parameters:
-x: used to enable filtering of query sequences
    containing low-complexity regions by using the SEG algorithm from the blast+
    package. Enabling this option is always recommended in order to avoid false
    positive matches caused by spurious matches due to simple repeat patterns or
    other sequencing noise.

### Trimmomatic:
--trim_minlen - minimum length read after trimming (default=50)
--trim_adapt - fasta file with Illumina adaptor sequences to allow trimming (default=False)
Set parameters:
ILUMINACLIP 2:30:10 <seed mismatches>:<palindrome threshold>:<simple clip threshold>
     seedMismatches: specifies the maximum mismatch count which will still allow a full match to be performed
     palindromeClipThreshold: specifies how accurate the match between the two 'adapter ligated'
        reads must be for PE palindrome read alignment.
        simpleClipThreshold: specifies how accurate the match between
	any adapter etc. sequence must be against a read.
LEADING:20 Specifies the minimum quality required to keep a base
TRAILING:20 Specifies the minimum quality required to keep a base

## Database construction parameters:
--output_dir - Output directory path where kraken and kaiju databases will be written, required')
--threads - number of threads on cluster (default=1)
--test - Make database for test_script.py
--host - Host tax ID (default=False)
--extra_files - List of file names added to "extra" dir (default=False)
--extra_taxids - List of tax ids corresponding to extra files (default=False)
--all_viruses - Build databases with all viruses (defult=plant viruses only)
--db_tag - Suffix for databases (default=none)

### Kraken database:
--kraken_kmer - Kraken kmer size type=int, (default=31)
--kraken_minimizer - Kraken minimizer size (default=15)

### ncbi-genome-download:
--download_parallel - number of genomes to download in parallel (default=4)
--no_download - Genomes have already been downloaded and are in output folder (default=False)


## Release History

| Version | Date       | Notes                                               |
| ------- | ---------- | --------------------------------------------------- |
| 0.0.1   | 2018-01-15 | Initial release for BioConda packaging              |

## Development

Kodoja is on GitHub, and has auotmated testing running on TravisCI, see special
file ``.travis.yml`` and webpage https://travis-ci.org/abaizan/kodoja/builds
for details.
