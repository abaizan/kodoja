[![BioConda package](https://img.shields.io/conda/vn/bioconda/kodoja.svg)](https://anaconda.org/bioconda/kodoja)
[![Build Status](https://travis-ci.org/abaizan/kodoja.svg?branch=master)](https://travis-ci.org/abaizan/kodoja)

# Kodoja: A workflow for virus detection in plants using k-mer analysis of RNA-sequencing data

Kodoja takes the raw data (either fasta or fastq) and uses Kraken, a k-mer-based tool,
and Kaiju, which used the Burrows–Wheeler transform, to detect viral sequences in RNA-seq or sRNA-seq data.

## Overview
There are three main scripts:

* ``kodoja_search.py`` - classify RNA-seq data.
* ``kodoja_build.py`` - download viral/host genomes and create new Kraken and Kaiju databases.
* ``kodoja_retrieve.py`` - pull out sequences of interest from kodoja_search.py results file.

Python files ``diagnostic_modules.py`` and ``database_modules.py`` contain the fuctions called by kodoja_search.py and kodoja_build.py
The ``.sh`` files are example script for submission to SGE cluster

## License

Kodoja is released under the MIT licence, see file ``LICENSE.txt`` for details.

## Dependencies

* FastQC v0.11.5
* Trimmomatic v0.36
* Kraken v1.0
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

IMPORTANT: do not put original data in the output directory when executing kodoja_search!

## kodoja_search.py parameters:
### General:
* ``--read1`` - path to the single-end or first paired-end file (required)
* ``--data_format`` - specify the file-type for file1 ("fasta" or "fastq" - default='fastq')
* ``--output_dir`` - path to the results folder (required)
* ``--threads`` - number of threads on cluster (default=1)
* ``--read2`` - path to second paired-end file (default=False)
*``--host_subset`` - tax id of host. Use this is a host genome was added to the
databases and you do not wish to see the number of reads classifed to this group
in the final table

### Kraken:
* ``--kraken_db`` - path to kraken database (required)
* ``--kraken_quick`` - Quick operation mode of Kraken, where instead of querying all
  k-mers in the database, it stops at nth k-mer hit preload (default=False)

### Kaiju:
* ``--kaiju_db`` - path to kaiju database, nodes.dmp and names.dmp files (required)
* ``--kaiju_minlen`` - minimun required fragment length length (default=15)
* ``--kaiju_mismatch`` - number of mismatches allowed by kaiju (default=1)
* ``--kaiju_score`` - minimum required match if mismatches introduced (default=85)

Set parameter for kaiju:
``-x`` -  used to enable filtering of query sequences
    containing low-complexity regions by using the SEG algorithm from the blast+
    package. Enabling this option is always recommended in order to avoid false
    positive matches caused by spurious matches due to simple repeat patterns or
    other sequencing noise.

### Trimmomatic:
* ``--trim_minlen`` - minimum length read after trimming (default=50)
* ``--trim_adapt`` - fasta file with Illumina adaptor sequences to allow trimming (default=False)

Set parameters for trimmomatic
``ILUMINACLIP 2:30:10`` (seed mismatches:palindrome threshold:simple clip threshold) -
     seedMismatches specifies the maximum mismatch count which will still allow a full match to be performed,
     palindromeClipThreshold specifies how accurate the match between the two 'adapter ligated',
        reads must be for PE palindrome read alignment,
        simpleClipThreshold: specifies how accurate the match between
	any adapter etc. sequence must be against a read.
``LEADING:20`` -  Specifies the minimum quality required to keep a base
``TRAILING:20`` - Specifies the minimum quality required to keep a base

## kodoja_build.py parameters:
### General parameters:
* ``--output_dir`` - Output directory path where kraken and kaiju databases will be written, required')
* ``--threads`` - number of threads on cluster (default=1)
* ``--host`` - NCBI tax id for the host genome to be downloaded from refseq and
added to the databases(default=False)
* ``--extra_files`` - List of file paths (default=False)
* ``--extra_taxids`` - List of tax ids corresponding to extra files (default=False)
* ``--all_viruses`` - Build databases with viruses from all hosts
* ``--db_tag`` - Suffix for databases (default=none)

### Kraken database:
* ``--kraken_kmer`` - Kraken kmer size type=int, (default=31)
* ``--kraken_minimizer`` - Kraken minimizer size (default=15)

### ncbi-genome-download:
* ``--download_parallel`` - number of genomes to download in parallel (default=4)
* ``--no_download`` - Genomes have already been downloaded and are in output folder (default=False)

## kodoja_retrieve.py parameters:
* ``--file_dir`` - Path to directory of kodoja_search results (required)
* ``--user_format`` - Sequence data format (default=fastq)
* ``--read1`` - Path to read 1 file (required)
* ``--read2`` - Path to read 2 file
* ``--taxID' - Virus tax ID for subsetting (default: All viral sequences)
* ``--genus`` - Include sequences classified at the genus level in subset file
* ``--stringent`` - Only subset sequences identified to same virus by both tools

## Release History

| Version | Date       | Notes                                               |
| ------- | ---------- | --------------------------------------------------- |
| 0.0.1   | 2018-01-15 | - Initial release for BioConda packaging            |
| 0.0.2   | 2018-01-22 | - Now tested under Python 3.6 as well as Python 2.7 |
| 0.0.3   | 2018-02-22 | - Include genus level counts in search results      |
|         |            | - Simplify internal renaming of sequencing reads    |
| 0.0.4   | *Pending*  | - Code style updates (no functional changes)        |
|         |            | - Provide cut-down NCBI taxonomy for tests cases    |
|         |            | - Additional database build testing                 |
|         |            | - Downloads virus files with HTTPS rather than FTP  |


## Development

Kodoja is on GitHub, and has auotmated testing running on TravisCI, see special
file ``.travis.yml`` and webpage https://travis-ci.org/abaizan/kodoja/builds
for details.

The release process includes:

1. Update version in ``diagnosticTool_scripts/diagnostic_modules.py``.
1. Update release history in this ``README.md`` file.
3. Commit changes.
4  Tag the commit with ``git tag kodoja-vX.Y.X``
5. Push commits and tags to github with ``git push origin master --tags``
6. Submit a pull request to BioConda to update the package, which usally
   just means bumping the verion and updating the checksum in ``meta.yaml``:
   https://github.com/bioconda/bioconda-recipes/tree/master/recipes/kodoja
