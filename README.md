[![BioConda package](https://img.shields.io/conda/vn/bioconda/kodoja.svg)](https://anaconda.org/bioconda/kodoja)
[![Build Status](https://travis-ci.org/abaizan/kodoja.svg?branch=master)](https://travis-ci.org/abaizan/kodoja)

# Kodoja: A workflow for virus detection in plants using k-mer analysis of RNA-sequencing data

Kodoja takes the raw data (either fasta or fastq) and uses Kraken, a k-mer-based tool,
and Kaiju, which used the Burrows–Wheeler transform, to detect viral sequences in RNA-seq or sRNA-seq data.

## Overview
There are three main scripts:

* ``kodoja_search.py`` - classify RNA-seq data.
* ``kodoja_build.py`` - download viral/host genomes and create new Kraken and Kaiju databases.
* ``kodoja_retrieve.py`` - pull out sequences of interest from ``kodoja_search.py`` results file.

Python files ``diagnostic_modules.py`` and ``database_modules.py`` contain the fuctions called by ``kodoja_search.py`` and ``kodoja_build.py``, and are not intended for public use.

The ``.sh`` files are example scripts for submission to an SGE cluster.

## License

Kodoja is released under the MIT licence, see file ``LICENSE.txt`` for details.

## Dependencies

The lower versions listed were those used in the initial development and/or
local testing of Kodoja. Later updates will likely work unless the tool makes
a backward incompatible change.

* FastQC v0.11.5
* Trimmomatic v0.36
* Kraken v1.0
* Kaiju v1.5.0

Python packages:

* numpy v1.9
* biopython v1.67
* pandas v0.14
* ncbi-genome-download v0.2.6

You can use Python 2.7 or Python 3, specifically Kodoja has been tested
on Python 3.6.

## Installation

A conda package has been prepared on the BioConda channel which will install
Kodoja and the dependencies, all with just:

```
$ conda install -c bioconda kodoja
```

For manual installation, you must install all the dependencies by hand and then add the main
scripts folder to your ``$PATH`` so that you can run ``kodoja_search.py`` etc at the command
line.

## Pre-built Databases

You can	use ``kodoja_build.py``	to make	your own databses, or download the
pre-built database as described here.

The kodojaDB v1.0 was released Sept 2018 under the CC-BY 4.0 license. It can
be downloaded and cited as https://doi.org/10.5281/zenodo.1406071 (where the
metadata describes how it was made). We suggest you install it as follows:

``` bash
$ cd /mnt/shared/data/
$ mkdir kodojaDB_v1.0
$ cd kodojaDB_v1.0
$ wget https://zenodo.org/record/1406071/files/kodojaDB_v1.0.tar.gz
$ tar -zxvf kodojaDB_v1.0.tar.gz
```

You would then use this with ``kodoja_search.py`` as follows:

``` bash
$ kodoja_search.py --kraken_db /mnt/shared/data/kodojaDB_v1.0/krakenDB \
                   --kaiju_db /mnt/shared/data/kodojaDB_v1.0/kaijuDB \
		   ...
```

## Usage

IMPORTANT: do not put original data in the output directory when executing kodoja_search!

## kodoja_search.py parameters:
### General:
* ``--read1`` - path to the single-end or first paired-end file (required)
* ``--read2`` - path to second paired-end file (default=False)
* ``--data_format`` - specify the file-type for file1 ("fasta" or "fastq" - default='fastq')
* ``--output_dir`` - path to the results folder (required)
* ``--threads`` - number of threads on cluster (default=1)
* ``--host_subset`` - tax id of host. Use this is a host genome was added to the
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
* ``--taxID`` - Virus tax ID for subsetting (default: All viral sequences)
* ``--genus`` - Include sequences classified at the genus level in subset file
* ``--stringent`` - Only subset sequences identified to same virus by both tools

## Release History

| Version | Date       | Notes                                               |
| ------- | ---------- | --------------------------------------------------- |
| 0.0.1   | 2018-01-15 | - Initial release for BioConda packaging            |
| 0.0.2   | 2018-01-22 | - Now tested under Python 3.6 as well as Python 2.7 |
| 0.0.3   | 2018-02-22 | - Include genus level counts in search results      |
|         |            | - Simplify internal renaming of sequencing reads    |
| 0.0.4   | 2018-08-22 | - Code style updates (no functional changes)        |
|         |            | - Provide cut-down NCBI taxonomy for tests cases    |
|         |            | - Additional database build testing                 |
|         |            | - Downloads virus files with HTTPS rather than FTP  |
| 0.0.5   | 2018-08-29 | - Refactor logging in ``kodoja_search.py``          |
|         |            | - Top level error handling, with logging in search  |
|         |            | - ``dictionary changed size during iteration`` bug  |
| 0.0.6   | 2018-09-04 | - Python 3 fix for ``kodoja_retrieve.py``           |
|         |            | - Automated testing of ``kodoja_retrieve.py``       |
|         |            | - Also test paired reads without /1 and /2 suffixes |
| 0.0.7   | 2018-09-07 | - Document installing prebuilt database from Zenodo |
|         |            | - Optimise sorting of pandas dataframes             |
|         |            | - Zero not blank in cols 6 and 7 of virus_table.txt |
|         |            | - Automated testing of pinned & latest dependencies |


## Development

Kodoja is on GitHub, and has auotmated testing running on TravisCI, see special
file ``.travis.yml`` and webpage https://travis-ci.org/abaizan/kodoja/builds
for details.

The release process includes:

1. Update version in ``diagnosticTool_scripts/diagnostic_modules.py``.
1. Update release history in this ``README.md`` file.
3. Commit changes.
4  Tag the commit with ``git tag kodoja-vX.Y.Z``
5. Push commits and tags to github with ``git push origin master --tags``
6. Submit a pull request to BioConda to update the package, which usally
   just means bumping the version and updating the checksum in ``meta.yaml``:
   https://github.com/bioconda/bioconda-recipes/tree/master/recipes/kodoja
