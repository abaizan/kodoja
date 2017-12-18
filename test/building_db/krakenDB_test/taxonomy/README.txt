This folder is a cut-down NCBI taxonomy for use with kraken.

names.dmp and nodes.dmp are symlinks to the files we download
once under TravisCI or for local testing, see README in the
folder test/taxonomy

nucl_*.accession2taxid are filtered versions of the very large
files from the NCBI using our three test viruses only.

accmap.dlflag, taxdump.dlflag and taxdump.untarflag are flag
files to trick kraken into not downloading and unzipping things.
