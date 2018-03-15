This folder is a cut-down NCBI taxonomy for use with kraken.

names.dmp and nodes.dmp etc are symlinks to the cut-down files
based on the NCBI taxonomy for the three test viruses only.

nucl_*.accession2taxid are filtered versions of the very large
files from the NCBI using our three test viruses only.

accmap.dlflag, taxdump.dlflag and taxdump.untarflag are flag
files to trick kraken into not downloading and unzipping things.
