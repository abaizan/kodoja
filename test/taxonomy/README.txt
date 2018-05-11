The full NCBI taxonomy files are large, but can be downloaded using

$ wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
$ tar zxvf taxdump.tar.gz

Or, using Kraken to download and setup a taxonomy folder as part of its
database building script - note the inclusion of flag files like
accmap.dlflag to avoid re-downloading the accesion maps again.

$ kraken-build --download-taxonomy --db .
$ ls taxonomy/ -h -l
total 25G
-rw-r--r-- 1 ae42909 cms    0 Dec 18 14:39 accmap.dlflag
-rw-r--r-- 1 ae42909 cms  16M Dec 18 14:20 citations.dmp
-rw-r--r-- 1 ae42909 cms 3.1M Dec 18 14:20 delnodes.dmp
-rw-r--r-- 1 ae42909 cms  442 Dec 18 14:20 division.dmp
-rw-r--r-- 1 ae42909 cms  15K Dec 18 14:20 gc.prt
-rw-r--r-- 1 ae42909 cms 4.5K Dec 18 14:20 gencode.dmp
-rw-r--r-- 1 ae42909 cms 865K Dec 18 14:20 merged.dmp
-rw-r--r-- 1 ae42909 cms 141M Dec 18 14:20 names.dmp
-rw-r--r-- 1 ae42909 cms 109M Dec 18 14:20 nodes.dmp
-rw-r--r-- 1 ae42909 cms 2.6G Dec 18 14:35 nucl_est.accession2taxid
-rw-r--r-- 1 ae42909 cms 4.6G Dec 18 14:37 nucl_gb.accession2taxid
-rw-r--r-- 1 ae42909 cms 1.4G Dec 18 14:37 nucl_gss.accession2taxid
-rw-r--r-- 1 ae42909 cms  17G Dec 18 14:39 nucl_wgs.accession2taxid
-rw-r----- 1 ae42909 cms 2.6K Jun 13  2006 readme.txt
-rw-r--r-- 1 ae42909 cms    0 Dec 18 14:39 taxdump.dlflag
-rw-r--r-- 1 ae42909 cms  40M Dec 18 14:39 taxdump.tar.gz
-rw-r--r-- 1 ae42909 cms    0 Dec 18 14:42 taxdump.untarflag

To avoid including these large files under version control, we
instead only include a subset filtered according to the limited
set of (viral) entries used in our test cases. First:

$ cd /tmp
$ wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
$ tar -zxvf taxdump.tar.gz

Then, returning to this directory:

$ ./filter_taxonomy.py 137758 946046 12227
Filtering NCBI taxonomy files /tmp/nodes.dmp and names.dmp etc
Will create ./nodes.dmp and ./names.dmp etc using just the given
3 entries and their parent nodes.
Loaded 1755515 entries from /tmp/nodes.dmp
Expanded 3 given TaxID to a list of 10 including ancestors
Filtering merged.dmp
Filtering names.dmp
Filtering nodes.dmp
Checking which genetic codes and divisions are needed
Only want these genetic codes: 1
Only want these divisions: 8,9
Filtering gencode.dmp
Filtering division.dmp
Filtering citations.dmp
Done
