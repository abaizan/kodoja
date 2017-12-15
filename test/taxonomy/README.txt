The test suite will download and cache NCBI taxonomy files here for use
with the small provided sample databases. This is to avoid including
the relatively large NCBI taxonomy files under git version control.

e.g.

$ wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
$ tar zxvf taxdump.tar.gz
