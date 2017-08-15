#!/bin/bash

#$ -cwd
#$ -j yes
#$ -pe smp 2

#Build database including complete genomes from NCBI RefSeq for bacteria (-r) and virus (-v) = ~14GB. You can also have a more complete set of bacterial genomes from proGenomes (http://progenomes.embl.de/, -p option) = ~13GB. You can also download the nr database that is used by NCBI BLAST and extract proteins belonging to Archaea, Bacteria and Viruses (-n), or as well as those for option "-n", you can include  proteins from fungi and microbial eukaryotes (-e) = ~43GB

/home/ae42909/Programs/Kaiju/kaiju-v1.5.0-linux-x86_64-static/bin/makeDB.sh -r -v

# After makeDB.sh is finished, only the files kaiju_db.fmi (or kaiju_db_nr.fmi / kaiju_db_nr_euk.fmi), nodes.dmp, and names.dmp are needed to run Kaiju. The remaining files and the genomes/ directory containing the downloaded genomes can be deleted.
