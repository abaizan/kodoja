Parameter testing:
--------------------
Optimisation of kraken and kaiju. There are three files to run different
  parameters
    _parameterValues.py - values for parameter testing
    _parameterTest.py - import parameter values and functions to run the tool
    _parameterTest.sh - script to run python script on the cluster
# Kraken
There are two '_parameterTest.py' scripts for kraken, on to make kraken
  databases ('krakenDB_') and one to run kraken ('kraken_'). They both import
  the values they will use from 'kraken_parameterValues.py'.

krakenDB parameters: kraken_kmer, kraken_minimizer
Each kraken database will be contanied within a folder with the k-mer size (k)
  and minimizer size (m) (e.g. krakenDB_k10_m3)

kraken parameters: quick_minhits, preload
When running kraken_parameterTest.py: the code will look for all directories
  that start with "kraken" in "kraken_db_dir" so add a 'no_' as the first 3
  characters of the directories you don't want to use

# Kaiju
kaiju parameters: kaiju_minlen, kaiju_mismatch, kaiju_score

# Results of parameter testing in
  /home/ae42909/Scratch/parameter_test/<tool>

Each run has a logfile and is either different parameters or different initial
  datafiles

kraken/
  run1
    input1 = /home/ae42909/Scratch/100_Potato_withViruses_1.fastq
    input2 = /home/ae42909/Scratch/100_Potato_withViruses_2.fastq
    krakenDB_parameters = {"kraken_kmer":[15,21,31], "kraken_minimizer":[5, 15]}
    kraken_parameters = {"quick_minhits":[False, 1, 2, 3, 4, 5],
                            "preload":[True, False]}
  run2
    input1 = /home/ae42909/synthetic_RNAseq/mappingRNAseq/
              concatenated_fastaFiles/Potato_withViruses_1.fastq
    input2 = /home/ae42909/synthetic_RNAseq/mappingRNAseq/
              concatenated_fastaFiles/Potato_withViruses_2.fastq
    krakenDB_parameters = {"kraken_kmer":[15,21,31], "kraken_minimizer":[5, 15]}
    kraken_parameters = {"quick_minhits":[False, 1, 2, 3, 4, 5],
                            "preload":[True, False]}
  run3
    input1 = '/home/ae42909/data_forTesting/Barerro_data/
                PB64-S1_clean.fq.gz_trim.fastq'
    krakenDB_parameters = {"kraken_kmer":[5,10,15], "kraken_minimizer":[2, 3, 5]}
    kraken_parameters = {"quick_minhits":[False], "preload":[False]}
  run4
    input1 = '/home/ae42909/data_forTesting/Barerro_data/
                  PB64-S1_clean.fq.gz_trim.fastq'
    krakenDB_parameters = {"kraken_kmer":[15, 16, 17, 18], "kraken_minimizer":[5]}
    kraken_parameters = {"quick_minhits":[False], "preload":[False]}
  run5
    input1 = '/home/ae42909/data_forTesting/Barerro_data/
                PB64-S2_clean.fq.gz_trim.fastq'
    krakenDB_parameters = {"kraken_kmer":[15, 16, 17, 18], "kraken_minimizer":[5]}
    kraken_parameters = {"quick_minhits":[False], "preload":[False]}
  run6
    input1 = '/home/ae42909/data_forTesting/Barerro_data/
                PB64-S3_clean.fq.gz_trim.fastq'
    krakenDB_parameters = {"kraken_kmer":[15, 16, 17, 18], "kraken_minimizer":[5]}
    kraken_parameters = {"quick_minhits":[False], "preload":[False]}
  run7
    input1 = '/home/ae42909/data_forTesting/Barerro_data/
                PB64-S4_clean.fq.gz_trim.fastq'
    krakenDB_parameters = {"kraken_kmer":[15, 16, 17, 18], "kraken_minimizer":[5]}
    kraken_parameters = {"quick_minhits":[False], "preload":[False]}
  run8
    input1 = '/home/ae42909/data_forTesting/Barerro_data/
                PB64-S5_clean.fq.gz_trim.fastq'
    krakenDB_parameters = {"kraken_kmer":[15, 16, 17, 18], "kraken_minimizer":[5]}
    kraken_parameters = {"quick_minhits":[False], "preload":[False]}
  run9
    input1 = '/home/ae42909/data_forTesting/Barerro_data/
                PB64-S7_clean.fq.gz_trim.fastq'
    krakenDB_parameters = {"kraken_kmer":[15, 16, 17, 18], "kraken_minimizer":[5]}
    kraken_parameters = {"quick_minhits":[False], "preload":[False]}
  run10
    input1 = '/mnt/shared/projects/virology/201609_BBSRC_Diagnostics/Data/
              synthetic/review-paper-test-datasets/SRR1123893Pepper.fastq'
    kraken k-mers = 15 (m5), 18 (m5), 21 (m15), 31 (m15)
    kraken_parameters = {"quick_minhits":[False], "preload":[False]}
  run11
    input1 = '/mnt/shared/projects/virology/201609_BBSRC_Diagnostics/Data/
              synthetic/review-paper-test-datasets/SRR1269627Pear.fastq'
    kraken k-mers = 15 (m5), 18 (m5), 21 (m15), 31 (m15)
    kraken_parameters = {"quick_minhits":[False], "preload":[False]}
  run12
    input1 = '/mnt/shared/projects/virology/201609_BBSRC_Diagnostics/Data/
              synthetic/review-paper-test-datasets/SRR1123893Pepper.fastq'
    kraken k-mers =  23 (m15), 27 (m15), 29 (m15)
    kraken_parameters = {"quick_minhits":[False], "preload":[False]}
  run13
    input1 = '/mnt/shared/projects/virology/201609_BBSRC_Diagnostics/Data/
              synthetic/review-paper-test-datasets/SRR1269627Pear.fastq'
    kraken k-mers = 23 (m15), 27 (m15), 29 (m15)
    kraken_parameters = {"quick_minhits":[False], "preload":[False]}
  run14
    input1/2 = '/home/ae42909/data_forTesting/Plum_data/PlumVirus_SRR1291170'
    kraken k-mers = 31 (m15)
    kraken_parameters = {"quick_minhits":[False], "preload":[False]}
  run15
    input1/2 = '/home/ae42909/data_forTesting/Plum_data/PlumVirus_SRR1291171'
    kraken k-mers = 31 (m15)
    kraken_parameters = {"quick_minhits":[False], "preload":[False]}
  run16
    input1/2 = '/home/ae42909/data_forTesting/Plum_data/PlumClean_SRR1289655'
    kraken k-mers = 31 (m15)
    kraken_parameters = {"quick_minhits":[False], "preload":[False]}
  run17
    input1/2 = '/home/ae42909/data_forTesting/Plum_data/PlumClean_SRR1291169'
    kraken k-mers = 31 (m15)
    kraken_parameters = {"quick_minhits":[False], "preload":[False]}

kaiju/
  run1
    input1 = /home/ae42909/Scratch/100_Potato_withViruses_1.fastq
    input2 = /home/ae42909/Scratch/100_Potato_withViruses_2.fastq
  run2
    input1 = /home/ae42909/synthetic_RNAseq/mappingRNAseq/
              concatenated_fastaFiles/Potato_withViruses_1.fastq
    input2 = /home/ae42909/synthetic_RNAseq/mappingRNAseq/
              concatenated_fastaFiles/Potato_withViruses_2.fastq
  run3
    input1 = '/home/ae42909/Barerro_data/PB64-S1_clean.fq.gz_trim.fastq'
    kaiju_parameters = {"kaiju_minlen":[2, 3, 4],"kaiju_mismatch":[False, 1],
                        "kaiju_score":[45, 65, 85]}
  run4
    input1 = '/home/ae42909/Barerro_data/PB64-S2_clean.fq.gz_trim.fastq'
    kaiju_parameters = {"kaiju_minlen":[2, 3, 4],"kaiju_mismatch":[False, 1],
                        "kaiju_score":[45, 65, 85]}
  run5
    input1 = '/home/ae42909/Barerro_data/PB64-S3_clean.fq.gz_trim.fastq'
    kaiju_parameters = {"kaiju_minlen":[2, 3, 4],"kaiju_mismatch":[False, 1],
                        "kaiju_score":[45, 65, 85]}
  run6
    input1 = '/home/ae42909/Barerro_data/PB64-S4_clean.fq.gz_trim.fastq'
    kaiju_parameters = {"kaiju_minlen":[2, 3, 4],"kaiju_mismatch":[False, 1],
                        "kaiju_score":[45, 65, 85]}
  run7
    input1 = '/home/ae42909/Barerro_data/PB64-S5_clean.fq.gz_trim.fastq'
    kaiju_parameters = {"kaiju_minlen":[2, 3, 4],"kaiju_mismatch":[False, 1],
                        "kaiju_score":[45, 65, 85]}
  run8
    input1 = '/home/ae42909/Barerro_data/PB64-S7_clean.fq.gz_trim.fastq'
    kaiju_parameters = {"kaiju_minlen":[2, 3, 4],"kaiju_mismatch":[False, 1],
                        "kaiju_score":[45, 65, 85]}
  run9
    input1 = '/mnt/shared/projects/virology/201609_BBSRC_Diagnostics/Data/
            synthetic/review-paper-test-datasets/SRR1123893Pepper.fastq'
    kaiju_parameters = {"kaiju_minlen":[3, 7, 11, 15],
                        "kaiju_mismatch":[False, 1, 3, 5],
                        "kaiju_score":[45, 65, 85]}
  run10
    input1 = '/mnt/shared/projects/virology/201609_BBSRC_Diagnostics/Data/
            synthetic/review-paper-test-datasets/SRR1269627Pear.fastq'
    kaiju_parameters = {"kaiju_minlen":[3, 7, 11, 15],
                        "kaiju_mismatch":[False, 1, 3, 5],
                        "kaiju_score":[45, 65, 85]}
  run11
    input1/2 = '/home/ae42909/data_forTesting/Plum_data/PlumVirus_SRR1291170'
    kaiju_parameters = {"kaiju_minlen":[3, 7, 11, 15],
                        "kaiju_mismatch":[False, 1, 3, 5],
                        "kaiju_score":[45, 65, 85]}

  run12
    input1/2 = '/home/ae42909/data_forTesting/Plum_data/PlumVirus_SRR1291171'
    kaiju_parameters = {"kaiju_minlen":[3, 7, 11, 15],
                        "kaiju_mismatch":[False, 1, 3, 5],
                        "kaiju_score":[45, 65, 85]}
  run13
    input1/2 = '/home/ae42909/data_forTesting/Plum_data/PlumClean_SRR1289655'
    kaiju_parameters = {"kaiju_minlen":[3, 7, 11, 15],
                        "kaiju_mismatch":[False, 1, 3, 5],
                        "kaiju_score":[45, 65, 85]}
  run14
    input1/2 = '/home/ae42909/data_forTesting/Plum_data/PlumClean_SRR1291169'
    kaiju_parameters = {"kaiju_minlen":[3, 7, 11, 15],
                        "kaiju_mismatch":[False, 1, 3, 5],
                        "kaiju_score":[45, 65, 85]}

  run15 (same as 9 but separate files for paired data)
    input1 = '/mnt/shared/projects/virology/201609_BBSRC_Diagnostics/Data/
            synthetic/review-paper-test-datasets/SRR1123893Pepper'
    kaiju_parameters = {"kaiju_minlen":[3, 7, 11, 15],
                        "kaiju_mismatch":[False, 1, 3, 5],
                        "kaiju_score":[45, 65, 85]}
  run16
      input1 = /home/ae42909/data_forTesting/Apple_data/AppleVirus_SRR1089477.fastq
      kaiju_parameters = {'kaiju_score': [45, 65, 85], 'kaiju_minlen': [3, 7, 11, 15], 'kaiju_mismatch': [False, 1, 3, 5]}

run17
    input1 = /home/ae42909/data_forTesting/Apple_data/AppleClean_SRR1089478.fastq
    kaiju_parameters = {'kaiju_score': [45, 65, 85], 'kaiju_minlen': [3, 7, 11, 15], 'kaiju_mismatch': [False, 1, 3, 5]}

run18

run19

run20

run21

run22

run23

run24

run25

run26

run27

run28

run29
    input1 = /home/ae42909/data_forTesting/Goosefoot_data/GfCMV1_SRR503604.fastq
    kaiju_parameters = {'kaiju_score': [45, 65, 85], 'kaiju_minlen': [3, 7, 11, 15], 'kaiju_mismatch': [False, 1, 3, 5]}
run30
    input1 = /home/ae42909/data_forTesting/Goosefoot_data/GfCMV2_SRR503605.fastq
    kaiju_parameters = {'kaiju_score': [45, 65, 85], 'kaiju_minlen': [3, 7, 11, 15], 'kaiju_mismatch': [False, 1, 3, 5]}
run31
    input1 = /home/ae42909/data_forTesting/Goosefoot_data/GfClean_SRR503601.fastq
    kaiju_parameters = {'kaiju_score': [45, 65, 85], 'kaiju_minlen': [3, 7, 11, 15], 'kaiju_mismatch': [False, 1, 3, 5]}
run32
    input1 = /home/ae42909/data_forTesting/Goosefoot_data/GfTMV1_SRR503602.fastq
    kaiju_parameters = {'kaiju_score': [45, 65, 85], 'kaiju_minlen': [3, 7, 11, 15], 'kaiju_mismatch': [False, 1, 3, 5]}
run33
    input1 = /home/ae42909/data_forTesting/Goosefoot_data/GfTMV2_SRR503603.fastq
    kaiju_parameters = {'kaiju_score': [45, 65, 85], 'kaiju_minlen': [3, 7, 11, 15], 'kaiju_mismatch': [False, 1, 3, 5]}
run34
    input1 = /home/ae42909/data_forTesting/Algae_data/AlgaeVirus20_SRR924349.fastq
    kaiju_parameters = {'kaiju_score': [45, 65, 85], 'kaiju_minlen': [3, 7, 11, 15], 'kaiju_mismatch': [False, 1, 3, 5]}
run35
    input1 = /home/ae42909/data_forTesting/Algae_data/AlgaeVirus14_SRR924347.fastq
    kaiju_parameters = {'kaiju_score': [45, 65, 85], 'kaiju_minlen': [3, 7, 11, 15], 'kaiju_mismatch': [False, 1, 3, 5]}
run36
    input1 = /home/ae42909/data_forTesting/Algae_data/AlgaeVirus7_SRR924343.fastq
    kaiju_parameters = {'kaiju_score': [45, 65, 85], 'kaiju_minlen': [3, 7, 11, 15], 'kaiju_mismatch': [False, 1, 3, 5]}
run37
    input1 = /home/ae42909/data_forTesting/Algae_data/AlgaeClean_SRR924142.fastq
    kaiju_parameters = {'kaiju_score': [45, 65, 85], 'kaiju_minlen': [3, 7, 11, 15], 'kaiju_mismatch': [False, 1, 3, 5]}
run38
    input1 = /home/ae42909/data_forTesting/Algae_data/AlgaeVirus40_SRR924350.fastq
    kaiju_parameters = {'kaiju_score': [45, 65, 85], 'kaiju_minlen': [3, 7, 11, 15], 'kaiju_mismatch': [False, 1, 3, 5]}
run39
    input1 = /home/ae42909/data_forTesting/Algae_data/AlgaeVirus60_SRR924352.fastq
    kaiju_parameters = {'kaiju_score': [45, 65, 85], 'kaiju_minlen': [3, 7, 11, 15], 'kaiju_mismatch': [False, 1, 3, 5]}
run40
    input1 = /home/ae42909/data_forTesting/PotatoVD_data/PotatoVirus7_SRR3152175.fastq
    kaiju_parameters = {'kaiju_score': [45, 65, 85], 'kaiju_minlen': [2, 3, 4], 'kaiju_mismatch': [False, 1]}
run41
    input1 = /home/ae42909/data_forTesting/PotatoVD_data/PotatoVirus3_SRR3152188.fastq
    kaiju_parameters = {'kaiju_score': [45, 65, 85], 'kaiju_minlen': [2, 3, 4], 'kaiju_mismatch': [False, 1]}
run42
    input1 = /home/ae42909/data_forTesting/PotatoVD_data/PotatoClean5_SRR3152186.fastq
    kaiju_parameters = {'kaiju_score': [45, 65, 85], 'kaiju_minlen': [2, 3, 4], 'kaiju_mismatch': [False, 1]}
run43
    input1 = /home/ae42909/data_forTesting/PotatoVD_data/PotatoVirus10_SRR3152166.fastq
    kaiju_parameters = {'kaiju_score': [45, 65, 85], 'kaiju_minlen': [2, 3, 4], 'kaiju_mismatch': [False, 1]}
run44
    input1 = /home/ae42909/data_forTesting/PotatoVD_data/PotatoVirus1_SRR3152214.fastq
    kaiju_parameters = {'kaiju_score': [45, 65, 85], 'kaiju_minlen': [2, 3, 4], 'kaiju_mismatch': [False, 1]}
run45
    input1 = /home/ae42909/data_forTesting/PotatoVD_data/PotatoVirus9_SRR3152173.fastq
    kaiju_parameters = {'kaiju_score': [45, 65, 85], 'kaiju_minlen': [2, 3, 4], 'kaiju_mismatch': [False, 1]}
run46
    input1 = /home/ae42909/data_forTesting/PotatoVD_data/PotatoVirus5_SRR3152182.fastq
    kaiju_parameters = {'kaiju_score': [45, 65, 85], 'kaiju_minlen': [2, 3, 4], 'kaiju_mismatch': [False, 1]}
run47
    input1 = /home/ae42909/data_forTesting/PotatoVD_data/PotatoClean6_SRR3152187.fastq
    kaiju_parameters = {'kaiju_score': [45, 65, 85], 'kaiju_minlen': [2, 3, 4], 'kaiju_mismatch': [False, 1]}
run48
    input1 = /home/ae42909/data_forTesting/PotatoVD_data/PotatoClean2_SRR3152171.fastq
    kaiju_parameters = {'kaiju_score': [45, 65, 85], 'kaiju_minlen': [2, 3, 4], 'kaiju_mismatch': [False, 1]}
run49
    input1 = /home/ae42909/data_forTesting/PotatoVD_data/PotatoClean4_SRR3152212.fastq
    kaiju_parameters = {'kaiju_score': [45, 65, 85], 'kaiju_minlen': [2, 3, 4], 'kaiju_mismatch': [False, 1]}
run50
    input1 = /home/ae42909/data_forTesting/PotatoVD_data/PotatoVirus4_SRR3152183.fastq
    kaiju_parameters = {'kaiju_score': [45, 65, 85], 'kaiju_minlen': [2, 3, 4], 'kaiju_mismatch': [False, 1]}
run51
    input1 = /home/ae42909/data_forTesting/PotatoVD_data/PotatoVirus2_SRR3152189.fastq
    kaiju_parameters = {'kaiju_score': [45, 65, 85], 'kaiju_minlen': [2, 3, 4], 'kaiju_mismatch': [False, 1]}
run52
    input1 = /home/ae42909/data_forTesting/PotatoVD_data/PotatoClean1_SRR3152184.fastq
    kaiju_parameters = {'kaiju_score': [45, 65, 85], 'kaiju_minlen': [2, 3, 4], 'kaiju_mismatch': [False, 1]}
run53
    input1 = /home/ae42909/data_forTesting/PotatoVD_data/PotatoVirus8_SRR3152174.fastq
    kaiju_parameters = {'kaiju_score': [45, 65, 85], 'kaiju_minlen': [2, 3, 4], 'kaiju_mismatch': [False, 1]}
run54
    input1 = /home/ae42909/data_forTesting/PotatoVD_data/PotatoClean3_SRR3152213.fastq
    kaiju_parameters = {'kaiju_score': [45, 65, 85], 'kaiju_minlen': [2, 3, 4], 'kaiju_mismatch': [False, 1]}
run55
    input1 = /home/ae42909/data_forTesting/PotatoVD_data/PotatoVirus6_SRR3152180.fastq
    kaiju_parameters = {'kaiju_score': [45, 65, 85], 'kaiju_minlen': [2, 3, 4], 'kaiju_mismatch': [False, 1]}

# Parameter test analysis scripts and graphs in /home/ae42909/viral_diagnostics/
    parameter_testing/parameterTest_analysis/<tool>
The script run_results.py (one for each tool in their own directories, to
  run on cluster using run_results.sh) will take the results data
  ('<kraken/kaiju>_table_<tag>.txt' and '<kraken/kaiju>_labels_<tag>.txt') and
  create the virus_table_<tag>.txt (like   in the diagnostic tool).
  It then calculates the % of sequences for each virus
  from the total sequences identified as virus sequences at the species level.
  It then subsets for viruses that we expect to be in the sample and adds their
  percents together. The set of parameters with the highest percent for a sample
  will be chosen as the ideal parameters.
The script parameterResults_analysis.py (used in ipython) can be run for either tool
  (tool = "kraken" or "kaiju"), and you specify which runs you want to compare
  for that tool. The result is a table (specify name using "result_filename"
  variable) with each sample as a column and each set of parameters as the row
  index.
