import pybedtools
import pandas as pd
# Tutorial https://daler.github.io/pybedtools/intersections.html
# Explanation of different arguments for intersectBED

# Get MACS data and TF data into bed-like format for pybedtools
def seqmonk_to_bedish (data_loc, outfile_name):
    seqmonk = pd.read_csv(data_loc, sep="\t", header=0)
    cols = seqmonk.columns.tolist()
    col_need = cols[1:4] + cols[0:1] + cols[5:8]
    bed = seqmonk[col_need]
    bed = bed[bed.Chromosome.str.contains("X" or "Y") == False] # Removing X and Y chromosome information as pybedtools.intersect trhows error
    bed.to_csv(outfile_name,sep='\t', header=False, index=False)

seqmonk_to_bedish(data_loc='/mnt/shared/scratch/ae42909/IL7Rko-paper/MACS/RAGko_rep1_MACSseqmonk.txt',
               outfile_name = '/mnt/shared/scratch/ae42909/IL7Rko-paper/peak_intersect/rag1_bed')

seqmonk_to_bedish(data_loc='/mnt/shared/scratch/ae42909/IL7Rko-paper/MACS/RAGko_rep2_MACSseqmonk.txt',
               outfile_name = '/mnt/shared/scratch/ae42909/IL7Rko-paper/peak_intersect/rag2_bed')

seqmonk_to_bedish(data_loc='/mnt/shared/scratch/ae42909/IL7Rko-paper/MACS/IL7RkoxRAGko_rep1_MACSseqmonk.txt',
               outfile_name = '/mnt/shared/scratch/ae42909/IL7Rko-paper/peak_intersect/il7ko1_bed')

seqmonk_to_bedish(data_loc='/mnt/shared/scratch/ae42909/IL7Rko-paper/MACS/IL7RkoxRAGko_rep2_MACSseqmonk.txt',
               outfile_name = '/mnt/shared/scratch/ae42909/IL7Rko-paper/peak_intersect/il7ko2_bed')

# Get ChIP-seq quantified data into the same format as above
def chip_to_bedish (data_loc, outfile_name):
    chip = pd.read_csv(data_loc, sep="\t", header= None, names=['Chromosome', 'Start', 'End','Name','5','6','7','8','9', '10'])
    bed = chip[chip.Chromosome.str.contains("X" or "Y") == False] # Removing X and Y chromosome information as pybedtools.intersect trhows error
    bed.to_csv(outfile_name,sep='\t', header=False, index=False)


chip_to_bedish(data_loc = '/mnt/shared/scratch/ae42909/IL7Rko-paper/E2A_peaks.narrowPeak',
               outfile_name = '/mnt/shared/scratch/ae42909/IL7Rko-paper/E2A_bed')

chip_to_bedish(data_loc = '/mnt/shared/scratch/ae42909/IL7Rko-paper/PAX5_peaks.narrowPeak',
               outfile_name = '/mnt/shared/scratch/ae42909/IL7Rko-paper/peak_intersect/PAX5_bed')

chip_to_bedish(data_loc = '/mnt/shared/scratch/ae42909/IL7Rko-paper/EBF1_peaks.narrowPeak',
               outfile_name = '/mnt/shared/scratch/ae42909/IL7Rko-paper/peak_intersect/EBF1_bed')

# Get intersecting peaks from each replicate and ones in both replicates. This function takes MACS peaks (identified from ATAC-seq data) from two replicates for each genotype and finds ChIP peaks from each TF (EBF1, PAX5 and E2A) that intersect with the MACS peaks. The result will be a list of ChIP peaks for each TF that seem to be present in the ATAC-seq data. However, 'intersecting peaks' could mean many things, and I haven't set a threshold as to how much they should intesect to be counted (if you want to change this see http://bedtools.readthedocs.io/en/latest/content/tools/intersect.html)
def intersecting_peaks (rep1_bed, rep2_bed, genotype, outdir):
    rep1 = pybedtools.BedTool(rep1_bed)
    rep2 = pybedtools.BedTool(rep2_bed)
    TF_list = ['EBF1', 'E2A', 'PAX5']
    peak_summary = pd.DataFrame(columns=['Rep1','Rep1_percent', 'Rep2', 'Rep2_percent', 'Rep1andRep2', 'Rep1andRep2_percent'], index=TF_list)

    for element in TF_list:
        TF =  pybedtools.BedTool('/mnt/shared/scratch/ae42909/IL7Rko-paper/peak_intersect/' + element + '_bed')
        outfile_rep1 = outdir + genotype + '_' + element + '_1'
        outfile_rep2 = outdir + genotype + '_' + element + '_2'
        outfile_joint = outdir + genotype + '_' + element + '_together'

       # TF_rep1 = TF.intersect(rep1) # get ChIP peaks that intersect with MACS peaks
        TF_rep1 = rep1.intersect(TF)  # get MACS peaks that intersect with ChIP peaks
        TF_rep1.saveas(outfile_rep1)
       # TF_rep2 = TF.intersect(rep2)
        TF_rep2 = rep2.intersect(TF)
        TF_rep2.saveas(outfile_rep2)
        TF_joint = TF_rep1.intersect(TF_rep2)
        TF_joint.saveas(outfile_joint)

        peak_summary.loc[element] = pd.Series({'Rep1': str(len(TF_rep1)),'Rep1_percent':str(round((float(len(TF_rep1))/float(len(TF)))*100, 2)), 'Rep2':str(len(TF_rep2)), 'Rep2_percent':str(round((float(len(TF_rep2))/float(len(TF)))*100, 2)), 'Rep1andRep2':str(len(TF_joint)), 'Rep1andRep2_percent':str(round((float(len(TF_joint))/float(len(TF)))*100, 2))})

    outfile_name = outdir + genotype + '_peakSummary'
    peak_summary.to_csv(outfile_name,sep='\t')



intersecting_peaks(rep1_bed = '/mnt/shared/scratch/ae42909/IL7Rko-paper/peak_intersect/rag1_bed',
                   rep2_bed = '/mnt/shared/scratch/ae42909/IL7Rko-paper/peak_intersect/rag2_bed',
                   genotype = 'rag',
                   outdir = '/mnt/shared/scratch/ae42909/IL7Rko-paper/peak_intersect/')

intersecting_peaks(rep1_bed = '/mnt/shared/scratch/ae42909/IL7Rko-paper/peak_intersect/il7ko1_bed',
                   rep2_bed = '/mnt/shared/scratch/ae42909/IL7Rko-paper/peak_intersect/il7ko2_bed',
                   genotype = 'il7ko',
                   outdir = '/mnt/shared/scratch/ae42909/IL7Rko-paper/peak_intersect/')

    

