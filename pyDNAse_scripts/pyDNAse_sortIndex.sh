#!/bin/bash

#$ -cwd
#$ -j yes
#$ -pe smp 4
THREADS=4

#samtools sort -o IL7RxRAGko_2_sorted.bam lane4033_AGGCAGAA_IL7RxRAGko_2_L001_R1_val_1.fq.gz_bowtie2.bam
#samtools index IL7RxRAGko_2_sorted.bam  IL7RxRAGko_2_sorted.bam.bai

samtools sort -o IL7RxRAGko_1_sorted.bam lane4033_TCCTGAGC_IL7RxRAGko_1_L001_R1_val_1.fq.gz_bowtie2.bam
samtools index IL7RxRAGko_1_sorted.bam  IL7RxRAGko_1_sorted.bam.bai

samtools sort -o RAGko_1_sorted.bam lane4033_TAGGCATG_RAG2ko_2_L001_R1_val_1.fq.gz_bowtie2.bam
samtools index RAGko_1_sorted.bam  RAGko_1_sorted.bam.bai

samtools sort -o RAGko_2_sorted.bam lane4033_GGACTCCT_RAG2ko_females_L001_R1_val_1.fq.gz_bowtie2.bam
samtools index RAGko_2_sorted.bam  RAGko_2_sorted.bam.bai
