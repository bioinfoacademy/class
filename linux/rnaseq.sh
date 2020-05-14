#!/bin/bash
#RNA-Seq analysis pipeline
#Numbers after the hashtag (#) sysmbol corresponds to sections in manual

#14 Data Preparation
#this source file tells the shell where our programs are
source /etc/profile.d/markcbm.sh

#Copy example data files
cd ~/Desktop/linux/advanced/rnaseq/

##14.1. Step1: Run FASTQC for our input files
fastqc fastq/*.fastq

##14.2. Step2: Build index - already created for you
cd index
bowtie-build mm9_chr1.fa mm9_chr1
cd ..

##14.3. Step 3:  Align fastq
#Run tophat2 with the following command
tophat2 -G mm9_chr1.gtf -o  tophat_wt/  index/mm9_chr1 fastq/myoblast_wt.fastq
tophat2 -G mm9_chr1.gtf -o  tophat_del/  index/mm9_chr1 fastq/myoblast_del.fastq

##14.4. Step4 : Check if the tophat_wt and tophat_del are found
#You should find accepted_hits.bam file which is the alignment file
ls -lh tophat_del
ls -lh tophat_wt

#Look at the alignment summary
cat tophat_del/align_summary.txt
cat tophat_wt/align_summary.txt

##14.5. Step6: Index aligned data
samtools19 index tophat_wt/accepted_hits.bam 
samtools19 index tophat_del/accepted_hits.bam 

##14.6. Step 7: Differential Expression
#Run CuffDiff with the following command, to calculate the differential expression between the wt and the del samples;
cuffdiff --no-update-check -o cuffdiff_out -L wt,del mm9_chr1.gtf tophat_wt/accepted_hits.bam tophat_del/accepted_hits.bam
ls -lh cuffdiff_out

## Free to use/modify/resuse for personal academic purposes
##© 2020 American Academy for Biomedical Informatics, Confidential and Proprietary 
