#!/bin/bash
#Commands to run RNA-Seq analysis
#Numbers after the hashtag (#) sysmbol corresponds to sections in manual

#2
##2.1

###1
cd Desktop
mkdir Trimming
cd Trimming
ln -s ~/Desktop/rnaseq/input/Trimming/* .

###2
mkdir ~/Desktop/fastqc
source /etc/profile.d/markcbm.sh
fastqc Tumor_RNAseq_R1.fastq -o ~/Desktop/fastqc
fastqc Tumor_RNAseq_R2.fastq -o ~/Desktop/fastqc

###3
cd ~/Desktop/fastqc
firefox Tumor_RNAseq_R1_fastqc.html

#If firefox doesnt open through command line, navigate to the ~/Desktop/fastqc results folder and open the results html file

###4
cd ~/Desktop/Trimming
java -jar /usr/local/molbiocloud/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 Tumor_RNAseq_R1.fastq Tumor_RNAseq_R2.fastq Tumor_RNAseq_R1_pe.fq Tumor_RNAseq_R1_se.fq Tumor_RNAseq_R2_pe.fq Tumor_RNAseq_R2_se.fq ILLUMINACLIP:adapt.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:70

###5
fastqc Tumor_RNAseq_R1_pe.fq -o ~/Desktop/fastqc

###6
cd ~/Desktop/fastqc
firefox Tumor_RNAseq_R1_pe_fastqc.html

##2.2

###1
cd ~/Desktop
mkdir hisat
cd hisat
ln -s ~/Desktop/rnaseq/input/hisat/* .

###2
mkdir fastqc 
fastqc *.fastq -o fastqc/
cd fastqc 
firefox rna_cntl_1_fastqc.html
cd ..

###3
mkdir hisat_index
hisat2-build chr16.fa chr16_Index
mv chr16_Index* hisat_index/

###4
hisat2 -p 8 --dta -x hisat_index/chr16_Index -U rna_cntl_1.fastq -S rna_cntl_1.sam --summary-file rna_cntl_1_alignStats.txt

###5
!!:gs/cntl_1/cntl_2
!!:gs/cntl_2/expt_1
!!:gs/expt_1/expt_2

###6
mkdir Sam_aligned
mv *.sam Sam_aligned/

###7
mkdir AlignedStats
mv *.txt AlignedStats/

###8
mkdir Bam_aligned 
samtools view -u Sam_aligned/rna_cntl_1.sam | samtools sort - -o Bam_aligned/rna_cntl_1.sorted.bam

###9
!!:gs/cntl_1/cntl_2
!!:gs/cntl_2/expt_1
!!:gs/expt_1/expt_2

###10
stringtie Bam_aligned/rna_cntl_1.sorted.bam -p 8 -o Assembled_transcripts/rna_cntl_1.gtf 
!!:gs/cntl_1/cntl_2
!!:gs/cntl_2/expt_1
!!:gs/expt_1/expt_2

###11
ls Assembled_transcripts/*.gtf > Assembled_transcripts/mergelist.txt
cat Assembled_transcripts/mergelist.txt
stringtie --merge -p 8 -o Assembled_transcripts/stringtie_merged.gtf Assembled_transcripts/mergelist.txt

###12
stringtie -e -B -p 8 -G Assembled_transcripts/stringtie_merged.gtf -o ballgown/rna_cntl_1/rna_cntl_1.gtf Bam_aligned/rna_cntl_1.sorted.bam

###13
!!:gs/cntl_1/cntl_2
!!:gs/cntl_2/expt_1
!!:gs/expt_1/expt_2

## Free to use/modify/resuse for personal academic purposes
##© 2020 American Academy for Biomedical Informatics, Confidential and Proprietary 
