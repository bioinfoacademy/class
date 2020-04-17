#!/bin/bash
#Commands to run RNA-Seq analysis
#Numbers after the hashtag (#) sysmbol corresponds to sections in manual

########  1. SET-UP ##

### Start a new MATE terminal to go to your home folder

# 1.1. Display working directory contents
ls

# Create a new directory named 'mydna'
mkdir mydna

# Change to the 'mydna' directory
cd mydna

## 1.2. copy necessary files
ls ~/Desktop/dnaseq
ln -s ~/Desktop/dnaseq/human* .
ln -s ~/Desktop/dnaseq/dad_1.fq.gz .
ln -s ~/Desktop/dnaseq/dad_2.fq.gz .
ln -s ~/Desktop/dnaseq/dbSNP* .

############# 2. QC & alignment ###

# 2.1. QC on fastq file

#test files available
#run fastqc on fastq files (generates .html file: with qc statistics) 
#should create dad_1_fastqc.html and dad_1_fastqc.zip
fastqc dad_1.fq.gz dad_2.fq.gz -o .  
firefox dad_1_fastqc.html  

#view results using x terminal (only works on Xterminal!!)

# 2.2. Index reference genome
#Already done here due to time!!!!
#bwa index human.hg38.chr16.fa

# 2.3. alignment

#1) using bwa, align fastq files, outputs stderr to log file, compresses the sam file
#for bwa-mem options: RG=read group: ID-id, LB-library, SM-Sample, PU-platform unit (ex.barcode), PL-platform (sequencer) 

#show options
bwa mem

#run alignment
bwa mem -R "@RG\\tID:dad\\tLB:dad\\tSM:dad\\tPU:FCC189PACXX\\tPL:ILLUMINA" -M human.hg38.chr16.fa dad_1.fq.gz dad_2.fq.gz | gzip > dad.sam.gz

#2)look at the content of the alignment file
less -S dad.sam.gz 

#3)convert SAM to sorted BAM- time consuming
picard SortSam INPUT=dad.sam.gz OUTPUT=dad.bam SORT_ORDER=coordinate  

#4)create index file
samtools index dad.bam 

#5)view the bam file, notice how it's aligned 
samtools view dad.bam | less -S 

#6)view the aligned reads using igv
source /etc/profile.d/markcbm.sh

#start igv from terminal, by running the command 'igv'
igv
#load alignment file dad.bam, load genome from file human.hg38.chr16.fa
#search for "chr16:400100"  

#--------now you have aligned reads----------------
#--------------------------------------------------

############ 3. PRE-PROCESSING ALIGNED READS  ##

# 3.1. PICARD Mark Duplicates
picard MarkDuplicates INPUT=dad.bam OUTPUT=dad.chr16.dedup.bam ASSUME_SORTED=true M=dad.chr16.metrics.txt

# 3.2. GATK Base recalibration
source /etc/profile.d/markcbm.sh
export PATH=/usr/lib/jvm/java-1.8.0-openjdk-amd64/bin/:$PATH

gatk BaseRecalibrator -I dad.bam -R human.hg38.chr16.fa --known-sites dbSNP151.chr16.vcf -O dad.recal_report_chr16.grp

# 3.3. Recalibrate the aligned bam file
gatk ApplyBQSR -R human.hg38.chr16.fa -I dad.chr16.dedup.bam -O dad.chr16.recal.bam -bqsr dad.recal_report_chr16.grp

############ 4. VARIANT DISCOVERY  ##

# 4.1. GATK Haplotype caller to derive likelihoods of genotypes, creates gVCF
gatk HaplotypeCaller -R human.hg38.chr16.fa -I dad.chr16.recal.bam -O dad.chr16.g.vcf -ERC GVCF

# 4.2. GATK GenotypeGVCFs 

#single sample!! 
gatk GenotypeGVCFs -R human.hg38.chr16.fa -V dad.chr16.g.vcf -O dad.raw.vcf

# 4.3. extract SNPs and apply hard filter, do the same for indels
# a) extracting SNPs
gatk SelectVariants -R human.hg38.chr16.fa -V dad.raw.vcf -select-type SNP -O dad_SNPs.vcf 

# b) applying hard filter
gatk VariantFiltration -R human.hg38.chr16.fa -V dad_SNPs.vcf -O dad.filtered_SNPs.vcf --filter-name "my_filter1" --filter-expression "DP < 3" --filter-name "my_filter2" --filter-expression "QD < 2.0" --filter-name "my_filter3" --filter-expression "MQ < 40.0" --filter-name "my_filter4" --filter-expression "FS > 60.0" --filter-name "my_filter5" --filter-expression "SOR > 3.0" --filter-name "my_filter6" --filter-expression "MQRankSum < -12.5" --filter-name "my_filter7" --filter-expression "ReadPosRankSum < -8.0"

# c) extracting INDELs
gatk SelectVariants -R human.hg38.chr16.fa -V dad.raw.vcf -select-type INDEL -O dad.indels.vcf

# d) applying hard filters to indels
gatk VariantFiltration -R human.hg38.chr16.fa -V dad.indels.vcf -O dad.filtered_indels.vcf --filter-name "indel_filter1" --filter-expression "QD < 2.0" --filter-name "indel_filter2" --filter-expression "ReadPosRankSum < -20.0" --filter-name "indel_filter3" --filter-expression "FS > 200.0" --filter-name "indel_filter4" --filter-expression "SOR > 10.0"

# 4.4. merge two files

picard MergeVcfs I=dad.filtered_SNPs.vcf I= dad.filtered_indels.vcf O=dad.filtered_SNP_indels.vcf

# 4.5. keep only selected variants
gatk SelectVariants -R human.hg38.chr16.fa -V dad.filtered_SNP_indels.vcf -O dad.final.vcf --exclude-filtered

## Free to use/modify/resuse for personal academic purposes
##© 2020 American Academy for Biomedical Informatics, Confidential and Proprietary 

