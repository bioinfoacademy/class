#!/bin/bash
#Commands to run RNA-Seq analysis
#Numbers after the hashtag (#) sysmbol corresponds to sections in manual
#Authors: Vijay Nagarajan PhD, Harold Smith PhD

#####################################
########  SET-UP ##################

###Start a new MATE terminal to go to your home folder
ls 				##display contents
mkdir mydna   	##replace the FOLDERNAME with a name
cd mydna 		## go to your folder


##copy necessary files
ls ~/Desktop/dnaseq
ln -s ~/Desktop/dnaseq/human* .
ln -s ~/Desktop/dnaseq/dad_1.fq.gz .
ln -s ~/Desktop/dnaseq/dad_2.fq.gz .
ln -s ~/Desktop/dnaseq/dbSNP* .


############################################
###########################################

#############   PART I. QC & alignment  ######################

#------step 1: QC on fastq file-------

#test files available
#run fastqc on fastq files (generates .html file: with qc statistics) 
#should create dad_1_fastqc.html and dad_1_fastqc.zip
fastqc dad_1.fq.gz dad_2.fq.gz -o .  
firefox dad_1_fastqc.html  

#view results using x terminal (only works on Xterminal!!)

#Already done here due to time!!!!
#bwa index human.hg38.chr16.fa

#-------step 2: alignment------------

#1) using bwa, align fastq files, outputs stderr to log file, compresses the sam file
#for bwa-mem options: RG=read group: ID-id, LB-library, SM-Sample, PU-platform unit (ex.barcode), PL-platform (sequencer) 
bwa mem  ##shows options
bwa mem -R "@RG\\tID:dad\\tLB:dad\\tSM:dad\\tPU:FCC189PACXX\\tPL:ILLUMINA" -M human.hg38.chr16.fa dad_1.fq.gz dad_2.fq.gz | gzip > dad.sam.gz

#2)look at the content of the alignment file
less -S dad.sam.gz 

#3)convert SAM to sorted BAM- time consuming
picard SortSam INPUT=dad.sam.gz OUTPUT=dad.bam SORT_ORDER=coordinate  

#create index file
samtools index dad.bam 

#view the bam file, notice how it's aligned 
samtools view dad.bam | less -S 

#view the aligned reads 
samtools tview dad.bam human.hg38.chr16.fa
#press g, then type "chr16:400100"  
#q to quit

#--------now you have aligned reads----------------
#--------------------------------------------------

############   PART II. PRE-PROCESSING ALIGNED READS  ####################

#1) PICARD Mark Duplicates
picard MarkDuplicates INPUT=dad.bam OUTPUT=dad.chr16.dedup.bam ASSUME_SORTED=true M=dad.chr16.metrics.txt

#2) GATK Base recalibration
source /etc/profile.d/markcbm.sh

gatk BaseRecalibrator -I dad.bam -R human.hg38.chr16.fa --known-sites dbSNP151.chr16.vcf -O dad.recal_report_chr16.grp

#Recalibrate the aligned bam file
gatk ApplyBQSR -R human.hg38.chr16.fa -I dad.chr16.dedup.bam -O dad.chr16.recal.bam -bqsr dad.recal_report_chr16.grp


############   PART III. VARIANT DISCOVERY  ####################

# 1) GATK Haplotype caller to derive likelihoods of genotypes, creates gVCF
gatk HaplotypeCaller -R human.hg38.chr16.fa -I dad.chr16.recal.bam -O dad.chr16.g.vcf -ERC GVCF


# 2) GATK GenotypeGVCFs 

#single sample!! 
gatk GenotypeGVCFs -R human.hg38.chr16.fa -V dad.chr16.g.vcf -O dad.raw.vcf

#extract SNPs and apply hard filter, do the same for indels
# extracting SNPs
gatk SelectVariants -R human.hg38.chr16.fa -V dad.raw.vcf -select-type SNP -O dad_SNPs.vcf 

# applying hard filter
gatk VariantFiltration -R human.hg38.chr16.fa -V dad_SNPs.vcf -O dad.filtered_SNPs.vcf --filter-name "my_filter1" --filter-expression "DP < 3" --filter-name "my_filter2" --filter-expression "QD < 2.0" --filter-name "my_filter3" --filter-expression "MQ < 40.0" --filter-name "my_filter4" --filter-expression "FS > 60.0" --filter-name "my_filter5" --filter-expression "SOR > 3.0" --filter-name "my_filter6" --filter-expression "MQRankSum < -12.5" --filter-name "my_filter7" --filter-expression "ReadPosRankSum < -8.0"

# extracting INDELs
gatk SelectVariants -R human.hg38.chr16.fa -V dad.raw.vcf -select-type INDEL -O dad.indels.vcf

#applying hard filters to indels
gatk VariantFiltration -R human.hg38.chr16.fa -V dad.indels.vcf -O dad.filtered_indels.vcf --filter-name "indel_filter1" --filter-expression "QD < 2.0" --filter-name "indel_filter2" --filter-expression "ReadPosRankSum < -20.0" --filter-name "indel_filter3" --filter-expression "FS > 200.0" --filter-name "indel_filter4" --filter-expression "SOR > 10.0"

##merge two files

picard MergeVcfs I=dad.filtered_SNPs.vcf I= dad.filtered_indels.vcf O=dad.filtered_SNP_indels.vcf

#keep only selected variants
gatk SelectVariants -R human.hg38.chr16.fa -V dad.filtered_SNP_indels.vcf -O dad.final.vcf --exclude-filtered

## Free to use/modify/resuse for personal academic purposes
##© 2019 American Academy for Biomedical Informatics, Confidential and Proprietary 

