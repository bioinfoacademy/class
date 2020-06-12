cd ~
mkdir sra
cd sra

#1. Then enter the following command:

esearch -db sra -q 'slr-seq AND mus_musculus[organism]' | efetch -format runinfo > slr_seq_runinfo.csv

#This will create a file named ‘slr_seq_runinfo.csv’. Use the unix command such as ‘less’ to look at the information present in this file. 

#2. Download the sra (SRR390728) run file for the cancer cell line RNA-Seq experiment SRX079566, using the below command;

prefetch -a "/usr/local/molbiocloud/aspera/cli/bin/ascp|/usr/local/molbiocloud/aspera/cli/etc/asperaweb_id_dsa.openssh" SRR390728

#-a option specifies the Aspera ascp command line utility and the required public key location.

#3. Run the sra-stat tools to look at the basic read statistics on the downloaded sra file

sra-stat --quick --xml SRR390728

# Extract data from sra file
#4. Run the fasterq-dump tool to download the fastq file, using the below command;

fasterq-dump SRR390728

#5. Dump SAM alignment from the sra file and write it in BAM format using the below command;

sam-dump SRR390728 | samtools view -bS - > SRR390728.bam

#6. Visualize the dumped bam file using the below command;
#Sort the bam file and create an index

samtools sort --threads 16 SRR390728.bam -o SRR390728.sorted.bam
samtools index SRR390728.sorted.bam

#use igv to view 1:89051831 SRR390728.sorted.bam Note reference is hg18

#Align reads using Magic-BLAST
#From NCBI Magic-BLAST documentation page:
#Magic-BLAST is a tool for mapping large next-generation RNA or DNA sequencing runs against a whole genome or transcriptome. Each alignment optimizes a composite score, taking into account simultaneously the two reads of a pair, and in case of RNA-seq, locating the candidate introns and adding up the score of all exons. This is very different from other versions of BLAST, where each exon is scored as a separate hit and read-pairing is ignored.

#We will align the reads from one of the covid 19 sequencing data to its reference genome
#In order to run magicblast, we will need a blastdb. Generate a blastdb for the reference genome as follows:

## create a directory in your home folder
cd ~/Desktop
mkdir covid
cd covid

#7. Fetch the reference genome sequence in FASTA format using Entrez Direct utility
efetch -db nuccore -id NC_045512.2 -format fasta > NC_045512.2.fasta

#8. Create a blastdb
makeblastdb -in NC_045512.2.fasta -dbtype nucl -parse_seqids -out covid -title 'covid'

#9. Run magicblast 
magicblast -sra SRR11652915 -db covid -num_threads 10 -out SRR11652915.aligns.sam -paired

#10. Look at the alignment
head SRR11652915.aligns.sam



