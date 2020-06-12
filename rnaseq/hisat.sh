#!/bin/bash
#################################
## FASTQC
#################################
# create working directory and make symlinks to the required files

mkdir ~/Desktop/FastQC
cd ~/Desktop/FastQC

source /etc/profile.d/markcbm.sh
ln -s ~/Desktop/rnaseqa/variant/Tumor_RNAseq_R1.fastq
ln -s ~/Desktop/rnaseqa/variant/Tumor_RNAseq_R2.fastq

# run fastqc
fastqc Tumor_RNAseq_R1.fastq -o ~/Desktop/FastQC
fastqc Tumor_RNAseq_R2.fastq -o ~/Desktop/FastQC

#################################
## TRIMMOMATIC
#################################
# create working directory and make symlinks to the required files
mkdir ~/Desktop/Trimming
cd ~/Desktop/Trimming
 
ln -s ~/Desktop/rnaseqa/variant/Tumor_RNAseq_R1.fastq .
ln -s ~/Desktop/rnaseqa/variant/Tumor_RNAseq_R2.fastq .

# run trimmomatic
TRIMMOMATIC=/usr/local/molbiocloud/Trimmomatic-0.39/trimmomatic-0.39.jar
TRUSEQ2_FA=/usr/local/molbiocloud/Trimmomatic-0.39/adapters/TruSeq2-PE.fa

java -jar ${TRIMMOMATIC} PE \
  -threads 4 \
  Tumor_RNAseq_R1.fastq Tumor_RNAseq_R2.fastq \
  Tumor_RNAseq_R1_paired.fq Tumor_RNAseq_R1_unpaired.fq \
  Tumor_RNAseq_R2_paired.fq Tumor_RNAseq_R2_unpaired.fq \
  ILLUMINACLIP:${TRUSEQ2_FA}:2:30:10 \
  LEADING:3 \
  TRAILING:3 \
  SLIDINGWINDOW:4:20 \
  MINLEN:70

#################################
## HISAT2
#################################
# create working directory and make symlinks to the required files
mkdir ~/Desktop/hisat
cd ~/Desktop/hisat

ln -s ~/Desktop/rnaseqa/tophat/* .

# generate hisat index
mkdir hisat_index
hisat2-build genome.fa hisat_index/SpIndex

# run hisat2
hisat2 \
 -p 8 \
 --rna-strandness RF \
 --dta \
 -x hisat_index/SpIndex \
 -1 Sp_ds.left.fq.gz \
 -2 Sp_ds.right.fq.gz \
 -S Sp_ds.sam \
 --summary-file Sp_ds_alignStats.txt

!!:gs/ds/hs
!!:gs/hs/plat
!!:gs/plat/log

# move sam files  
mkdir Aligned_sam
mv *.sam Aligned_sam/

# move alignment stats files
mkdir AlignedStats
mv *.txt AlignedStats

# convert sam to bam and create index files
mkdir Aligned_bam 
samtools view -Su Aligned_sam/Sp_ds.sam \
  | samtools sort - -o Aligned_bam/Sp_ds.sorted.bam 

!!:gs/ds/hs
!!:gs/hs/plat
!!:gs/plat/log

cd ~/Desktop/hisat/Aligned_bam/
for f in *.bam ; do 
   samtools index $f ; 
done 

# launch igv 
# DO THIS IN A NEW TERMINAL
source /etc/profile.d/markcbm.sh
igv -g ~/Desktop/hisat/genome.fa \
   ~/Desktop/hisat/genes.bed \
   ~/Desktop/hisat/Aligned_bam/*.bam

# RETURN TO THE OLD TERMINAL, by closing IGV
# run stringtie
cd ..
stringtie \
  Aligned_bam/Sp_ds.sorted.bam \
  -p 8 \
  --rf \
  -o Assembled_transcripts/Sp_ds.gtf 

!!:gs/ds/hs
!!:gs/hs/plat
!!:gs/plat/log

# merge stringtie output
ls Assembled_transcripts/*.gtf > mergelist.txt
cat mergelist.txt
stringtie --merge -p 8 -o stringtie_merged.gtf mergelist.txt

# Now add the `stringtie_merged.gtf` file to the open igv window
# to view the transcripts assembled by stringtie

# generate ballgown files 
stringtie -e -B -p 8 \
  -G stringtie_merged.gtf \
  -o ballgown/Sp_ds/Sp_ds.gtf \
  Aligned_bam/Sp_ds.sorted.bam

!!:gs/ds/hs
!!:gs/hs/plat
!!:gs/plat/log

#################################
## BALLGOWN
#################################
# Now launch RStudio and open the ~/Desktop/Scripts/ballgown.r file

## Free to use/modify/resuse for personal academic purposes
##© 2020 American Academy for Biomedical Informatics, Confidential and Proprietary 
