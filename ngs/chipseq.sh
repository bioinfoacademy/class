#!/bin/bash
#Commands to run ChIP-Seq analysis
#Numbers after the hashtag (#) sysmbol corresponds to sections in manual

#4
##1

###1.1. change to your desktop, create macs2 directory, copy input files
cd Desktop
mkdir macs2
cp chipseq/* macs2

###1.2. change to macs2 working directory
cd macs2

###look at input files
head Input_tags.bed 
head Treatment_tags.bed
wc -l Input_tags.bed
wc -l Treatment_tags.bed

###1.3. Run macs2 program
macs2 callpeak -t Treatment_tags.bed -c Input_tags.bed -f BED -g hs -n FoxA1 -B -p 0.005

##2
### list to see the results
ls

### peaks
head FoxA1_peaks.narrowPeak 

### peak summits
head FoxA1_summits.bed 

###
head FoxA1_treat_pileup.bdg

###instructor alone run this
###less FoxA1_treat_pileup.bdg
###q to quit if u did the above less command

### total number of peaks identified
wc -l FoxA1_peaks.narrowPeak

###grep to find how many peaks in chr1
grep -w 'chr1' FoxA1_summits.bed | wc -l

###rename peaks to bed file
cp FoxA1_peaks.narrowPeak FoxA1_peaks.bed

###run enrichr with foxa1 peak bed file, hg38, 2000 genes

###2.1. Extract top 1000 peaks 
head -n 1000 FoxA1_summits.bed > 1000summits.bed

###head 4 lines of new summit file
head -n4 FoxA1_summits.bed

###2.2. Convert 5 column macs bed to ucsc 4 column bed
cut -f 1-4 1000summits.bed > 1000summits4.bed

###head 4 lines of 4 column bed
head -n4 1000summits4.bed

###2.2. Extract FASTA format sequences from uCSC
###google for ucsc table browser
###click add custom tracks button
###use browse button to choose 1000summits4.bed file and click submit
###choose view in Table Browser option and click go
###in table browser choose output format - sequence, name output as 'sequence.fasta'
###click getoutput, add 10 bp up and 10 bp down and click get sequence button
###save the file (it gets saved in Downloads folder)
###move sequence.fasta to your macs2 working directory
cp ~/Downloads/sequence.fasta .
head sequence.fasta

###2.3. Run meme-chip for two motifs, zoops options
source /etc/profile.d/markcbm.sh
meme-chip -meme-mod zoops -meme-nmotifs 2 sequence.fasta

## Free to use/modify/resuse for personal academic purposes
##© 2020 American Academy for Biomedical Informatics, Confidential and Proprietary 
