#!/bin/bash
## Create working directory and change to it
mkdir ~/Desktop/asarptest
cd ~/Desktop/asarptest

## Make symbolic links to input data
ln -s /usr/local/molbiocloud/ASARP/* .

## Check files...
ls

## Run test perl
perl testR.pl

## Create output directories
mkdir preproc
mkdir preproc/R1
mkdir preproc/R2
mkdir preproc/merged
mkdir demo.results

## Removing PCR Duplicates
ls demo/data/chr*.sam
export PERL5LIB="/usr/local/molbiocloud/ASARP"

perl rmDup.pl demo/data/chr1.R1.sam preproc/R1/chr1.rmDup.sam 1 >> rmDup.R1.log
perl rmDup.pl demo/data/chr5.R1.sam preproc/R1/chr5.rmDup.sam 1 >> rmDup.R1.log
perl rmDup.pl demo/data/chr10.R1.sam preproc/R1/chr10.rmDup.sam 1 >> rmDup.R1.log
perl rmDup.pl demo/data/chr1.R2.sam preproc/R2/chr1.rmDup.sam 1 >> rmDup.R2.log 
perl rmDup.pl demo/data/chr5.R2.sam preproc/R2/chr5.rmDup.sam 1 >> rmDup.R2.log
perl rmDup.pl demo/data/chr10.R2.sam preproc/R2/chr10.rmDup.sam 1 >> rmDup.R2.log

## Merge SAM Replicate Files
## First, create the replicate folders list
echo "preproc/R1" >> rep.demo.lst
echo "preproc/R2" >> rep.demo.lst

perl mergeSam.pl rep.demo.lst chr rmDup.sam preproc/merged > preproc/mergeSam.log

## Generate Allelic Reads
head -n2 demo/dna.snv.demo.lst
perl procReads.pl chr5 preproc/merged/chr5.sam demo/dna.snv.demo.lst preproc/chr5.snv preproc/chr5.bed 1 2
perl procReads.pl chr1 preproc/merged/chr1.sam demo/dna.snv.demo.lst preproc/chr1.snv preproc/chr1.bed 1 2
perl procReads.pl chr10 preproc/merged/chr10.sam demo/dna.snv.demo.lst preproc/chr10.snv preproc/chr10.bed 1 2

## Merge the SNV read counts and generate plots
perl mergeSnvs.pl preproc/ snv mono=0 preproc/rna.snv.demo.lst 1        

## Allele Specific Expression Analysis
perl asarp.pl demo.results/asarp.results demo/demo.config demo/demo.param


## Free to use/modify/resuse for personal academic purposes
##© 2020 American Academy for Biomedical Informatics, Confidential and Proprietary 
