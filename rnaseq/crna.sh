#!/bin/bash
## Create working directory and change to it
mkdir ~/Desktop/circ
cd ~/Desktop/circ

## link all the CIRCexplorer input files to your ‘circ’ folder
ln -s ~/Desktop/rnaseqa/circ/* .

## Look at the contents of the input files
head fusion_junction.bed
head ref.txt
head chr21.fa

## annotate known junctions and then identify circular fusion junctions
CIRCexplorer2 annotate -r ref.txt -g chr21.fa -b fusion_junction.bed

## Perform transcript assembly using cufflinks
CIRCexplorer2 assemble -r ref.txt -m tophat

## identify novel circRNA and perform alternative splicing of circular RNAs
CIRCexplorer2 denovo -r ref.txt -g chr21.fa -b fusion_junction.bed -d assemble --as as -m tophat -n paplus_tophat -o denovo

## Free to use/modify/resuse for personal academic purposes
##© 2020 American Academy for Biomedical Informatics, Confidential and Proprietary 
