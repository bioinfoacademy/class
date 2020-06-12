#!/bin/bash
#Author: Vijay Nagarajan PhD
#This script takes the IRFinder suite trimmed and aligned data for ERP001458 - an RNA-Seq study of HNRNPC knockdown in HeLa cells.
#raw fastq files downloaded from https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-1147/samples/?full=true

#source the program location
source /etc/profile.d/markcbm.sh

#make a working directory
mkdir ~/Desktop/irfinder

#change to your irfinder directory
cd ~/Desktop/irfinder

#make directories for storing kd and control data
mkdir -p kd/kd1
mkdir -p kd/kd2
mkdir -p control/control1
mkdir -p control/control2

echo "Directories setup complete"

#link the IRFinder reference
ln -s ~/Desktop/rnaseqa/ir/Human-hg19-release75 .

#link to the unsorted alignment data for kd replicates and control replicates
ln -s ~/Desktop/rnaseqa/ir/kd/kd1/Unsorted.bam kd/kd1/.
ln -s ~/Desktop/rnaseqa/ir/kd/kd2/Unsorted.bam kd/kd2/.
ln -s ~/Desktop/rnaseqa/ir/control/control1/Unsorted.bam control/control1/.
ln -s ~/Desktop/rnaseqa/ir/control/control2/Unsorted.bam control/control2/.

#link to the unsorted merged replicated alignment data
ln -s ~/Desktop/rnaseqa/ir/kd/kdmergedunsorted.bam kd/.
ln -s ~/Desktop/rnaseqa/ir/control/controlmergedunsorted.bam control/.
echo "Input example RNA-Seq alignment bam files linked"

#Quantify Intron retention using aligned reads for the kd data
IRFinder -m BAM -d kd/kd1 -r Human-hg19-release75 kd/kd1/Unsorted.bam
rm kd/kd1/unsorted.frag.bam
echo "Intron Retention quantified for kd1"

IRFinder -m BAM -d kd/kd2 -r Human-hg19-release75 kd/kd2/Unsorted.bam
rm kd/kd2/unsorted.frag.bam
echo "Intron Retention quantified for kd2"

IRFinder -m BAM -d kd/pooled -r Human-hg19-release75 kd/kdmergedunsorted.bam
rm  kd/pooled/unsorted.frag.bam
echo "Intron Retention quantified for merged kd"

#Quantify Intron retention using aligned reads for the control data
IRFinder -m BAM -d control/control1 -r Human-hg19-release75 control/control1/Unsorted.bam
rm control/control1/unsorted.frag.bam
echo "Intron Retention quantified for control1"

IRFinder -m BAM -d control/control2 -r Human-hg19-release75 control/control2/Unsorted.bam
rm control/control2/unsorted.frag.bam
echo "Intron Retention quantified for control2"

IRFinder -m BAM -d control/pooled -r Human-hg19-release75 control/controlmergedunsorted.bam
rm control/pooled/unsorted.frag.bam
echo "Intron Retention quantified for merged control"

#Run differential Intron retention analysis between kd and control data
analysisWithLowReplicates.pl -A kd/pooled/IRFinder-IR-nondir.txt kd/kd1/IRFinder-IR-nondir.txt kd/kd2/IRFinder-IR-nondir.txt -B control/pooled/IRFinder-IR-nondir.txt control/control1/IRFinder-IR-nondir.txt control/control2/IRFinder-IR-nondir.txt > kd-v-control.tab

echo "Differential intron retention analysis between kd and control complete. Check the results in kd-v-control.tab file"

## Free to use/modify/resuse for personal academic purposes
##© 2020 American Academy for Biomedical Informatics, Confidential and Proprietary 


