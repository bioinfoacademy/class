#!/bin/bash
################################################################################
########### STAR - Fusion
## To Identify and annotate Fusion genes, we will be using STAR-Fusion.
## It is a component of the Trinity Cancer Transcriptome Analysis Toolkit (CTAT).
## STAR-Fusion (for identifying fusion events), Fusion Inspector
## (In-Silico validation of predicted Fusion transcripts) and Trinity (The de-novo assembler).
## STAR-Fusion also takes advantage of the STAR Alignment tool’s ability to determine
## splice junction information and uses  discordant read alignment information
## to predict gene fusion events.
################################################################################

## We will be using the demo data distributed with the STAR-Fusion package for
## this exercise. Go ahead and make a directory called Fusion on your Desktop
## and copy all the required files for the demo into that folder.

mkdir ~/Desktop/Fusion
cd ~/Desktop/Fusion
 
#link the example data files

ln -s ~/Desktop/rnaseqa/Fusion/rnaseq_1.fastq.gz
ln -s ~/Desktop/rnaseqa/Fusion/rnaseq_2.fastq.gz
#data coming from rnaseq snapshot
ln -s /home/manager/usr/local/mbcdata/rnaseqa/Fusion/ctat_genome_lib_build_dir

source /etc/profile.d/markcbm.sh
## Now we’re ready to predict fusion transcripts using STAR-Fusion as follows:
STAR-Fusion --left_fq rnaseq_1.fastq.gz --right_fq rnaseq_2.fastq.gz --genome_lib_dir ctat_genome_lib_build_dir

## The outputs are stored in a folder called ‘STAR-Fusion_outdir’. Let’s take a look

cd ~/Desktop/Fusion/STAR-Fusion_outdir
head star-fusion.fusion_predictions.abridged.tsv | column -t

#validate the fusions
cd ~/Desktop/Fusion
pip install igv-reports requests
STAR-Fusion --left_fq rnaseq_1.fastq.gz --right_fq rnaseq_2.fastq.gz --genome_lib_dir ctat_genome_lib_build_dir --FusionInspector validate

ls -ltr STAR-Fusion_outdir/FusionInspector-validate


#look at output, html result in;
#~/Desktop/Fusion/STAR-Fusion_outdir/FusionInspector-validate/finspector.fusion_inspector_web.html

## Free to use/modify/resuse for personal academic purposes
##© 2020 American Academy for Biomedical Informatics, Confidential and Proprietary 
