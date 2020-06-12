#!/bin/bash
#make working directory
mkdir ~/Desktop/viral

#change to the viral directory
cd ~/Desktop/viral

#link input files
ln -s /home/manager/Downloads/vigen/* .
rm output
mkdir output
ln -s /home/manager/Downloads/vigen/output/* output/.
rm -rf SRR1946637_un

#make sure all required files are in there
#look at the input files and their content
#raw data is from hepatocarcinoma https://www.ebi.ac.uk/ena/data/view/PRJNA279878
#raw data has been already aligned to the human reference genome and unaligned reads saved
ls
#unaligned reads are here
ls output/SRR1946637.temp/

#open the pipeline script viral.pipeline_public_final.sh in pluma and have a look at the workflow
#Run ViGen Search
bash viral.pipeline_public_final.sh

#look at the alignment stats
cat SRR1946637_un/log.bowtie2.txt

#the output file is *idxstats.txt in the output/samplename folder

## Free to use/modify/resuse for personal academic purposes
##© 2020 American Academy for Biomedical Informatics, Confidential and Proprietary 
