#!/bin/bash
#create a working directory
mkdir dapars

#Change directory to your dapars folder
cd dapars

#Link the test dataset to your dapars folder
ln -s ~/Desktop/rnaseqa/apa/* .

#make sure you have the necessary files
#Look at the contents of the wig files, UTR file, refseq id file, gene bed file.
ls

#Run DaPars
source /etc/profile.d/markcbm.sh
source activate dapars
python2 /usr/local/molbiocloud/dapars-0.9.1/src/DaPars_main.py DaPars_test_data_configure.txt 
conda deactivate

#visualize a significantly differentially expressed alternative polyadenylation usage using igv


## Free to use/modify/resuse for personal academic purposes
##© 2020 American Academy for Biomedical Informatics, Confidential and Proprietary 
