#!/bin/bash

################################################################################
########### Trinity
## This workflow will take us through the Trinity package and we can get a feel
## for using another tool to complete RNAseq work. Also we are using the same exact
## data so we can see how using a different workflow can change the results.
## In this instance we will also perform a de novo  transcriptome assembly
## as part of the exercise.
################################################################################

## Getting the data
source /etc/profile.d/markcbm.sh
python -m pip install numpy

cd ~/Desktop
mkdir Trinity
cd Trinity

## First we need to combine all the left and right reads into single files.
## For left reads
cat ~/Desktop/rnaseqa/tophat/*left.fq.gz > ALL.LEFT.fq.gz
## For right reads
cat ~/Desktop/rnaseqa/tophat/*right.fq.gz > ALL.RIGHT.fq.gz

## Run Trinity to assemble the transcriptome for these files
Trinity --seqType fq --SS_lib_type RF --left ALL.LEFT.fq.gz --right ALL.RIGHT.fq.gz --CPU 4 --max_memory 4G

## Let's take a look at the results from transcriptome assembly
head trinity_out_dir/Trinity.fasta

## Look at some metrics for transcriptome assembly
TrinityStats.pl trinity_out_dir/Trinity.fasta > Trinity_stats.txt

##open txt file or use:
less Trinity_stats.txt
#### How many genes did you assemble?
#### How many total transcripts?

##############################################################################
## Transcript alignment using gmap
##############################################################################
## build gmap index 
gmap_build -d genome -D . -k 13 ~/Desktop/rnaseqa/tophat/genome.fa

## align transcripts using gmap 
gmap -n 0 -D . -d genome trinity_out_dir/Trinity.fasta -f samse > trinity_gmap.sam

## convert sam file to sorted bam and create index 
samtools view -Su trinity_gmap.sam | samtools sort - -o trinity_gmap.sorted.bam

samtools index trinity_gmap.sorted.bam

## now load the bam file in igv, along with the transcripts assembled by stringtie
## to compare the transcripts from yesterday
igv -g ~/Desktop/hisat/genome.fa ~/Desktop/hisat/genes.bed trinity_gmap.sorted.bam ~/Desktop/hisat/stringtie_merged.gtf


##############################################################################
## Abundance Estimation - Trinity has some built-in pipelines for these workflows
##############################################################################
align_and_estimate_abundance.pl --prep_reference --seqType fq --est_method RSEM --aln_method bowtie --left ~/Desktop/rnaseqa/tophat/Sp_ds.left.fq.gz --right ~/Desktop/rnaseqa/tophat/Sp_ds.right.fq.gz --transcripts trinity_out_dir/Trinity.fasta --output_dir Sp_ds

### Process the other files immediately after using this nice cheat

!!:gs/ds/hs
!!:gs/hs/plat
!!:gs/plat/log

### Once finished, RSEM will have generated two files: ‘Sp_ds.isoforms.results’
### and ‘Sp_ds.genes.results’. These files contain the Trinity transcript and
### component (the Trinity analogs to Isoform and gene) rna-seq fragment counts
### and normalized expression values.
### Take a quick look at the ouput like this
head Sp_ds/RSEM.isoforms.results | column -t
head Sp_ds/RSEM.genes.results | column -t

### To run edgeR and identify differentially expressed transcripts, we need a
### data table containing the raw RNA-seq fragment counts for each transcript
### and sample analyzed.  We can combine the RSEM-computed isoform fragment counts
### into a matrix file.
abundance_estimates_to_matrix.pl --est_method RSEM Sp_ds/RSEM.isoforms.results Sp_hs/RSEM.isoforms.results Sp_log/RSEM.isoforms.results Sp_plat/RSEM.isoforms.results --gene_trans_map trinity_out_dir/Trinity.fasta.gene_trans_map --name_sample_by_basedir --out_prefix Trinity_trans

### You should find a matrix file called 'Trinity_trans.gene.counts.matrix', which
### contains the counts of RNA-Seq fragments mapped to each transcript.
### Examine the first few lines of the counts matrix:

head -n 20 Trinity_trans.gene.TMM.EXPR.matrix
head -n 20 Trinity_trans.isoform.TMM.EXPR.matrix

### To detect differentially expressed transcripts, run the Bioconductor package
### edgeR using our isoform counts matrix:
run_DE_analysis.pl --matrix Trinity_trans.isoform.counts.matrix --method edgeR --dispersion 0.1 --output isoforms_DE

### To detect differentially expressed genes, run the Bioconductor package
### edgeR using our gene counts matrix:
run_DE_analysis.pl --matrix Trinity_trans.gene.counts.matrix --method edgeR --dispersion 0.1 --output genes_DE

### The files '*.DE_results' contain the output from running EdgeR to identify
### differentially expressed transcripts in each of the pairwise sample comparisons.
### Trinity includes scripts for extracting transcripts that are above some
### statistical significance (FDR threshold) and fold-change in expression, and
### generating figures such as heatmaps and other useful plots,
cd ~/Desktop/Trinity/isoforms_DE

analyze_diff_expr.pl --matrix ../Trinity_trans.isoform.TMM.EXPR.matrix --output DE_isoforms

cd ~/Desktop/Trinity/genes_DE

analyze_diff_expr.pl --matrix ../Trinity_trans.gene.TMM.EXPR.matrix --output DE_genes

## Free to use/modify/resuse for personal academic purposes
##© 2020 American Academy for Biomedical Informatics, Confidential and Proprietary 

