## I highly recommend reading the Bioconductor Vignette for Ballgown available at this webpage:
## http://bioconductor.org/packages/release/bioc/vignettes/ballgown/inst/doc/ballgown.html

## load libraries
library(ballgown)
library(genefilter)
library(dplyr)
library(gplots)
library(ggfortify)

## set the data directory
data_directory = ('~/Desktop/hisat/ballgown')

## read outputs from stringtie and assign it to a variable
bg = ballgown(dataDir=data_directory, samplePattern='S', meas='all')

## let's take a peek at the data
bg

## three important slots in the bg object:
### structure -- genomic structure of the features (exons, introns and transcripts)
### expr -- expression data for the features
### indexes -- index info connecting the features 

## but first, let's add metadata information to the bg object
## NOTE: Ballgown doesn't perform DE analysis without replicates. 
## For convenience, we're combining Sp_ds & Sp_hs as replicates and Sp_plat & Sp_log as replicates.
pData(bg) = data.frame(id=sampleNames(bg), group=factor(c(1,1,2,2)))

## let's see what's inside the structure slot
names(structure(bg))

## now, let's look at the individual features within the structure slot
structure(bg)$exon
structure(bg)$intron
structure(bg)$trans

## You can isolate FPKM values for each transcript like this:
transcript_fpkm = texpr(bg, 'FPKM')
head(transcript_fpkm)

# You could also isolate coverage information for each transcript like this: 
# Note: changing the command to eexpr or iexpr will give you similar information for exons and introns.
transcript_cov = texpr(bg, 'cov')
head(transcript_cov)

## Extracting all available stats for each transcript
whole_tx_table = texpr(bg, 'all')
head(whole_tx_table)

## Isolate multi-map-corrected average per-base read coverage for each exon
exon_mcov = eexpr(bg, 'mcov')
head(exon_mcov)

## Isolate read count for each intron
junction_rcount = iexpr(bg)
head(junction_rcount)

## Isolate all available stats for each intron
whole_intron_table = iexpr(bg, 'all')
head(whole_intron_table)

## Isolate gene-level expression
gene_expression = gexpr(bg)
head(gene_expression)

## generate a boxplot of transcript expression by FPKM for each of the samples
transcript_fpkm_log2 = log2(transcript_fpkm+1)
head(transcript_fpkm_log2)
par(mar=c(10,4,4,2))
boxplot(transcript_fpkm_log2,col=as.numeric(pData(bg)$group),
        las=2,ylab='log2(FPKM+1)')
dev.off()

## Plot a heatmap of transcript FPKM across all samples
cols = colorpanel(99,"blue","black","yellow")
heatmap.2(transcript_fpkm_log2,density.info = "none",trace ="none", margins = c(10,15), cexCol = 1.2,keysize = 1,col = cols,scale = "none")
dev.off()

## PCA
pca = prcomp(t(transcript_fpkm_log2))
autoplot(pca,data=pData(bg),colour="group",label=T,label.label="id")
dev.off()

## Extract indexes from ballgown object
names(indexes(bg))

## mapping exons to transcripts
head(indexes(bg)$e2t)

## mapping introns to transcripts
head(indexes(bg)$i2t)

## mapping transcripts to genes
head(indexes(bg)$t2g)

## Mapping exons to transcripts and transcripts to gene.
exon_transcript_table = indexes(bg)$e2t
transcript_gene_table = indexes(bg)$t2g
head(transcript_gene_table)

## Performing statistical testing for differential expression by transcript & gene and writing the results to a table.
  results_transcripts = stattest(bg, feature='transcript', meas='FPKM',covariate='group', getFC=TRUE)
  head(results_transcripts)
  write.table(results_transcripts,
              file="ballgown_transcript_results.txt", 
              sep="\t",quote=F,row.names = F)

results_genes = stattest(bg, feature='gene', meas='FPKM',covariate='group', getFC=TRUE)
head(results_genes)
write.table(results_genes,file="ballgown_gene_results.txt", sep="\t",quote=F,row.names = F)

## Add geneIDs to the differentially expressed transcripts
results_transcripts = data.frame(geneIDs=ballgown::geneIDs(bg), results_transcripts)
head(results_transcripts)
results_transcripts = arrange(results_transcripts,pval)
results_genes = arrange(results_genes,pval)

## Plot a map of this transcript across all the samples and color by transcript
head(results_transcripts)
plotTranscripts(gene='MSTRG.135', gown=bg,samples = sampleNames(bg),  meas='FPKM', colorby='transcript', main='transcripts from gene MSTRG.135:  FPKM')
dev.off()

## boxplot showing log2 FPKM of most differentially expressed transcript
idx = which(geneIDs(bg)=="MSTRG.135")
idx
plot(transcript_fpkm_log2[idx,] ~ pData(bg)$group, border=c(1,2),
     main=paste(ballgown::geneIDs(bg)[idx]),pch=19,xlab="Group",
     ylab='log2(FPKM+1)')
points(transcript_fpkm_log2[idx,] ~ jitter(as.numeric(pData(bg)$group)),
       col=as.numeric(pData(bg)$group),pch=20)
dev.off()
