#!/usr/bin/env Rscript

## author: Emily Tsang

## order transcripts by usage within each gene
## and gather their tpm info

## The environment variables should be set globally, e.g. in bashrc ##
rnadir = Sys.getenv('GTEX_RNAv7')
dir = Sys.getenv('EGTEXDIR')
v6dir = Sys.getenv('GTEXv6')
######################################################################

library(data.table)
library(reshape2)
library(stringr)
library(ggplot2)
library(dplyr)

## NOTE: The expression files were updated since this was run.
## The output used downstream is from these outdated files.
## Since the newer files just have a slightly different set of individuals, it should only have a minor effect on results.
## If you want to runt hings with the newer files, uncomment the two lines below and comment out the two after that.
#percfile = paste0(rnadir,'/GTEx_Analysis_2016-01-15_v7_RSEMv1.2.22_transcript_isopct.txt.gz')
#tpmfile = paste0(rnadir, '/GTEx_Analysis_2016-01-15_v7_RSEMv1.2.22_transcript_tpm.txt.gz')
percfile = paste0(rnadir,'/GTEx_Analysis_2016-09-07_RSEMv1.2.22_transcript_isopct.txt.gz')
tpmfile = paste0(rnadir, '/GTEx_Analysis_2016-09-07_RSEMv1.2.22_transcript_tpm.txt.gz')
gencodefile = paste0(v6dir, '/reference_files/gencode.v19.genes.v6p.patched_contigs.gtf.gz')
mapfile = paste0(dir, '/sampleInfo/master.platemap.txt')


## get input and output files passed from the command line
args = commandArgs(trailingOnly = TRUE)
suffix = args[1]
suffix = ifelse(is.na(suffix), '', suffix)
genefile = paste0(dir, '/geneSelection/selected.genes', suffix, '.txt')
outfile = paste0(dir, '/selectedSites/transcripts', suffix, '.txt')

cat(outfile, genefile, '\n')
stop('')

### FUNCTIONS ##############

## function to format tissue names without spaces
## to make them match the format used in other gtex files
format.tissue.names = function(tissue){
    tissue = gsub(" - ", "_", tissue, fixed = T)
    tissue = gsub("(", "", tissue, fixed = T)
    tissue = gsub(")", "", tissue, fixed = T)
    tissue = gsub(" ", "_", tissue, fixed = T)
    return(tissue)
}

## function to fix treanscript file column names of isoform percentages
fix.indids = function(indid) {
    indid.split = strsplit(indid, '-')[[1]]
    return(paste(indid.split[1:(length(indid.split)-2)], collapse = '-'))
}

## function to process transcript files into summary stat
## get median (percentage usage or tpm) across gtex individuals
## for each tissue, transcript, gene triplet
process.transcript.file = function(filename, remove.unexpressed = F) {
    dt = fread(paste("zcat", filename))
    colnames(dt)[3:ncol(dt)] = sapply(colnames(dt)[3:ncol(dt)], fix.indids)
    ## reduce the data table to the samples in egtex and make in long
    dt = dt[, colnames(dt) %in% c('transcript_id','gene_id',map$Collaborator.Sample.ID), with = F]
    dt.melted = melt(dt, id.vars = c(1,2), variable.name = "ind", value.name = "stat")
    rm(dt)
    ## add tissue names
    dt.melted[, tissue := map[match(dt.melted$ind, map$Collaborator.Sample.ID), 'Tissue.Site.Detail']]
    ## get median percentages for each transcript, gene, tissue triplet
    if (remove.unexpressed) {
        ## first removing samples that have no expression for a gene (all transcript percentages set to zero)
        dt.median = dt.melted %>% group_by(gene_id, tissue, ind) %>% filter(sum(stat) > 0) %>% ungroup() %>%
            group_by(tissue, gene_id, transcript_id) %>% summarize(median = median(stat))
    } else {
        dt.median = dt.melted %>% group_by(tissue, gene_id, transcript_id) %>% summarize(median = median(stat))
    }
    return(dt.median)
}

############################

## process map to the set of columns that are useful and reformat tissue.names
map = read.table(mapfile, header = T, stringsAsFactors = F, sep = "\t")
map = map[, c('Collaborator.Sample.ID', 'Tissue.Site.Detail')]
map$Tissue.Site.Detail = format.tissue.names(map$Tissue.Site.Detail)

## process the transcript files into medians across individuals for each transcript, gene, tissue triplet
perc = process.transcript.file(percfile, remove.unexpressed = T)
tpm = process.transcript.file(tpmfile)

## merge them and only keeped the merged one
setnames(perc, "median", "perc")
setnames(tpm, "median", "tpm")
perc = data.table(perc)
tpm = data.table(tpm)
setkey(perc, transcript_id, gene_id, tissue)
setkey(tpm, transcript_id, gene_id, tissue)
transcripts = tpm[perc]
transcripts$tpm = ifelse(is.na(transcripts$tpm), 0, transcripts$tpm)
rm(perc)
rm(tpm)

## also take medians (of the medians) across tissues
transcripts.summary = transcripts[, .(median.tpm = median(tpm),
                                      mad.tpm = mad(tpm),
                                      tpm.gt0 = sum(tpm > 0),
                                      median.perc = median(perc),
                                      mad.perc = mad(perc),
                                      perc.gt0 = sum(perc > 0)), by = .(gene_id, transcript_id)]
setkey(transcripts.summary, gene_id)

## restrict to genes of interest
## get gencode, hgnc mapping: this code was taken from gene_selection/select_eqtl_genes.R
gencode = read.table(gzfile(gencodefile), header = F, stringsAsFactors = F, sep = "\t")
gene.map = data.table(str_split_fixed(gencode[,9], ";? ", n = 11))
stopifnot(nrow(unique(gene.map[,c(1,5,9), with = F])) == 1) # sanity checks
gene.map = unique(gene.map[, c(2,6,10), with = F])
setnames(gene.map, c("ENSG", "GENE_TYPE", "HGNC"))
## get selected genes
genes = read.table(genefile, header = T, stringsAsFactors = F)[,1]
selected = gene.map[HGNC %in% genes, ]

transcripts.summary.subset = transcripts.summary[gene_id %in% selected$ENSG,]

## adding hgnc names and transcript rank
transcripts.summary.subset[, HGNC := gene.map$HGNC[match(gene_id, gene.map$ENSG)]]
transcripts.summary.subset[, rank := rank(-median.perc, ties.method = "random"), by = .(gene_id, HGNC)]

## write subsetted summary to file
write.table(transcripts.summary.subset, outfile, quote = F, sep = "\t", row.names = F, col.names = T)

## small updates for nicer plotting
## make rank a factor with correct level order
transcripts.summary.subset$rank = factor(transcripts.summary.subset$rank,
                                         levels = sort(unique(transcripts.summary.subset$rank), decreasing = T))
## remove transcripts with median percentage <= 1
transcripts.summary.subset = transcripts.summary.subset[median.perc > 1,]

save.image(paste0('process_transcripts', suffix, '.RData'))

## plot transcripts per gene
maxrank = length(levels(transcripts.summary.subset$rank))
colors = c('mediumorchid3','dodgerblue2','springgreen3','goldenrod2','orangered2',rep('darkgrey', maxrank - 5))
names(colors) = as.character(c(1:maxrank))

pdf(paste0(dir, '/figures/transcript.usage', suffix, '.pdf'), width = 8, height = 8)

selected.hgnc = sort(selected$HGNC)

for (i in seq(1,length(selected.hgnc),10)) {
    end = min(i+9, length(selected.hgnc))
    plot.df = transcripts.summary.subset[HGNC %in% selected.hgnc[i:end], ]
    p = ggplot(plot.df, aes(x = HGNC, y = median.perc, colour = rank, group = rank)) +
        geom_pointrange(aes(ymax = median.perc + mad.perc, ymin = median.perc - mad.perc),
                        position = position_dodge(width = 0.7)) +
        geom_vline(aes(xintercept = as.numeric(factor(HGNC)) + 0.5), colour = 'darkgrey') +
        xlab('') + ylab('Transcript usage (%, median +/- MAD across tissues)') +
        scale_colour_manual(values = colors) +
        guides(colour = FALSE) + coord_flip() + theme_bw() +
        theme(panel.grid.major.y = element_blank())
    print(p)
}

dev.off()

