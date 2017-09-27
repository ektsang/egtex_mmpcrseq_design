#!/usr/bin/env Rscript

## author: Emily Tsang

## make some plots of the final primers included in the design

library(ggplot2)
library(plyr)

dir = Sys.getenv('EGTEXDIR')

parse.vep = function(sitedf) {
    vepdf = data.frame(varname = character(), HGNC = character(), VEP = character(), stringsAsFactors = F)
                                        # yes, this is probably slow, but it works...
    for (i in 1:nrow(sitedf)) {
        var = sitedf[i, 'varname']
        gene = sitedf[i, 'HGNC']
        vep = sitedf[i, 'VEP']
        vepSplit = unique(strsplit(vep, "[,&]")[[1]])
        n = length(vepSplit)
        newdf = data.frame(varname = rep(var, n), HGNC = rep(gene, n), VEP = vepSplit)
        vepdf = rbind(vepdf, newdf)
    }
    return(vepdf)
}

get.reason = function(columns) {
    ## if there are multiple reasons, pick the first encountered
    reason = NA
    index = which(!is.na(columns))[2]
    if (index == 2){
        reason = 'Cancer'
    } else if (index == 3) {
        reason = 'eQTL'
    } else if (index == 4) {
        reason = 'Other'
    }
    return(reason)
}


## read in pools
pools = read.table(paste0(dir, '/primerDesign/final.pools.txt'), header = F, stringsAsFactors = F)
colnames(pools) = c('pool','site','fseq','ftemp','rseq','rtemp','amplicon','length','annotation')

## read in reasons for gene selection
genes = read.table(paste0(dir, '/geneSelection/selected.genes.txt'), header = T, stringsAsFactors = F)
genes = genes[, -3]
genes2 = read.table(paste0(dir, '/geneSelection/selected.genes.v2.txt'), header = T, stringsAsFactors = F)
colnames(genes2) = c('gene','n.eqtl.effects')
genes2 = cbind(genes2, data.frame(sources = NA, requested = NA))
genes = rbind(genes, genes2)
rm(genes2)

## read in input sites then subset them to the ones that were selected
sites = read.table(paste0(dir, '/selectedSites/genes.intersect.WGS79.OMNI86.selected.variants.ordered.firstDesign.txt'),
                   header = T, stringsAsFactors = F)
sites = rbind(sites,
              read.table(paste0(dir, '/selectedSites/genes.v2.intersect.WGS79.OMNI86.selected.variants.ordered.firstDesign.txt'),
                         header = T, stringsAsFactors = F))
sites = sites[sites$varname %in% pools$site, ]

sites.summary = ddply(sites, .(HGNC), summarize, nsites = length(varname))

## identify the reason each gene was chosen and intersect with the sites data frames
genes$reason = apply(genes, 1, get.reason)
genes$colour = genes$reason
genes$colour[genes$reason == 'eQTL' & genes$n.eqtl.effects == 1] = 'eQTL (tissue-specific)'
genes$colour[genes$reason == 'eQTL' & genes$n.eqtl.effects == 44] = 'eQTL (ubiquitous)'
sites.summary = merge(sites.summary, genes, by.x = 'HGNC', by.y = 'gene')
sites = merge(sites, genes, by.x = 'HGNC', by.y = 'gene')

pdf(paste0(dir,'/figures/selected.variants.with.primers.characteristics.pdf'), height = 5, width = 5)

mytheme = theme_bw() + theme(axis.text = element_text(size = 11), axis.title = element_text(size = 12))

## gene selection (by gene count)
gene.colors = c('mediumorchid4','dodgerblue','dodgerblue3','dodgerblue4','springgreen4')
ggplot(sites.summary, aes(x = reason, fill = colour)) + geom_bar() +
    mytheme + xlab('') + ylab('Number of genes') +
    scale_fill_manual(values = gene.colors, name = '') +
    theme(legend.key = element_blank(),
          legend.position = c(0.7, 0.8))

## gene selection (by site count)
ggplot(sites, aes(x = reason, fill = colour)) + geom_bar() +
    mytheme + xlab('') + ylab('Number of sites') +
    scale_fill_manual(values = gene.colors, name = '') +
    theme(legend.key = element_blank(),
          legend.position = c(0.7, 0.8))

## amplicon length
ggplot(pools, aes(x = length)) + geom_histogram() +
    mytheme + xlim(150,430) + xlab('Amplicon length') +
    ylab('Number of sites')

## exon vs transcript design
pools$anno.short = 'exon'
pools$anno.short[grep('transcript', pools$annotation)] = 'transcript'
ggplot(pools, aes(x = anno.short)) + geom_bar() +
    mytheme + xlab('Primers designed within') +
    ylab('Number of sites')

## transcript type
ggplot(sites, aes(x = type)) + geom_bar() +
    mytheme + xlab('Transcript type') +
    ylab('Number of genes') +
    theme(axis.text.x = element_text(angle = 10, hjust = 1))

## number of sites per gene
ggplot(sites.summary, aes(x = nsites)) + geom_bar() +
    mytheme + xlab('Number of sites per gene') +
    ylab('Number of genes')

## numer of sites per gene, colored by selection reason
ggplot(sites.summary, aes(x = nsites, fill = colour)) + geom_bar() +
    mytheme + xlab('Number of sites per gene') +
    ylab('Number of genes') +
    scale_fill_manual(values = gene.colors, name = '') +
    theme(legend.key = element_blank(),
          legend.position = c(0.7, 0.8))

## number of hets per site
percData = sum(sites$nhet)/(86*nrow(sites))
ggplot(sites, aes(x = nhet)) + geom_histogram() +
    mytheme + xlab('Number of heterozygours individuals') +
    ylab('Number of sites') +
    annotate("text", x = 20, y = 75,
             label = paste0("Proportion of all site,individual pairs\nwhere the individual is het: ",
                            signif(percData, 2)))

## genotyping?
sites$gt[is.na(sites$gt)] = 'none'
ggplot(sites, aes(x = gt)) + geom_bar() +
    mytheme + xlab('Genotyping in addition to WGS') +
    ylab('Number of sites')

## VEP annotations
sites.vep = parse.vep(sites)
ggplot(sites.vep, aes(x = VEP)) + geom_bar() +
    mytheme + xlab('VEP annotations') +
    ylab('Number of sites (each >= 1 annotation)') +
    theme(axis.text.x = element_text(angle = 55, hjust = 1))

## transcript expression
ggplot(sites, aes(median.TPM+ 0.01)) + geom_histogram() +
    scale_x_log10() + mytheme +
    xlab('Median TPM across tissues') + ylab('Number of sites')

## transcript usage
ggplot(sites, aes(median.perc)) + geom_histogram() +
    mytheme + xlab('Median % transcript usage across tissues') +
    ylab('Number of sites')

dev.off()
