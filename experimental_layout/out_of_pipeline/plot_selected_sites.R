#!/usr/bin/env Rscript

## author: Emily Tsang

library(ggplot2)
library(plyr)
library(cowplot)

dir = '/srv/scratch/etsang/projects/egtex_mmPCR_2015'
sitesfile = paste0(dir, '/selectedSites/genes.intersect.WGS79.OMNI86.selected.variants.ordered.v4.txt')

## FUNCTIONS
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

get.closest = function(sitedf) {
    closest = rep(NA, nrow(sites))
    ## first
    closest[1] = ifelse(sitedf[1, 'var.chr'] == sitedf[2, 'var.chr'],
                        yes = sitedf[2, 'var.pos0'] - sitedf[1, 'var.pos0'],
                        no = -1)
    for (i in 2:(nrow(sitedf)-1)){
        opt1 = ifelse(sitedf[i - 1, 'var.chr'] == sitedf[i, 'var.chr'],
                      yes = sitedf[i, 'var.pos0'] - sitedf[i-1, 'var.pos0'], no = -1)
        opt2 = ifelse(sitedf[i, 'var.chr'] == sitedf[i+1, 'var.chr'],
                      yes = sitedf[i+1, 'var.pos0'] - sitedf[i, 'var.pos0'], no = -1)
        closest[i] = ifelse(min(opt1, opt2) != -1, yes = min(opt1, opt2), no = max(opt1, opt2))
    }
    ## last
    closest[nrow(sitedf)] = ifelse(sitedf[nrow(sitedf) - 1, 'var.chr'] == sitedf[nrow(sitedf), 'var.chr'],
                        yes = sitedf[nrow(sitedf), 'var.pos0'] - sitedf[nrow(sitedf) - 1, 'var.pos0'],
                        no = -1)
    return(closest)
}

###############################################

sites = read.table(sitesfile, header = T, stringsAsFactors = F, sep = "\t")
sites$gt[is.na(sites$gt)] = "not gt"
sites$dist.closest = get.closest(sites) # add distance to nearest variant (assumes the file is sorted by chr, pos)

sites.vep = parse.vep(sites)

## max
sites.max = ddply(sites, .(HGNC), summarise, max.nhet = max(nhet), min.rank = min(rank),
                  max.TPM = max(median.TPM), max.perc = max(median.perc), nsites = length(varname),
                  types = paste(sort(unique(type)), collapse = ","))

pdf(paste0(dir, '/figures/selected.variants.characteristics.v4.pdf'), width = 9, height = 9)

## basic plots
p1 = ggplot(sites, aes(x = nhet)) + geom_histogram() + theme_bw() +
    xlab('Number of heterozygous individuals') + ylab('Number of sites') +
    ggtitle('Heterozygosity of sites passing filter')

p2 = ggplot(sites, aes(x = log2(median.TPM+2))) + geom_histogram() + theme_bw() +
    xlab('log2(Median TPM across tissues + 2)') + ylab('Number of sites') +
    ggtitle('Transcript expression')

p3 = ggplot(sites, aes(x = median.perc)) + geom_histogram() + theme_bw() +
    xlab('Median transcript percentage across tissues') + ylab('Number of sites') +
    ggtitle('Transcript usage')

p4 = ggplot(sites, aes(x = rank)) + geom_bar() + theme_bw() +
    xlab('Rank') + ylab('Number of sites') +
    ggtitle('Transcript rank')

ggdraw() +
    draw_plot(p1, 0, 1/2, 1/2, 1/2) +
    draw_plot(p2, 1/2, 1/2, 1/2, 1/2) +
    draw_plot(p3, 0, 0, 1/2, 1/2) +
    draw_plot(p4, 1/2, 0, 1/2, 1/2)

## picking the best site per gene for various metrics
pm1 = ggplot(sites.max, aes(max.nhet)) + geom_histogram() + theme_bw() +
    xlab('Number of heterozygous individuals') + ylab('Number of genes') +
    ggtitle('Heterozygosity of the best site per gene')

pm2 = ggplot(sites.max, aes(x = log2(max.TPM+2))) + geom_histogram() + theme_bw() +
    xlab('log2(Median TPM across tissues + 2)') + ylab('Number of genes') +
    ggtitle('Highest transcript expression per gene')

pm3 = ggplot(sites.max, aes(x = max.perc)) + geom_histogram() + theme_bw() +
    xlab('Median transcript percentage across tissues') + ylab('Number of genes') +
    ggtitle('Highest transcript usage per gene')

pm4 = ggplot(sites.max, aes(x = min.rank)) + geom_bar() + theme_bw() +
    xlab('Rank') + ylab('Number of genes') +
    ggtitle('Lowest transcript rank')

ggdraw() +
    draw_plot(pm1, 0, 1/2, 1/2, 1/2) +
    draw_plot(pm2, 1/2, 1/2, 1/2, 1/2) +
    draw_plot(pm3, 0, 0, 1/2, 1/2) +
    draw_plot(pm4, 1/2, 0, 1/2, 1/2)


p5 = ggplot(sites.max, aes(x = nsites)) + geom_histogram() + theme_bw() +
    xlab('Number of sites') + ylab('Number of genes') +
    ggtitle('Number of sites per gene')

p6 = ggplot(sites.vep, aes(x = VEP)) + geom_bar() + theme_bw() +
    xlab('') + ylab('Number of sites') +
    ggtitle('Types of variants') +
    theme(axis.text.x = element_text(angle = 55, hjust = 1))

p7 = ggplot(sites, aes(x = dist.closest/1000)) + geom_histogram() + theme_bw() +
    scale_x_log10() + ggtitle('Distance between sites') +
    xlab('Distance to the nearest site (kb)') + ylab('Number of sites')

ggdraw() +
    draw_plot(p5, 0, 1/2, 1/2, 1/2) +
    draw_plot(p6, 1/2, 1/2, 1/2, 1/2) +
    draw_plot(p7, 0, 0, 1/2, 1/2) #+
#    draw_plot(p8, 1/2, 0, 1/2, 1/2)

b1 = ggplot(sites, aes(x = gt)) + geom_bar() + theme_bw() +
    xlab('') + ylab('Number of sites') +
    ggtitle('Genotyping in addition to WGS')

b2 = ggplot(sites, aes(x = factor(intra.exon.plausible))) + geom_bar() + theme_bw() +
    xlab('Plausible intra-exon design') + ylab('Number of sites') +
    ggtitle('Primer design')

b3 = ggplot(sites, aes(x = !is.na(overlap.same))) + geom_bar() + theme_bw() +
    xlab('Exon of interest overlaps other exons in the same gene') +
    ylab('Number of sites') + ggtitle('Intra-gene overlap')

b4 = ggplot(sites, aes(x = !is.na(overlap.other))) + geom_bar() + theme_bw() +
    xlab('Exon of interest overlaps other exons in another gene') +
    ylab('Number of sites') + ggtitle('Inter-gene overlap')

ggdraw() +
    draw_plot(b1, 0, 1/2, 1/2, 1/2) +
    draw_plot(b2, 1/2, 1/2, 1/2, 1/2) +
    draw_plot(b3, 0, 0, 1/2, 1/2) +
    draw_plot(b4, 1/2, 0, 1/2, 1/2)

## some plots about transcript usage
mycols = c('orangered3','goldenrod3','springgreen4','dodgerblue3','mediumorchid4', 'grey30')
type.colors = mycols[1:5]

t1 = ggplot(sites, aes(x = type)) + geom_bar(aes(fill = type)) +
    xlab('') + ylab('Number of sites') + theme_bw() +
    guides(fill = FALSE) + scale_fill_manual(values = type.colors) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

type.names = sort(unique(sites.max$types))
type.comb.colors = rep(mycols[6], length(type.names))
type.comb.colors[grep("lincRNA", type.names, fixed = T)] = mycols[1]
type.comb.colors[grep("protein_coding", type.names, fixed = T)] = mycols[4]
names(type.comb.colors) = type.names

t2 = ggplot(sites.max, aes(x = types)) + geom_bar(aes(fill = types)) +
    xlab('') + ylab('Number of genes') + theme_bw() +
    guides(fill = FALSE) + scale_fill_manual(values = type.comb.colors) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))

t3 = ggplot(sites, aes(x = type, y = log2(median.TPM+2))) + geom_violin(aes(fill = type)) +
    geom_boxplot(width = 0.05) + guides(fill = FALSE) +
    scale_fill_manual(values = type.colors) +
    xlab('') + ylab('log2(Median TPM across tissues + 2)') + theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

t4 = ggplot(sites, aes(x = type, y = median.perc)) + geom_violin(aes(fill = type)) +
    geom_boxplot(width = 0.05) + guides(fill = FALSE) +
    scale_fill_manual(values = type.colors) +
    xlab('') + ylab('Median transcript percentage across tissues') + theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggdraw() +
    draw_plot(t1, 0, 1/2, 1/2, 1/2) +
    draw_plot(t2, 1/2, 0.4, 1/2, 0.6) +
    draw_plot(t3, 0, 0, 1/2, 1/2) +
    draw_plot(t4, 1/2, 0, 1/2, 1/2)

dev.off()

