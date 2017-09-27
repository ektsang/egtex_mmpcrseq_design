#!/usr/bin/env Rscript

# author: Emily Tsang

# some initial plotting of the high frequency het sites that fall within the selected genes

library(plyr)
library(ggplot2)
library(RColorBrewer)

mybar = function(data, ...) {
    barplot(data, col = 'dodgerblue3', border = 'white', ...)
}

qc.plots.helper = function(data, ngenes, title) {
    data = data[data$count >0, ]
    # number of genes without any sites
    nosites = ngenes - nrow(data)
    cat('Number of genes with no sites:', nosites, '\n')
    
    # get the number of sites per gene
    pergene = table(data$count)

    # summarize that data
    counts = c('0' = nosites, pergene)
    # make one category 15+
    countSummary = c(head(counts,15), "15+"=sum(counts[c(16:length(counts))]))

    mybar(countSummary, main = title, xlab = 'Number of sites', ylab = 'Number of genes')
    return(countSummary)
}

qc.plots = function(data, ngenes) {
    any = ddply(data, c('gene'), summarize, count = sum(count))
    # > 10 hets
    min11 = ddply(data[data$hetGroup %in% c('11-20','21-40','>40'), ], c('gene'), summarize, count = sum(count))
    # > 20 hets
    min21 = ddply(data[data$hetGroup %in% c('21-40','>40'), ], c('gene'), summarize, count = sum(count))

    qc.plots.helper(any, ngenes, 'Number of sites per gene')
    qc.plots.helper(min11, ngenes, 'Number of sites per gene > 10 hets')
    qc.plots.helper(min21, ngenes, 'Number of sites per gene > 20 hets')
}

## requires the data frame picked.subset generated below
plot.num.sites = function(geneset) {
    cum.cols = brewer.pal(5, 'Spectral')
    p = ggplot(picked.subset[picked.subset$gene %in% geneset,], aes(x = gene, y = cum.num.people, fill = factor(pick.order))) +
        geom_bar(stat = "identity", position = "dodge") + geom_hline(yintercept = 79, linetype = 2) +
        theme_bw() + ylab('Cumulative number of het individuals') + xlab('Genes') +
        scale_fill_manual(values = rev(cum.cols), name = "Number of sites") +
        theme(legend.key = element_blank(),
              legend.text = element_text(size = textsize),
              legend.title = element_text(size = textsize),
              axis.text = element_text(size = textsize),
              axis.title = element_text(size = textsize))
    return(p)
}

#######################################################################

# read in site by exon intersection
dir = '/srv/scratch/etsang/projects/egtex_mmPCR_2015/selectedSites/'
outdir = '/srv/scratch/etsang/projects/egtex_mmPCR_2015/figures/'
hetCounts = read.table(paste0(dir,'genes.intersect.WGS.79ind.hetcounts.byGene.txt'), header = T, stringsAsFactors = F)
picked = read.table(paste0(dir,'genes.intersect.WGS.79ind.orderedSites.byGene.txt'), header = T, stringsAsFactors = F)

genedf = read.table('/srv/scratch/etsang/projects/egtex_mmPCR_2015/geneSelection/selected.genes.txt',
                    header = T, stringsAsFactors = F)
genes = genedf$gene

ngenes = length(genes)

pdf(paste0(outdir, 'genes.intersect.WGS.79ind.histograms.pdf'), height = 5, width = 7.5)
qc.plots(hetCounts, ngenes)
dev.off()

# comparison with OMNI
hetCountsOmni = read.table(paste0(dir,'genes.intersect.OMNI.86ind.hetcounts.byGene.txt'), header = T, stringsAsFactors = F)
omni.gt10 = ddply(hetCountsOmni[hetCountsOmni$hetGroup %in% c('11-20','21-40','>40'), ], c('gene'), summarize, count = sum(count))
wgs.gt10 = ddply(hetCounts[hetCounts$hetGroup %in% c('11-20','21-40','>40'), ], c('gene'), summarize, count = sum(count))
omni.table = qc.plots.helper(omni.gt10, ngenes, 'OMNI')
wgs.table = qc.plots.helper(wgs.gt10, ngenes, 'WGS')
omni.table = c(omni.table[c(1:10)], "10 +"=sum(omni.table[c(11:length(omni.table))]))
wgs.table = c(wgs.table[c(1:10)], "10 +"=sum(wgs.table[c(11:length(wgs.table))]))
comp.data = data.frame(Categ = c(names(omni.table),names(wgs.table)),
    Count = c(omni.table, wgs.table),
    Platform = rep(c("OMNI","WGS"), each = 11))

pdf(paste0(outdir, 'compare.OMNI.WGS.pdf'), height = 5, width = 5)

comp.data$Categ = factor(comp.data$Categ, levels = comp.data$Categ[1:11])
ggplot(comp.data, aes(x = Categ, y = Count, fill = Platform)) + geom_bar(stat = "identity", position = "dodge", colour = "white") +
    theme_classic() + scale_fill_manual(values = c("dodgerblue", "dodgerblue4")) +
    xlab('Number of heterozygous sites') + ylab('Number of genes') + ggtitle("Number of sites per gene with > 10 hets")
    theme(legend.position = c(0.2, 0.7),
          text = element_text(size = 14),
          legend.text = element_text(size = 14),
          axis.text = element_text(size = 14))

dev.off()

# make gene and hetGroup factors for nice plotting
hetCounts = hetCounts[hetCounts$count > 0, ]
hetCounts = hetCounts[order(hetCounts$gene), ]

hetCounts$hetGroup = factor(hetCounts$hetGroup, levels = c('1','2-10','11-20','21-40','>40'))
hetCounts$gene = factor(hetCounts$gene, levels = unique(hetCounts$gene[order(hetCounts$hetGroup, hetCounts$count, decreasing = T)]))
hetCounts = hetCounts[order(hetCounts$hetGroup, decreasing = T), ]

textsize = 14

pdf(paste0(outdir, 'genes.WGS.num.hets.pdf'), height = 6.5, width = 13)

colors = brewer.pal(5, "Spectral")
ggplot(hetCounts, aes(x = gene, y = count, fill = hetGroup)) + geom_bar(stat = "identity") +
    theme_bw() + scale_fill_manual(values = rev(colors), name = "Number of hets") +
    xlab('Genes') + ylab('Number of sites') +
    theme(legend.position = c(0.9, 0.67),
          legend.key = element_blank(),
          legend.background = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.border = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_text(margin = margin(-15,0,0,0)),
          legend.text = element_text(size = textsize + 6),
          legend.title = element_text(size = textsize + 6),
          axis.text.y = element_text(size = textsize + 6),
          axis.title = element_text(size = textsize + 6))

dev.off()

picked.subset = picked[picked$pick.order <= 5,]
picked.genes = sort(unique(picked.subset$gene))


pdf(paste0(outdir, 'cumulative.sites.pdf'), height = 6, width = 10)

for (i in seq(1,(length(picked.genes) - 9),10)) {
    print(plot.num.sites(picked.genes[i:(i+9)]))
}
## then deal with the leftovers
last = tail(seq(1,(length(picked.genes) - 9),10), 1) + 9
if (length(picked.genes) > last) {
    print(plot.num.sites(picked.genes[(last+1):length(picked.genes)]))
}


dev.off()
