#!/usr/bin/env Rscript

## author: Emily Tsang

# script to pick cancer genes to include
setwd('~/projects/egtex_mmPCR/cancer_genes/')

# cosmic database of cureated genes
cosmic = read.table('cosmic.txt', header = F, col.names = c('gene', 'nsamples', 'nmutations', 'npapers'), stringsAsFactors = F)

# read in data from a variety of cancer panels and keep genes that appear more than once
panels = read.table('cancer.wtrusight.wampliseq.panels.txt', header = F, stringsAsFactors = F, col.names = c('gene','n.panels'))
panel.genes = panels$gene[panels$n.panels > 1]

# read in foundation medicine panels
fone.heme = scan('foundationOneHeme.txt', what = character())
fone = scan('foundationOne.txt', what = character())

# read in ongen ASE genes
ongen = scan('ongen.excessASE.txt', what = character())

# read in stamp genes (Stanford's cancer panel)
stamp = scan('stampv2.txt', what = character())

# read in stanford's gols standard panel from J. Merker
stanford = read.table('stanford.cancer.gold.txt', header = T, stringsAsFactors = F, sep = "\t")[,1]

# read in gene lists from TCGA breast cancer paper and from multiple cancer paper
breast.all = scan('TCGAbreastCancer_SMGall.txt', what = character())
breast.subtype = scan('TCGAbreastCancer_SMGsubtypes.txt', what = character())
multi.top = scan('multiplatformAnalysis12cancerTypes_top40SMGs.txt', what = character())
multi.other = scan('multiplatformAnalysis12cancerTypes_otherSMGs.txt', what = character())

other = unique(c(panel.genes, fone.heme, fone, ongen, stamp, stanford, breast.all, breast.subtype, multi.top, multi.other))

# add rows for the genes in the lists above that aren't in cosmic
cosmic$cosmic = 1
extra = other[!(other %in% cosmic$gene)]
cosmic = rbind(cosmic, data.frame(gene = extra, nsamples = NA, nmutations = NA, npapers = NA, cosmic = 0, stringsAsFactors = F))

# add columns to indicate presence in each of the lists above
cosmic$panels = ifelse(cosmic$gene %in% panel.genes, yes = 1, no = 0)
cosmic$fone.heme = ifelse(cosmic$gene %in% fone.heme, yes = 1, no = 0)
cosmic$fone = ifelse(cosmic$gene %in% fone, yes = 1, no = 0)
cosmic$ongen = ifelse(cosmic$gene %in% ongen, yes = 1, no = 0)
cosmic$stamp = ifelse(cosmic$gene %in% stamp, yes = 1, no = 0)
cosmic$stanford = ifelse(cosmic$gene %in% stanford, yes = 1, no = 0)
cosmic$breast.all = ifelse(cosmic$gene %in% breast.all, yes = 1, no = 0)
cosmic$breast.subtype = ifelse(cosmic$gene %in% breast.subtype, yes = 1, no = 0)
cosmic$multi.top = ifelse(cosmic$gene %in% multi.top, yes = 1, no = 0)
cosmic$multi.other = ifelse(cosmic$gene %in% multi.other, yes = 1, no = 0)

# add score
cosmic$score = rowSums(cosmic[,5:15])

# hack to prioritize likely regulatory genes higher
# add 4 to score of each of the STAMP genes that are annotated as "promoter"
reggenes = c('TERT', 'SDHD', 'DPH3', 'PLEKHS1')
cosmic$score[cosmic$gene %in% reggenes] = cosmic$score[cosmic$gene %in% reggenes] + 4 

cosmic$pass = as.numeric(cosmic$score > 1)

# reorder cosmic by score (within that, alphabetically)
cosmic = cosmic[order(cosmic$gene), ]
cosmic = cosmic[order(cosmic$score, decreasing = T), ]

# for every pair, count the overlap
overlap = function(indices) {
  i1 = indices[1]
  i2 = indices[2]
  cat(colnames(cosmic)[i1], ":", sum(cosmic[, i1], na.rm = T), "\n")
  cat(colnames(cosmic)[i2], ":", sum(cosmic[, i2], na.rm = T), "\n")
  overlap = sum(cosmic[, i1] & cosmic[, i2])
  cat('overlap: ', overlap , 'jaccard:', overlap/(sum(cosmic[, i1] | cosmic[, i2])), "\n\n")
}

pairs = combn(c(17, 5:15), 2)
apply(pairs, 2, overlap)

# add column with sources in comma-separated list
cosmic$sources = apply(cosmic[,c(5:15)], 1, function(x) paste(which(x == 1), collapse = ","))

# first write comments to the top of the file
infostring = "# sources
# 1- cosmic
# 2- at least two of the following panels: GeneDx Breat and Ovarian or comprehensive cancer panels, Illumina TruSeq and TruSight cancer panels, Mayo Clinic Solid Tumor Targeted Cancer Gene Panel, University of Washington BROCA and ColoSeq panels, Ion AmpliSeq Comprehensive Cancer Panel
# 3- Foundation One Heme panel
# 4- Foundation One Solid Tumor panel
# 5- Genes with significantly more ASE in colorectal cancer (Ongen et al. Nature, 2014)
# 6- STAMP v2 panel from Stanford
# 7- Stanford Cancer gold standard panel
# 8- Significantly mutated genes in breast cancer (TCGA. Nature, 2012)
# 9- Significantly mutated genes in specific breast cancer subtypes (TCGA. Nature, 2012)
# 10- Top 40 significantly mutated genes from analysis of 12 cancer types (Hoadley et al. Cell, 2014)
# 11- Other significantly mutated gens from 12 cancer types (Hoadley et al. Cell, 2014)"

write(infostring, 'genes.ordered.txt')

# write to file
write.table(cosmic[, c('gene','pass','sources')], 'genes.ordered.txt', col.names = T, row.names = F, quote = F, sep = "\t", append = T)
