#!/usr/bin env Rscript

## author: Emily Tsang

## combine the various gene lists and maintain the information about why each gene was chosen

## The environment variables should be set globally, e.g. in bashrc ##
baseDir = Sys.getenv('EGTEXDIR')
v6dir = Sys.getenv('GTEXv6')
#########################################################################

library(stringr)

## selected gene files
dir = paste0(baseDir, '/geneSelection')
cancerfile = paste0(dir, '/cancer/genes.ordered.txt')
eqtlfile = paste0(dir, '/other/eqtl.selected.txt')
requestedfile = paste0(dir, '/other/requested.txt')
## gene annotation file
gencodefile = paste0(v6dir, '/reference_files/gencode.v19.genes.v6p.patched_contigs.gtf.gz')

## small function to get one string to describe the eqtl choice
make.eqtl.string = function(x) {
    eqtl.string = ''
    if (x[1] == 1) {
        eqtl.string = paste0(eqtl.string, ",cancer")
    }
    if (x[2] == 1) {
        eqtl.string = paste0(eqtl.string, ",outlier")
    }
    if (x[3] == 1) {
        eqtl.string = paste0(eqtl.string, ",specific")
    }
    if (x[3] == 44) {
        eqtl.string = paste0(eqtl.string, ",ubiquitous")
    }
    if (x[1] == 0 & x[2] == 0 & x[3] != 1 & x[3] != 44) {
        eqtl.string = ",other"
    }
    ## strip leading comma
    eqtl.string = substr(eqtl.string, 2, nchar(eqtl.string))
    return(eqtl.string)
}

## process gencode file to get list of autosomal, protein-coding genes
gencode = read.table(gzfile(gencodefile), header = F, stringsAsFactors = F, sep = "\t")
gene.map = str_split_fixed(gencode[,9], ";? ", n = 11)
genemap = data.frame(cbind(gene.map[, c(10,6)], gencode[, 1]))
colnames(genemap) = c('gene', 'genetype', 'chrom')
genes.keep = genemap[genemap$genetype %in% c('protein_coding', 'lincRNA') &
                     !(genemap$chrom %in% c('X','Y','MT')), 'gene']

## read in the files and process them as necessary
cancer = read.table(cancerfile, header = T, stringsAsFactors = F)
cancer = cancer[cancer$pass == 1, c('gene','sources')]
## remove gene that is already there under another name
cancer = cancer[cancer$gene != "MLL3",]

eqtl = read.table(eqtlfile, header = T, stringsAsFactors = F)
eqtl = eqtl[, colnames(eqtl) != "selected"]
eqtl$eqtl.source = apply(eqtl[,c('cancer.other','outlier.gene','NEFFECTS')], 1, make.eqtl.string)
eqtl = eqtl[, c('HGNC','eqtl.source','NEFFECTS')]
colnames(eqtl) = c('gene','eqtl.source', 'n.eqtl.effects')
# remove HLA genes
eqtl = eqtl[!(eqtl$gene %in% c("HLA-A","HLA-B","HLA-C","HLA-DPA1")), ]

requested = read.table(requestedfile, header = F, stringsAsFactors = F)[,1]

## combine all gene lists
combined = cancer
combined = merge(combined, eqtl, by = 'gene', all = T)
combined$requested = NA
for (new in requested[!(requested %in% combined$gene)]) {
    combined = rbind(combined, c(new,NA,NA,NA,NA))
}
combined$requested[combined$gene %in% requested] = 1

## order genes alphabetically
combined = combined[order(combined$gene),]

## remove genes that aren't protein coding and autosomal
print('removed:\n')
print(combined[!(combined$gene %in% genes.keep), ])
combined = combined[combined$gene %in% genes.keep, ]


write.table(combined, paste0(dir, '/selected.genes.txt'), col.names = T, row.names = F, quote = F, sep = "\t")
