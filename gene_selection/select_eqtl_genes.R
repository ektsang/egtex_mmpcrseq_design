#!/usr/bin/env Rscript

## author: Emily Tsang

## Read in metasoft output and select a diverse set of genes to assay with mmPCR-seq.
## Include a combination of tissue-specific, ubiquitous, and heterogeneous eQTLs.
## Also, include some where Meta-soft is not confident.

## The environment variables should be set globally, e.g. in bashrc ##
dir = Sys.getenv('GTEXv6')
resdir = Sys.getenv('EGTEXDIR')
######################################################################

library(stringr)
library(data.table)

eqtldir = paste0(dir, '/eqtl_updated_annotation')
gencodefile = paste0(dir, '/reference_files/gencode.v19.genes.v6p.patched_contigs.gtf.gz')
metasoftfile = paste0(eqtldir, '/Metasoft_Output_v6p.txt')
variantfile = paste0(dir, '/eqtl_data/eQTLInputFiles/GTEx_Analysis_2015-01-12_eQTLInputFiles_snpPositions.txt')
effectsizefile = paste0(dir, '/lappalainen_lab/eqtl_effect_size_estimates.txt')

## Cancer & other gene list
genedir = paste0(resdir, '/geneSelection')
cancerfile = paste0(genedir, '/cancer/genes.ordered.txt')
outlierfile = paste0(genedir, '/other/medz.outliers.txt')

### Load files that are needed for translating gene names and variant ids
gencode = read.table(gzfile(gencodefile), header = F, stringsAsFactors = F, sep = "\t")
gene.map = data.table(str_split_fixed(gencode[,9], ";? ", n = 11))
stopifnot(nrow(unique(gene.map[,c(1,5,9), with = F])) == 1) # sanity checks
gene.map = unique(gene.map[, c(2,6,10), with = F])
setnames(gene.map, c("ENSG", "GENE_TYPE", "HGNC"))

variant.map = fread(variantfile, header = T)

### Load effect size data
effectsize = fread(effectsizefile, header = T)

### Load metasoft results into a data table
## process header
metasoft.header = str_split_fixed(readLines(metasoftfile, n = 1), "\t", n = 17)[1,-17]
## add tissue names
tissue.names = scan(paste0(eqtldir, '/Metasoft_tissue_order.txt'), what = character())
metasoft.header = c(metasoft.header, paste0("PVALUE.", tissue.names), paste0("MVALUE.", tissue.names))
## actually read in metasoft results
metasoft = fread(metasoftfile, sep = "\t", skip = 1, col.names = metasoft.header, drop = c(105))

### Add columns with the chr, pos, gene name (HGNC), and remove unnecessary ones
## get variant is, gene pairs
rsid.gene = data.table(str_split_fixed(metasoft$RSID, ",", n = 2))
nrow.before = nrow(rsid.gene)
setnames(rsid.gene, c("ID","ENSG"))
## merge with variant positions
setkey(rsid.gene, ID)
setkey(variant.map, ID)
rsid.gene = variant.map[rsid.gene]
## merge with gene names
setkey(rsid.gene, ENSG)
setkey(gene.map, ENSG)
rsid.gene = gene.map[rsid.gene]
stopifnot(nrow(rsid.gene) == nrow.before)
## remove those that did not appear in the variant map (all chr X, see commented out line below)
## also remove non protein-coding/lincRNA genes
#unique(str_split_fixed(tmp[is.na(CHROM),ID], "_", n=2)[,1]) # checking that they are all chr X
rsid.gene = rsid.gene[!is.na(CHROM) & GENE_TYPE %in% c("protein_coding", "lincRNA"), ]
rsid.gene[, RSID := paste(ID, ENSG, sep = ",")]
## merge with main metasoft data table and remove unnecessary columns
setkey(rsid.gene, RSID)
setkey(metasoft, RSID)
keepindx = c(grep("RE2", colnames(metasoft)), grep("MVALUE", colnames(metasoft)))
cols2keep = c(colnames(rsid.gene)[-ncol(rsid.gene)], colnames(metasoft)[keepindx]) # remove rsid column
metasoft = metasoft[rsid.gene, cols2keep, with = F]

### Mark SNP-gene pairs with the number of tissues that pass m-value >= 0.9
mval.cols = grep("MVALUE", colnames(metasoft))
metasoft[, NEFFECTS := rowSums(metasoft[, mval.cols, with = F] >= 0.9, na.rm = T)]

### Add rank of metasoft results by RE2 p-value (from most to least significant)
## both at the gene-snp pair level and at the gene level
metasoft[, RANK := order(metasoft[,PVALUE_RE2])]
generanks = metasoft[, .(MIN_RANK = min(.SD[,RANK])), by = ENSG]
generanks[, GENE_RANK := rank(MIN_RANK)]
setkey(metasoft, ENSG)
setkey(generanks, ENSG)
metasoft = generanks[metasoft,]
setkey(metasoft, GENE_RANK)

### Some plotting
my.tissue.hist = function(vals, title = '') {
    hist(vals, breaks = seq(-0.99, 44.99, 1), xlab = 'Number of tissues',
         main = title, col = 'dodgerblue3', border = 'white', xlim = c(0, 44))
}

pdf(paste0(resdir, '/figures/metasoft.smile.plots.pdf'), height = 5, width = 5)

my.tissue.hist(metasoft$NEFFECTS, title = 'All autosomal protein-coding/lincRNA genes')
my.tissue.hist(metasoft[GENE_RANK <= 500, NEFFECTS], title = 'Top 500 protein-coding/lincRNA genes')
my.tissue.hist(metasoft[GENE_RANK <= 1000, NEFFECTS], title = 'Top 1,000 protein-coding/lincRNA genes')
my.tissue.hist(metasoft[GENE_RANK <= 2000, NEFFECTS], title = 'Top 2,000 protein-coding/lincRNA genes')
my.tissue.hist(metasoft[GENE_RANK <= 5000, NEFFECTS], title = 'Top 5,000 protein-coding/lincRNA genes')
my.tissue.hist(metasoft[GENE_RANK <= 8000, NEFFECTS], title = 'Top 8,000 protein-coding/lincRNA genes')


### Look at overlap with cancer and outlier genes
cancer = read.table(cancerfile, header = T, stringsAsFactors = F)
cancer.picked = cancer$gene[cancer$pass == 1]
cancer.notpicked = cancer$gene[cancer$pass == 0]
outlier = read.table(outlierfile, header = F, stringsAsFactors = F,
                     col.names = c('ENSG','IND','NTISSUE','Z'))
outliers = outlier$ENSG
## Add columns to metasoft data table indicating if the respective genes were in these lists
metasoft[, cancer.gene := ifelse(HGNC %in% cancer.picked, yes = 1, no = 0)]
metasoft[, cancer.other := ifelse(HGNC %in% cancer.notpicked, yes = 1, no = 0)]
metasoft[, outlier.gene := ifelse(ENSG %in% outliers, yes = 1, no = 0)]

## Look at distribution of gene ranks for each gene group
hist(unique(metasoft[cancer.gene == 1, GENE_RANK]), xlab = 'Gene rank', main = 'Cancer genes',
     col = 'mediumorchid4', border = 'white', breaks = 20)
abline(v = 2000, col = 'red')
hist(unique(metasoft[cancer.other == 1, GENE_RANK]), xlab = 'Gene rank', main = 'Cancer genes not picked',
     col = 'mediumorchid4', border = 'white', breaks = 20)
abline(v = 2000, col = 'red')
hist(unique(metasoft[outlier.gene == 1, GENE_RANK]), xlab = 'Gene rank', main = 'Outlier genes',
     col = 'mediumorchid4', border = 'white', breaks = 20)
abline(v = 2000, col = 'red')
## Look at the distribution of the number of tissues with an effect (limiting to genes in the top 2000)
my.tissue.hist(metasoft[GENE_RANK <= 2000 & cancer.gene == 1, NEFFECTS], title = 'Cancer genes in the top 2000 eQTL genes')
my.tissue.hist(metasoft[GENE_RANK <= 2000 & outlier.gene == 1, NEFFECTS], title = 'Outlier genes in the top 2000 eQTL genes')

dev.off()

### Restrict metasoft results to the top SNP per gene and the top 2000 genes
metasoft.subset = metasoft[GENE_RANK <= 2000 & RANK == MIN_RANK,]
## change the chr files to have the "chr" (for downstream matching)
metasoft.subset$CHROM = paste0('chr', metasoft.subset$CHROM)
## Write that list of genes
metasoft.gene.file = paste0(genedir, '/other/metasoft_scratch/top2000.metasoft.txt')
write.table(metasoft.subset$HGNC, metasoft.gene.file, quote = F, row.names = F, col.names = F)

### Run commands to get the information on the most commmon het per gene.
### These take a couple minutes to run and this only needs to be repeated if the top2000.metasoft list changes
processed.metasoft.genes = paste0(genedir, '/other/metasoft_scratch/top2000.metasoft.exons.txt')
mycommand1 = paste('python ../experimental_layout/retrieve_genes.py',
                   '/mnt/lab_data/montgomery/shared/annotations/gencode.v19.annotation.gtf',
                   metasoft.gene.file, '| sort -k1,1 -k2,2n >', processed.metasoft.genes)
system(mycommand1)
# NOTE: used an old version of the get_sites_in_exons.py script. just use the output file that already exists
#intersection.prefix = paste0(genedir, '/other/metasoft_scratch/top2000.metasoft.exons.intersect.WGS.79ind')
#genotypeFile = paste0(resdir,'/sampleInfo/genotypes/WGS_V7_79ind_SNP_PASS_ONLY.bed')
#mycommand2 = paste('python ../experimental_layout/get_sites_in_exons.py -o',
#                   intersection.prefix,
#                   genotypeFile,
#                   processed.metasoft.genes)
#system(mycommand2)

### Read in results to get the most common het per gene
gene.hets = fread(paste0(intersection.prefix, '.orderedSites.byGene.txt'))
common.hets = gene.hets[pick.order == 1, .(gene, site.chr, site.pos)]
## Merge them with the subsetted metasoft results
setkey(common.hets, gene)
setkey(metasoft.subset, HGNC)
metasoft.subset = metasoft.subset[common.hets]

### Count the nubmer of individuals that are het both for the top eQTL SNP and the most common coding SNP
## Read in the genotype data
genotypes = fread(genotypeFile)
setnames(genotypes, c('chrom', 'pos0', 'pos1', 'nhets', 'hetinds'))
setkey(genotypes, chrom, pos1)
## actually do counting and add results into a new column
hetcount.vector = rep(NA, nrow(metasoft.subset))
for (i in 1:nrow(metasoft.subset)) {
    inds.eqtl = genotypes[.(metasoft.subset[i, CHROM], metasoft.subset[i, POS]), hetinds]
    inds.coding = genotypes[.(metasoft.subset[i, site.chr], metasoft.subset[i, site.pos]), hetinds]
    inds.eqtl.vector = strsplit(inds.eqtl, ',')[[1]]
    inds.coding.vector = strsplit(inds.coding, ',')[[1]]
    numinds = sum(!is.na(intersect(inds.eqtl.vector, inds.coding.vector)))
    hetcount.vector[i] = numinds
}
metasoft.subset[, num.double.hets := hetcount.vector]

### Gene selection
rowindx = metasoft.subset$num.double.hets >= 10 & (metasoft.subset$cancer.other == 1 | metasoft.subset$outlier.gene == 1 | metasoft.subset$NEFFECTS >= 43)
rowindx2 = metasoft.subset$num.double.hets >= 25 & metasoft.subset$NEFFECTS == 1

### Some plotting to get a sense of what's left
pdf(paste0(resdir, '/figures/eqtl.gene.selection.qc.pdf'), height = 5, width = 5)

## distribution of the number of individuals het for both the best eQTL SNP and the most common (in terms of hets) coding SNP
hist(metasoft.subset$num.double.hets, breaks = 20, main = '', xlab = 'Number of individuals het for best eQTL and coding SNPs')
hist(metasoft.subset[rowindx | rowindx2, num.double.hets], breaks = 20,
     main = '', xlab = 'Number of selected individuals het for best eQTL and coding SNPs')
## distribution of effects
my.tissue.hist(metasoft.subset[rowindx, NEFFECTS],
               title = 'eQTL genes with >= 10 hets\n(other) cancer or outlier genes or 43+ effects')
my.tissue.hist(metasoft.subset[rowindx | rowindx2, NEFFECTS],
               title = 'eQTL genes with >= 10 hets\n(other) cancer or outlier genes or 43+ effects\nOR >= 25 hets and 1 effect')
## gene ranks
hist(metasoft.subset[rowindx, GENE_RANK], xlab = 'Gene rank',
     main = 'eQTL genes with >= 10 hets\n(other) cancer or outlier genes or 43+ effects')
hist(metasoft.subset[rowindx | rowindx2, GENE_RANK], xlab = 'Gene rank',
     main = 'eQTL genes with >= 10 hets\n(other) cancer or outlier genes or 43+ effects\nOR >= 25 hets and 1 effect')

dev.off()

### Number of total genes in different classes
cat("Number of selected cancer genes with at least 10 double hets:",
    sum(metasoft.subset$cancer.gene & metasoft.subset$num.double.hets >= 10), "\n")
selected = rowindx | rowindx2
cat("Number of selected eQTL genes that were not already selected cancer genes:",
    sum(selected & !(metasoft.subset$cancer.gene)), "\n")
cat("Number of selected eQTL genes that are other cancer genes:",
    sum(selected & metasoft.subset$cancer.other == 1), "\n")
cat("Number of selected eQTL genes that are outlier genes:",
    sum(selected & metasoft.subset$outlier.gene == 1), "\n")
cat("Number of selected eQTL genes with tissue specific effects:",
    sum(selected & metasoft.subset$NEFFECTS == 1), "\n")
cat("Number of selected eQTL genes with ubiquitous effects:",
    sum(selected & metasoft.subset$NEFFECTS == 44), "\n")
cat("Number of selected eQTL genes with effects in 43 tissues:",
    sum(selected & !(metasoft.subset$cancer.other == 1 |
                     metasoft.subset$outlier.gene == 1 |
                     metasoft.subset$NEFFECTS == 1 |
                     metasoft.subset$NEFFECTS == 44)), "\n")

### Print out list of selected genes with extra columns indicating why they were selected
selected.genes = metasoft.subset[selected, .(HGNC, selected = 1, cancer.gene, cancer.other, outlier.gene, NEFFECTS)]
write.table(selected.genes, paste0(genedir, '/other/eqtl.selected.txt'), col.names = T, row.names = F, quote = F, sep = "\t")


### Look at effect size data to make sure we include some positive controls
effectsize = read.table(effectsizefile, header = T, stringsAsFactors = F, sep = "\t")
#keep = !is.na(effectsize$effect_size_estimate) & !is.na(effectsize$ase_effect_size_estimate)
#plot(effectsize$effect_size_estimate[keep], effectsize$ase_effect_size_estimate[keep])
effectsize = effectsize[effectsize$gene %in% metasoft.subset[selected, ENSG],]

## how many genes have effect size estimates for ase only, eqtl only, or both
cat('Number of genes with only ASE effect sizes:\n',
    length(unique(effectsize$gene[!is.na(effectsize$ase_effect_size_estimate) & is.na(effectsize$effect_size_estimate)])), '\n')
cat('Number of genes with only eQTL effect sizes:\n',
    length(unique(effectsize$gene[!is.na(effectsize$effect_size_estimate) & is.na(effectsize$ase_effect_size_estimate)])), '\n')

## for the subset of genes with both effect size estimates, are there cases for which they are highly correlated?
dev.new()
keep = !is.na(effectsize$effect_size_estimate) & !is.na(effectsize$ase_effect_size_estimate)
plot(effectsize$effect_size_estimate[keep], effectsize$ase_effect_size_estimate[keep])

save.image('select_eqtl_genes.RData')

## new genes for when the first set of genes didn't yield enough primers (16.12.21)
rowindx3 = !(selected | metasoft.subset$cancer.gene) & metasoft.subset$num.double.hets >= 20 & metasoft.subset$NEFFECTS == 1
rowindx4 = !(selected |metasoft.subset$cancer.gene | rowindx3) & metasoft.subset$num.double.hets >= 25
set.seed(1)
rowindx4[sample(which(rowindx4), sum(rowindx4) - 31)] = FALSE

selected2 = rowindx3 | rowindx4
selected.genes2 = metasoft.subset[selected2, .(HGNC, NEFFECTS)]
write.table(selected.genes2, paste0(genedir, '/selected.genes.v2.txt'), col.names = T, row.names = F, quote = F, sep = "\t")
