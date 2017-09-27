#!/usr/bin/env Rscript

## author: Emily Tsang

## read in the designed primer pools and format into an order for IDT
## also write to file the sites for which primers were designed, split by pool

dir = Sys.getenv('EGTEXDIR')

poolfile = paste0(dir, '/primerDesign/final.pools.adaptors.txt')
pools = read.table(poolfile, header = F, stringsAsFactors = F)
colnames(pools) = c('pool','name','primer')

pools$plate = "Plate 1"
plateColumns = paste0(rep(LETTERS[1:8], 12), rep(1:12, each = 8))
pools$well = rep(plateColumns, each = 22)

## reorder columns and drop the pool names
pools = pools[, c('plate','well','name','primer')]
write.table(pools, paste0(dir,'/primerDesign/final.pools.formatted.txt'), quote = F, sep = '\t', row.names = F, col.names = F)

## get sites belonging to pool A, pool B
pools$sitename = sub('_[FR]', '', pools$name)
poolA = unique(pools$sitename[1:(nrow(pools)/2)])
poolB = unique(pools$sitename[(nrow(pools)/2 + 1):nrow(pools)])

## read in sites
sites = read.table(paste0(dir, '/selectedSites/genes.intersect.WGS79.OMNI86.selected.variants.ordered.firstDesign.txt'),
                   header = T, stringsAsFactors = F)
sites = rbind(sites,
              read.table(paste0(dir, '/selectedSites/genes.v2.intersect.WGS79.OMNI86.selected.variants.ordered.firstDesign.txt'),
                         header = T, stringsAsFactors = F))

sites = sites[sites$varname %in% c(poolA, poolB), ]
sites$pool = 'A'
sites$pool[sites$varname %in% poolB] = 'B'

write.table(sites[order(sites$varname), ],
            paste0(dir, '/selectedSites/genes.intersect.WGS79.OMNI86.selected.variants.primerpools.txt'),
            row.names = F, quote = F, sep = '\t')
