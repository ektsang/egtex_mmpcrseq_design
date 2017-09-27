#!/usr/bin/env Rscript

# Author: Emily Tsang

# depict plate map to get a sense of sample layout

# shapes
# colors

library(ggplot2)

#############
# Functions #
#############

# function that adds a column with the row and the column
# assumes there is a column called Position with the plate position (e.g. A01)
add.row.column = function(df){
	pos = as.character(df$Position)
	df$row = factor(sapply(pos, substr, start=1, stop=1), levels=LETTERS[8:1])
	df$col = as.numeric(sapply(pos, substr, start=2, stop=3))
	return (df)
}

# function that takes care of plotting
# prints the plots instead of returning them
plot.plate = function(df, color.col, mapnum, colors=NULL) {
	
	for (plate in unique(df$Container)) {
		df.subset = df[df$Container==plate,]
		# tissue
		p = ggplot(df.subset, aes_string(x='col', y='row', color=color.col, size=2)) + geom_point()
		p = p + theme_bw() + xlab('') + ylab('')
		p = p + scale_x_discrete(breaks=seq(1:12), labels=c(1:12))
		p = p + guides(color=guide_legend(ncol=2)) + guides(size=F)
		p = p + ggtitle(paste("Shipment",mapnum,"Plate",plate))
		if (!is.null(colors)) {
			p = p + scale_colour_manual(values=colors)
		}
		print(p)
	}
}

########################################################################################

# read in data
map1 = read.table('~/projects/egtex_mmPCR/plate_maps/processed/eGTEx_RNA_CORE_set1_Jin_Billy_Li_Plate_Map.txt', sep="\t", header=T)

map2 = read.table('~/projects/egtex_mmPCR/plate_maps/processed/BillyLi_RNA_platemap_04_06_16.txt', sep='\t', header=T)

# subset to useful columns
map1 = map1[, c('Participant.ID','Container','Position','Position..ColWise.','Tissue.Site.Detail')]
map2 = map2[, c('Participant.ID.s.','Container','Position','Position..ColWise.','Tissue.Site.Detail')]
colnames(map2)[1] = 'Participant.ID'

# add numeric row and col columns
map1 = add.row.column(map1)
map2 = add.row.column(map2)

# actual plotting
pdf('~/projects/egtex_mmPCR/plate_maps/processed/plate.map.plots.bytissue.pdf', width=12, height=6, bg='white')
plot.plate(map1, 'Tissue.Site.Detail', '1')
plot.plate(map2, 'Tissue.Site.Detail', '2')
dev.off()

c16 = c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "palegreen2", "#FDBF6F", "gray70", "maroon", "orchid1", "darkturquoise", "darkorange4", "brown") 
pdf('~/projects/egtex_mmPCR/plate_maps/processed/plate.map.plots.byparticipant.pdf', width=8, height=5, bg='white')
plot.plate(map1, 'Participant.ID', '1', colors = c16)
plot.plate(map2, 'Participant.ID', '2')
dev.off()
