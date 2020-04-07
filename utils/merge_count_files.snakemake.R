log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library(data.table)

w <- snakemake@input[["w_counts"]]
c <- snakemake@input[["c_counts"]]

w.counts <- data.table()
c.counts <- data.table()
for (i in 1:length(w)){
	cell.name <- gsub("manual_segments_", "", basename(w[[i]]))
	cell.name <- gsub("_w_counts.txt", "", cell.name)
	cell.w.counts <- fread(w[[i]], col.names=c('chrom', 'start', 'end', 'W'))
	cell.w.counts[, cell:=cell.name]
	w.counts <- rbind(w.counts, cell.w.counts)
	cell.c.counts <- fread(c[[i]], col.names=c('chrom', 'start', 'end', 'C'))
	cell.c.counts[, cell:=cell.name]
	c.counts <- rbind(c.counts, cell.c.counts)
}

# TODO: add the class column later if needed
counts <- merge(c.counts, w.counts, by=c("chrom","start","end","cell"))
counts[, sample:=snakemake@wildcards[["sample"]]]

fwrite(counts[, .(chrom,start,end,sample,cell,C,W)], snakemake@output[[1]], sep="\t", row.names=F)
