log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library(data.table)
library(assertthat)
source("utils/simulate_counts.R")

sv <- fread(snakemake@input[["variants"]])
sce <- fread(snakemake@input[["sce"]])
info <- fread(snakemake@input[["info"]])
alpha <- snakemake@params[["alpha"]]
bin.size <- as.numeric(snakemake@wildcards[["window_size"]])
seed <- snakemake@wildcards[["seed"]]

# FIXME: fix this sample naming problem in the simulation
if (unique(sv[, sample]) != paste0("simulation", seed)) {
	print(paste0("WARNING!!!! sample name in SV file is not equal to simulation", seed))
	print("correcting the sample name")
	sv[, sample:=paste0("simulation", seed)]
	fwrite(sv, snakemake@input[["variants"]], sep="\t", quote=F, row.names=F)
}
if (unique(sce[, sample]) != paste0("simulation", seed)) {
	print(paste0("WARNING!!!! sample name in SCE file is not equal to simulation", seed))
	print("correcting the sample name")
	sce[, sample:=paste0("simulation", seed)]
	fwrite(sce, snakemake@input[["sce"]], sep="\t", quote=F, row.names=F)
}
if (unique(info[, sample]) != paste0("simulation", seed)) {
	print(paste0("WARNING!!!! sample name in info file is not equal to simulation", seed))
	print("correcting the sample name")
	info[, sample:=paste0("simulation", seed)]
	fwrite(info, snakemake@input[["info"]], sep="\t", quote=F, row.names=F)
}

counts <- simulateCounts(sv, sce, info, alpha, bin.size, seed)
setkey(counts, sample, cell, chrom, start, end)

# convert to non scientific mode
start_no_sci = gsub(" ", "", format(counts[, start], scientific=F))
end_no_sci = gsub(" ", "", format(counts[, end], scientific=F))
counts[, `:=`(start=NULL, end=NULL)]
counts[, `:=`(start=start_no_sci, end=end_no_sci)]

setcolorder(counts, c("chrom", "start", "end", "sample", "cell", "c", "w", "class"))

# set all CW classes to WC (because mosaicatcher does not accept them)
counts[class=="CW", class:="WC"]

gz <- gzfile(snakemake@output[["counts"]], "w")
write.table(counts, gz, sep="\t", quote=F, row.names=F)
close(gz)


