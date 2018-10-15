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

counts <- simulateCounts(sv, sce, info, alpha, bin.size, seed)
setkey(counts, sample, cell, chrom, start, end)

# convert to non scientific mode
start_no_sci = gsub(" ", "", format(counts[, start], scientific=F))
end_no_sci = gsub(" ", "", format(counts[, end], scientific=F))
counts[, `:=`(start=NULL, end=NULL)]
counts[, `:=`(start=start_no_sci, end=end_no_sci)]

setcolorder(counts, c("chrom", "start", "end", "sample", "cell", "c", "w", "class"))

gz <- gzfile(snakemake@output[["counts"]], "w")
write.table(counts, gz, sep="\t", quote=F, row.names=F)
close(gz)


