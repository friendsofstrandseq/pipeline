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

gz <- gzfile(snakemake@output[["counts"]], "w")
write.csv(counts, gz, row.names=F)
close(gz)


